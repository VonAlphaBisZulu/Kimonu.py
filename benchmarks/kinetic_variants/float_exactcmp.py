"""Exact kernel, tolerance comparison in exact arithmetic: the third kernel arm of float_four.py.

Both routes compute their kernels exactly, as Kimonu does, in Kimonu's own exact (RREF) bases, and keep every
entry as a fraction. Only the comparison of rows takes a tolerance tau, absolute and per element, as Kimonu's
readout would with "equal" replaced by "within tau":

    ker(W)            two complexes coupled when max_k |r_ik - r_jk| <= tau (no scaling: u = 1 lies in ker(W), so
                      coupled complexes have equal rows, not only proportional ones)
    ACR               max_k |r_ik| <= tau
    ACRR, modules     max_k |r_ik - r_jk| <= tau (rows of B, and of Y^T B), no scaling

Every decision is exact: the entries are fractions, and |x| <= tau is tested as x^2 <= tau^2, which is rational
on the half-decade grid. Floats only propose candidate pairs (a sort along a random direction w with |w|_1 = 1,
since |w.(a - b)| <= max_k |a_k - b_k|); a pair is decided in fractions whenever its float distance is near any
threshold. Truly equal rows have distance exactly zero, so this arm can only make relations up, never miss one;
at tau = 0 it is Kimonu. Groups are single-linkage, ACRR pairs pairwise, as in float_four.py. Detached uptakes.

The same element-wise rule is also applied to the floating-point (SVD) kernels of float_four.py (parts with the
prefix svd_), in the orthonormal bases the SVD returns, so that both kernels are compared by one rule; only the
basis differs (RREF for the exact kernel, orthonormal for the SVD).

    python float_exactcmp.py <model> <part>    part: "tangent", "iterated:k" (k-th of four tolerance chunks) or "iterated:tN" (tolerance N alone),
                                                "svd_tangent" or "svd_iterated:k"
                                                -> $KF_OUT/parts/<model>__<arm>.json, merged by make_float_four.py
"""
import hashlib, json, os, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from concurrent.futures import ThreadPoolExecutor
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.acr import _combined, _difference_columns, _kernel_rows, free_directions, readout
from kimonu.io import zero_complexes
from kimonu.kinetic import steady_state_map
from kimonu.linalg import exact, null_space_rows, positive_flux
from run_kimonu_crn import triplets_to_sparse, to_exact

OUT = os.environ.get("KF_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_float_four")
TAUS = [float(10.0 ** e) for e in np.arange(-14, -0.5 + 1e-9, 0.5)]
TAU2 = [Fraction(10) ** int(round(2 * e)) for e in np.arange(-14, -0.5 + 1e-9, 0.5)]     # tau^2, exact
EXACT_TAU, NCHUNK, MEMO = 1e-200, 4, int(os.environ.get("KF_MEMO", 8))
M, PART = sys.argv[1], sys.argv[2]
SVD = PART.startswith("svd_"); PART = PART[4:] if SVD else PART
ARM = "svdabs" if SVD else "exactcmp"
EPS = np.finfo(float).eps
THREADS = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))


def svd(Mx, full):
    try:
        return np.linalg.svd(Mx, full_matrices=full)
    except np.linalg.LinAlgError:
        import scipy.linalg as sl
        return sl.svd(Mx, full_matrices=full, lapack_driver="gesvd", check_finite=False)


def svd_kernel(Mx):
    r, c = Mx.shape
    if r == 0: return np.eye(c)
    _, s_, Vt = svd(Mx, True)
    rank = int((s_ > max(r, c) * EPS * s_[0]).sum()) if s_.size and s_[0] > 0 else 0
    return Vt[rank:].T


def svd_left_kernel(Mx):
    m_, c = Mx.shape
    if c == 0: return np.eye(m_)
    U, s_, _ = svd(Mx, c < m_)
    rank = int((s_ > max(m_, c) * EPS * s_[0]).sum()) if s_.size and s_[0] > 0 else 0
    return U[:, rank:]
log = lambda *a: print(f"[{time.strftime('%H:%M:%S')}] {M} {PART}:", *a, flush=True)


def as_tuple(r):
    return r if isinstance(r, tuple) else tuple(sorted(r.items()))


class Rows:
    """Exact rows (tuples of (column, Fraction)) with exact tolerance tests. Equal rows are merged first."""

    def __init__(self, rows):
        if isinstance(rows, np.ndarray):                       # floating point (SVD): no exact refinement
            self.exact = False
            keys = [r.tobytes() for r in rows]
            self.uniq, self.inv = {}, np.empty(len(rows), dtype=np.int64)
            first = []
            for i, k in enumerate(keys):
                if k not in self.uniq: self.uniq[k] = len(self.uniq); first.append(i)
                self.inv[i] = self.uniq[k]
            self.F = rows[first] if len(first) else rows[:0]
            self.uniq = None
            self.u = list(range(len(first)))
            self.scale = np.abs(self.F).max(axis=1) if self.F.shape[1] else np.zeros(len(first))
            return
        self.exact = True
        rows = [as_tuple(r) for r in rows]
        self.uniq, self.inv = {}, np.empty(len(rows), dtype=np.int64)
        for i, r in enumerate(rows):
            self.inv[i] = self.uniq.setdefault(r, len(self.uniq))
        self.u = list(self.uniq)
        cols = sorted({j for r in self.u for j, _ in r}); idx = {j: k for k, j in enumerate(cols)}
        ii = [i for i, r in enumerate(self.u) for _ in r]
        jj = [idx[j] for r in self.u for j, _ in r]
        vv = [float(v) for r in self.u for _, v in r]
        self.F = sp.csr_matrix((vv, (ii, jj)), shape=(len(self.u), len(cols)))   # exact rows are very sparse
        self.nnz = np.array([len(r) for r in self.u])
        self.scale = np.array([max((abs(float(v)) for _, v in r), default=0.0) for r in self.u])

    # ---- one row: max |entry|, decided exactly near thresholds
    def norm_levels(self):
        """for every unique row, the number of tau (ascending) with max|r| > tau, i.e. its row is 'zero' for the
        taus from that index on; exact."""
        if not self.exact:                                     # -1: exactly zero; else count of tau < max|r|
            return np.where(self.scale == 0, -1, np.searchsorted(TAUS, self.scale, side="left"))
        out = np.empty(len(self.u), dtype=np.int64)
        for i, r in enumerate(self.u):
            if not r:
                out[i] = -1; continue
            mx2 = max(v * v for _, v in r)
            out[i] = sum(1 for t2 in TAU2 if mx2 > t2)
        return out                                     # row i is zero at tau_k exactly when k >= out[i]

    # ---- pairs: max |a - b|, decided exactly
    def pairs(self, tau_max):
        """int32 array of (i, j, level) for unique-row pairs with max|a - b| <= tau_max; level = first tau index at
        which they are within tau (they are within tau_k for k >= level)."""
        n = len(self.u)
        none = np.empty((0, 3), dtype=np.int32)
        if n < 2: return none
        if self.exact:
            w = np.random.default_rng(0).standard_normal(self.F.shape[1]); w /= np.abs(w).sum()
            proj = self.F @ w
        else:                                                  # dense rows: sort along the column with the widest spread,
            spread = self.F.std(axis=0)                        # whose window is tau itself; a dense random direction with
            top = np.argsort(spread)[::-1]                     # |w|_1 = 1 would shrink distances by about sqrt(dim)
            proj = self.F[:, top[0]]
            probe = self.F[:, top[1:65]]                       # the next widest columns reject most candidates cheaply
        order = np.argsort(proj); ps = proj[order]
        T = tau_max * (1 + 1e-9)
        # window per row. |w.(a - b)| <= max_k |a_k - b_k| <= tau, so a pair within tau has projections within
        # tau + eta_a + eta_b, eta being the rounding error of a projection, and hence within the wider of the two
        # windows tau + 2 eta. eta follows each row's own size: a few rows have entries near 1e15, and a global
        # margin would make every pair a candidate.
        if self.exact:
            W = T + 2 * 8 * EPS * (self.nnz[order] + 2) * self.scale[order]
        else:
            W = np.full(n, T + 1e-12 * (np.abs(ps).max() + 1))
        ends = np.searchsorted(ps, ps + W, side="right")
        tmax2 = max(t2 for t2, t in zip(TAU2, TAUS) if t <= T)
        taus = np.array(TAUS); out = [none]; step = 200000 if self.exact else 5000

        def float_batch(c):                                    # NumPy releases the GIL here, so these run side by side
            c = c[np.abs(probe[c[:, 0]] - probe[c[:, 1]]).max(axis=1) <= W[0]] if probe.shape[1] else c
            if not len(c):
                return none
            d = np.abs(self.F[c[:, 0]] - self.F[c[:, 1]]).max(axis=1)
            err = 8e-16 * np.maximum(self.scale[c[:, 0]], self.scale[c[:, 1]])
            ok = d - err <= T
            return np.column_stack([c[ok], np.searchsorted(taus, d[ok], side="left")]).astype(np.int32)

        def decide(p0, p1):                                    # sorted positions of candidate pairs
            if not self.exact:
                batches = (np.stack([order[p0[s1:s1 + step]], order[p1[s1:s1 + step]]], 1) for s1 in range(0, len(p0), step))
                if THREADS > 1:
                    with ThreadPoolExecutor(THREADS) as pool:
                        out.extend(pool.map(float_batch, batches))
                else:
                    out.extend(map(float_batch, batches))
                return
            for s1 in range(0, len(p0), step):
                c = np.stack([order[p0[s1:s1 + step]], order[p1[s1:s1 + step]]], 1)
                d = abs(self.F[c[:, 0]] - self.F[c[:, 1]]).max(axis=1).toarray().ravel()
                err = 8e-16 * np.maximum(self.scale[c[:, 0]], self.scale[c[:, 1]])
                ok = d - err <= T
                c, d, err = c[ok], d[ok], err[ok]
                lev = np.searchsorted(taus, d, side="left")    # first tau >= d (float); decided in fractions unless clearly between two tau
                hi, lo = taus[np.minimum(lev, len(taus) - 1)], taus[np.maximum(lev - 1, 0)]
                gap = np.minimum(np.abs(d - hi), np.abs(d - lo))
                keep = (gap > 2 * err) & (gap > 1e-6 * lo)
                for i, j in c[~keep]:
                    ra, rb = dict(self.u[i]), dict(self.u[j])
                    d2 = max(((ra.get(k, 0) - rb.get(k, 0)) ** 2 for k in set(ra) | set(rb)), default=Fraction(0))
                    if d2 <= tmax2:
                        out.append(np.array([[i, j, sum(1 for t2 in TAU2 if d2 > t2)]], dtype=np.int32))
                c, lev = c[keep], lev[keep]
                out.append(np.column_stack([c, lev]).astype(np.int32))

        counts = ends - np.arange(n) - 1                       # forward: the window of the lower row
        starts = np.flatnonzero(counts > 0)
        for s0 in range(0, len(starts), 256):
            sa = starts[s0:s0 + 256]
            decide(np.repeat(sa, counts[sa]), np.concatenate([np.arange(a + 1, ends[a]) for a in sa]))
        if self.exact:                                         # backward: the window of the upper row, where wider
            begin = np.searchsorted(ps, ps - W, side="left")
            for a in np.flatnonzero(begin < np.arange(n)):
                b = np.arange(begin[a], a)
                b = b[ends[b] <= a]
                if len(b): decide(b, np.full(len(b), a))
        return np.concatenate(out)


def groups_at(n_items, pairs, k, items):
    """single-linkage groups of items at tau index k (k = -1: exact equality only)"""
    sel = pairs[pairs[:, 2] <= k] if k >= 0 else pairs[:0]
    graph = sp.coo_matrix((np.ones(len(sel), dtype=bool), (sel[:, 0], sel[:, 1])), shape=(n_items, n_items))
    labels = connected_components(graph, directed=False)[1]
    out = {}
    for a, l_ in enumerate(labels):
        out.setdefault(int(l_), set()).add(items[a])
    return list(out.values())


def sha(x):
    b = x.astype(np.int64).tobytes() if isinstance(x, np.ndarray) else repr(x).encode()
    return hashlib.sha1(b).hexdigest()[:16]


hmods = lambda sets: sha(sorted(sorted(int(c) for c in s) for s in sets if len(s) > 1))


def score(R, lev_norm, pairs, k, acr_x, E, free, iu):
    """ACR and ACRR errors at tau index k (k = -1: exact equality), as float_four.score; a true ACRR pair whose
    members are both made ACR is counted as made up, not as missed."""
    acr_u = (lev_norm <= k) if k >= 0 else (lev_norm == -1)
    acr_f = acr_u[R.inv]
    sel = pairs[pairs[:, 2] <= k] if k >= 0 else pairs[:0]
    U = np.eye(len(R.u), dtype=bool)                        # equal rows have distance 0
    U[sel[:, 0], sel[:, 1]] = True; U[sel[:, 1], sel[:, 0]] = True
    F = U[R.inv[:, None], R.inv]
    keep = ~acr_f
    found = F[iu] & keep[iu[0]] & keep[iu[1]]
    Ex = E[iu]; both = acr_f[iu[0]] & acr_f[iu[1]]
    ff = free[iu[0]] & free[iu[1]]
    fp, fn = found & ~Ex, Ex & ~found & ~both
    return dict(acr_fp=int((acr_f & ~acr_x).sum()), acr_fn=int((acr_x & ~acr_f).sum()),
                acrr_fp=int(fp.sum()), acrr_fn=int(fn.sum()),
                acr_fp_free=int((acr_f & ~acr_x & free).sum()), acr_fn_free=int((acr_x & ~acr_f & free).sum()),
                acrr_fp_free=int((fp & ff).sum()), acrr_fn_free=int((fn & ff).sum()),
                h_acr=sha(np.flatnonzero(acr_f)), h_acrr=sha(np.flatnonzero(found)))


def pair_matrix(m, groups):
    E = np.zeros((m, m), bool)
    for g in groups:
        for a in g:
            for b in g:
                if a < b: E[a, b] = True
    return E


crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
mets = [m_ if isinstance(m_, str) else m_.get("id") for m_ in crn["metabolites"]]
free = np.zeros(nm, bool); free[:mets.index("E_1") if "E_1" in mets else nm] = True
skip = zero_complexes(Y); cx = [c for c in range(nc) if c not in skip]
Ycsc = sp.csc_matrix(Y); Yf = sp.csc_matrix(Y, dtype=float); iu = np.triu_indices(nm, 1)
out = dict(model=M, m=nm, n=nc, taus=TAUS)
tau_max = TAUS[-1]


def mods_from(Rc, pc, k):
    """modules: complexes (cx order) whose unique rows are grouped at tau index k"""
    gu = groups_at(len(Rc.u), pc, k, list(range(len(Rc.u))))
    of = {}
    for gi, g in enumerate(gu):
        for u in g: of[u] = gi
    mods = {}
    for a, c in enumerate(cx):
        mods.setdefault(of[int(Rc.inv[a])], set()).add(c)
    return list(mods.values())


if PART == "tangent":
    t0 = time.time()
    flux = positive_flux(N, 0, scale=10**6)
    W, (_, ncol) = steady_state_map(A, Y, N, flux, True)
    J = {}
    for (i, c), w in W.items():
        if c >= nc:
            J[(i, nm + c - nc)] = w; continue
        for k_, y in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
            J[(i, int(k_))] = J.get((i, int(k_)), Fraction(0)) + w * Fraction(str(y))
    rows = null_space_rows(exact(((i, k_, v) for (i, k_), v in J.items()), (nm, nm + ncol - nc)))
    rows_m = [as_tuple(r) for r in rows[:nm]]
    ex = readout(Y, rows_m, skip)
    acr_x = np.zeros(nm, bool); acr_x[list(ex.acr)] = True
    E = pair_matrix(nm, ex.acrr)
    out["tangent_exact_relations"] = dict(acr=int(acr_x.sum()), acrr_pairs=int(E.sum()),
                                          h_acr=sha(np.flatnonzero(acr_x)), h_acrr=sha(np.flatnonzero(E[iu])),
                                          h_mod=hmods(ex.modules))
    if SVD:
        Jf = np.zeros((nm, nm + ncol - nc))
        for (i, c), w in W.items():
            if c >= nc: Jf[i, nm + c - nc] += float(w); continue
            for k_, y in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
                Jf[i, int(k_)] += float(w) * float(y)
        Tm = np.ascontiguousarray(svd_kernel(Jf)[:nm]); del Jf
        Rm = Rows(Tm); Rc = Rows(np.ascontiguousarray(np.asarray(Yf.T @ Tm)[cx]))
    else:
        Rm = Rows(rows_m); Rc = Rows([_combined(Ycsc, rows_m, c) for c in cx])
    lev_m = Rm.norm_levels(); pm = Rm.pairs(tau_max); pc = Rc.pairs(tau_max)
    at = lambda k: dict(score(Rm, lev_m, pm, k, acr_x, E, free, iu), h_mod=hmods(mods_from(Rc, pc, k)))
    out[f"tangent_{ARM}"] = dict(sweep=[at(k) for k in range(len(TAUS))], check=at(-1), seconds=time.time() - t0)
    log("check", {k: out[f"tangent_{ARM}"]["check"][k] for k in ("acr_fp", "acr_fn", "acrr_fp", "acrr_fn")},
        "modules right:", out[f"tangent_{ARM}"]["check"]["h_mod"] == out["tangent_exact_relations"]["h_mod"])
else:
    which = PART.split(":")[1]
    if which.startswith("t"):                                  # one tolerance alone, for the slowest network
        ks = [int(which[1:])]
    else:
        chunk = int(which)
        ks = [k for k in range(len(TAUS)) if k * NCHUNK // len(TAUS) == chunk] + ([-1] if chunk == 0 else [])
    tau_max = max(TAUS[k] for k in ks if k >= 0)      # the pairs for larger tolerances are never used in this chunk
    t0 = time.time()
    km = kimonu.kinetic_modules(A, Y, N, detach_uptakes=True, seed=0)
    Rx = readout(Y, free_directions(Y, km.modules), skip)
    acr_x = np.zeros(nm, bool); acr_x[list(Rx.acr)] = True
    E = pair_matrix(nm, Rx.acrr)
    out["iterated_exact_relations"] = dict(acr=int(acr_x.sum()), acrr_pairs=int(E.sum()),
                                           h_acr=sha(np.flatnonzero(acr_x)), h_acrr=sha(np.flatnonzero(E[iu])),
                                           h_mod=hmods(km.modules))
    Wx, (mW, ncolI) = steady_state_map(A, Y, N, km.flux, True)
    memo_W, memo_L = {}, {}
    if SVD:
        Wf = np.zeros((mW, ncolI))
        for (i, c), w in Wx.items(): Wf[i, c] = float(w)
        K1 = svd_kernel(Wf); del Wf

    def kerW(equal):
        if equal not in memo_W and SVD:
            if len(memo_W) >= MEMO: memo_W.clear()     # keep memory bounded; only recomputation is lost
            X = K1 @ svd_kernel(np.array([K1[c] - K1[d] for c, d in equal])) if equal else K1
            R = Rows(np.ascontiguousarray(X[cx]))
            memo_W[equal] = (R, R.pairs(tau_max))
        if equal not in memo_W:
            if len(memo_W) >= MEMO: memo_W.clear()
            Wq = dict(Wx)
            for k_, (c, d) in enumerate(equal):
                Wq[(mW + k_, c)], Wq[(mW + k_, d)] = Fraction(1), Fraction(-1)
            rows = null_space_rows(exact(((i, j, v) for (i, j), v in Wq.items()), (mW + len(equal), ncolI)))
            R = Rows([as_tuple(rows[c]) for c in cx])
            memo_W[equal] = (R, R.pairs(tau_max))
        return memo_W[equal]

    def kerL(sets):
        key = frozenset(frozenset(s) for s in sets if len(s) > 1)
        if key not in memo_L and SVD:
            if len(memo_L) >= MEMO: memo_L.clear()
            cols = []
            for ms in sorted(map(sorted, key)):                # one order, so the basis does not depend on the path
                cols += [(Yf[:, ms[0]] - Yf[:, o]).toarray().ravel() for o in ms[1:]]
            L = svd_left_kernel(np.array(cols).T if cols else np.zeros((nm, 0)))
            Rm = Rows(np.ascontiguousarray(L)); Rc = Rows(np.ascontiguousarray(np.asarray(Yf.T @ L)[cx]))
            memo_L[key] = (Rm, Rm.norm_levels(), Rm.pairs(tau_max), Rc, Rc.pairs(tau_max))
        if key not in memo_L:
            if len(memo_L) >= MEMO: memo_L.clear()
            Lrows = _kernel_rows(_difference_columns(Y, [set(s) for s in sorted(map(sorted, key))], ()), nm)
            Rm = Rows(Lrows); lev = Rm.norm_levels()
            Rc = Rows([_combined(Ycsc, Lrows, c) for c in cx])
            memo_L[key] = (Rm, lev, Rm.pairs(tau_max), Rc, Rc.pairs(tau_max))
        return memo_L[key]

    def run(k):
        equal, rnd = (), 0
        while True:
            rnd += 1
            Rw, pw = kerW(equal)
            coupled = mods_from(Rw, pw, k)
            Rm, lev, pm, Rc, pc = kerL(coupled)
            merged = mods_from(Rc, pc, k)
            if len(merged) == len(coupled) or rnd >= 20: break
            equal = tuple(sorted((min(s), c) for s in merged for c in s if c != min(s)))
        r = score(Rm, lev, pm, k, acr_x, E, free, iu)
        r.update(rounds=rnd, modules=sum(len(s) > 1 for s in merged), h_mod=hmods(merged))
        return r

    arm = dict(sweep=[], taus=[TAUS[k] for k in ks if k >= 0])
    for k in ks:
        r = run(k)
        if k < 0:
            arm["check"] = r
            log("check", {x: r[x] for x in ("acr_fp", "acr_fn", "acrr_fp", "acrr_fn")},
                "modules right:", r["h_mod"] == out["iterated_exact_relations"]["h_mod"])
        else:
            arm["sweep"].append(r); log(f"tau {TAUS[k]:.1e}", {x: r[x] for x in ("acr_fp", "acr_fn", "acrr_fp", "acrr_fn", "rounds")})
    arm["seconds"] = time.time() - t0
    out[f"iterated_{ARM}"] = arm

dest = f"{OUT}/parts/{M}__{f'tangent_{ARM}' if PART == 'tangent' else f'iterated_{ARM}_' + PART.split(':')[1]}.json"
os.makedirs(os.path.dirname(dest), exist_ok=True)
json.dump(out, open(dest, "w"))
log("wrote", dest, round(time.time() - t0, 1), "s")
