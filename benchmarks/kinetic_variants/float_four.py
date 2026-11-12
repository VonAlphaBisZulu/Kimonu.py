"""Float64 against exact, four ways: both routes, each with a float (SVD) kernel and with an exact kernel.

float_tradeoff.py asks how the tangent route fares in float64. This script separates the two places where
floating point needs a tolerance, and does so for both routes (detached uptakes throughout):

    kernel      SVD: the kernel comes from a dense SVD, numpy's default rank tolerance (max(shape) * eps
                relative to the largest singular value). Exact: the kernel is computed exactly and only then
                converted to an orthonormal float basis, equal exact rows to bitwise equal float rows.
    comparison  in both, rows count as equal when they differ by less than tau times the larger row norm,
                and as zero when their norm is below tau, for tau on a half-decade grid.

With an exact kernel, truly equal rows differ by exactly zero, so nothing can be missed, and every error is
a relation made up by the tolerance. With an SVD kernel, a wrong rank decision can also break relations.

The tangent route reads ACR and ACRR off T = ker(W Y^T) at r* (kimonu.interior_flux, as in tangent.py). The
iterated route alternates ker(W) (plus the feedback rows u_c - u_d) with L = ker((Y Delta)^T) as in
iterated.py, at the same r*; in floats every grouping in the loop (equal rows of
ker(W), equal rows of Y^T L) uses the same tau, so a wrong merge feeds back into the next round. Each arm is
scored against its own exact route: ACR metabolites and ACRR pairs, missed (false negatives) and made up
(false positives), over all metabolites and over the free ones.

Distances between nearly equal rows are taken directly, not from the Gram matrix: |a|^2 + |b|^2 - 2 a.b loses
everything below about 1e-8 relative, which would count truly equal rows as different at tight tolerances.

    python float_four.py <model> [part]   -> $KF_OUT/<model>.json, or $KF_OUT/parts/<model>__<part>.json

part: "tangent" (both tangent arms), or "iterated_exact:k" / "iterated_svd:k" for the k-th of four contiguous
chunks of the tolerance grid; without it, everything. figures/make_float_four.py merges the parts.
"""
import hashlib, json, os, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp, scipy.linalg as sl
from scipy.sparse.csgraph import minimum_spanning_tree
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.readout import complex_directions, free_directions, readout
from kimonu.io import zero_complexes
from kimonu.linalg import entries, exact, null_space_rows
from kimonu.steady_state import interior_flux, rate_columns
from run_kimonu_crn import triplets_to_sparse

OUT = os.environ.get("KF_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_float_four")
TAUS = [float(10.0 ** e) for e in np.arange(-14, -0.5 + 1e-9, 0.5)]
EXACT_TAU = 1e-200          # only exact equality passes: the check that the exact arms reproduce kimonu
REFINE = 1e-6               # pairs closer than this (from the Gram matrix) get their distance taken directly
EPS = np.finfo(float).eps


# ---------------------------------------------------------------- kernels
def svd(M, full):
    """numpy's SVD (LAPACK gesdd); where gesdd does not converge, which happens on some of these matrices,
    the slower QR-iteration driver gesvd. Recorded in FALLBACKS."""
    try:
        return np.linalg.svd(M, full_matrices=full)
    except np.linalg.LinAlgError:
        FALLBACKS.append(list(M.shape))
        return sl.svd(M, full_matrices=full, lapack_driver="gesvd", check_finite=False)


FALLBACKS = []


def svd_kernel(M):
    """Orthonormal basis (columns) of ker(M), numpy's default rank."""
    r, c = M.shape
    if r == 0:
        return np.eye(c)
    _, s, Vt = svd(M, True)
    rank = int((s > max(r, c) * EPS * s[0]).sum()) if s.size and s[0] > 0 else 0
    return Vt[rank:].T


def svd_left_kernel(M):
    """Orthonormal basis (columns) of ker(M^T) for M (m x c), numpy's default rank."""
    m, c = M.shape
    if c == 0:
        return np.eye(m)
    U, s, _ = svd(M, c < m)
    rank = int((s > max(m, c) * EPS * s[0]).sum()) if s.size and s[0] > 0 else 0
    return U[:, rank:]


def as_tuple(r):
    return r if isinstance(r, tuple) else tuple(sorted(r.items()))


class ExactBasis:
    """An exact basis in row form as an orthonormal float basis T = K R^-1 (R from the QR of K).

    Each distinct exact row is converted once, so equal exact rows give bitwise equal float rows, and
    ``apply`` maps further exact rows (such as rows of Y^T L) through the same R."""

    def __init__(self, rows):
        rows = [as_tuple(r) for r in rows]
        self.cols = sorted({j for r in rows for j, _ in r})
        self.idx = {j: k for k, j in enumerate(self.cols)}
        K = self._dense(rows)
        self.R = np.linalg.qr(K, mode="r") if self.cols else None
        self.T = self.apply(rows)

    def _dense(self, rows):
        K = np.zeros((len(rows), len(self.cols)))
        for i, r in enumerate(rows):
            for j, v in r:
                K[i, self.idx[j]] = float(v)
        return K

    def apply(self, rows):
        rows = [as_tuple(r) for r in rows]
        if not self.cols:
            return np.zeros((len(rows), 0))
        uniq, inv = {}, np.empty(len(rows), dtype=np.int64)
        for i, r in enumerate(rows):
            inv[i] = uniq.setdefault(r, len(uniq))
        U = self._dense(list(uniq))
        Tu = sl.solve_triangular(self.R, U.T, trans="T", lower=False, check_finite=False).T
        return Tu[inv]


# ---------------------------------------------------------------- row comparisons
def rel_dist(X):
    """||x_i - x_j|| / max(|x_i|, |x_j|) for all row pairs, and the row norms; close pairs taken directly."""
    nrm = np.linalg.norm(X, axis=1)
    D = X @ X.T
    D *= -2
    D += nrm[:, None] ** 2
    D += nrm[None, :] ** 2
    np.maximum(D, 0, out=D)
    np.sqrt(D, out=D)
    den = np.maximum(nrm[:, None], nrm[None, :])
    den[den == 0] = 1
    D /= den
    ii, jj = np.nonzero(np.triu(D < REFINE, 1))
    for a in range(0, len(ii), 20000):
        i, j = ii[a:a + 20000], jj[a:a + 20000]
        d = np.linalg.norm(X[i] - X[j], axis=1) / den[i, j]
        D[i, j] = d
        D[j, i] = d
    np.fill_diagonal(D, 0)
    return D, nrm


def mst(D):
    """Minimum spanning forest of the distance matrix over edges shorter than the largest tau; the
    single-linkage groups at any tau are its edges shorter than tau."""
    # built sparse: from a dense array scipy masks every weight below 1e-8 as "no edge"
    ii, jj = np.nonzero(np.triu(D < TAUS[-1], 1))
    W = sp.csr_matrix((np.maximum(D[ii, jj], 1e-300), (ii, jj)), shape=D.shape)
    T = minimum_spanning_tree(W).tocoo()
    return T.row, T.col, T.data


def groups(k, edges, tau, items):
    """Groups of items (single linkage at tau), as sets of item labels."""
    parent = list(range(k))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b, w in zip(*edges):
        if w < tau:
            ra, rb = find(int(a)), find(int(b))
            if ra != rb:
                parent[ra] = rb
    out = {}
    for a in range(k):
        out.setdefault(find(a), set()).add(items[a])
    return list(out.values())


def score(D, nrm, tau, acr_x, E, free, iu):
    """ACR and ACRR errors at tau, against the exact ACR set and the exact pair matrix E (upper triangle)."""
    tau = np.float64(tau)                               # a strong float64: D may be stored as float32
    acr_f = nrm < tau
    keep = ~acr_f
    found = (D[iu] < tau) & keep[iu[0]] & keep[iu[1]]
    Ex = E[iu]
    ff = free[iu[0]] & free[iu[1]]
    fp, fn = found & ~Ex, Ex & ~found
    return dict(acr_fp=int((acr_f & ~acr_x).sum()), acr_fn=int((acr_x & ~acr_f).sum()),
                acrr_fp=int(fp.sum()), acrr_fn=int(fn.sum()),
                acr_fp_free=int((acr_f & ~acr_x & free).sum()), acr_fn_free=int((acr_x & ~acr_f & free).sum()),
                acrr_fp_free=int((fp & ff).sum()), acrr_fn_free=int((fn & ff).sum()),
                h_acr=sha(np.flatnonzero(acr_f)), h_acrr=sha(np.flatnonzero(found)))


def sha(x):
    """A short fingerprint of a relation set, to compare the two routes' answers at one tolerance."""
    b = x.astype(np.int64).tobytes() if isinstance(x, np.ndarray) else repr(x).encode()
    return hashlib.sha1(b).hexdigest()[:16]


def hmods(sets):
    return sha(sorted(sorted(int(c) for c in s) for s in sets if len(s) > 1))


def exact_hashes(acr, E, modules):
    return dict(h_acr=sha(np.flatnonzero(acr)), h_acrr=sha(np.flatnonzero(E[iu])), h_mod=hmods(modules))


def pair_matrix(m, acrr_groups):
    E = np.zeros((m, m), bool)
    for g in acrr_groups:
        for a in g:
            for b in g:
                if a < b:
                    E[a, b] = True
    return E


# ---------------------------------------------------------------- the model
M = sys.argv[1]
PART = sys.argv[2] if len(sys.argv) > 2 else "all"
NCHUNK = 4
DO_TANGENT = PART in ("all", "tangent")
ITER = [("exact", None), ("svd", None)] if PART == "all" else \
       [(PART.split(":")[0].split("_")[1], int(PART.split(":")[1]))] if PART.startswith("iterated") else []
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
net = kimonu.crn_from_matrices(A, Y)                 # N = Y·A, the S of the crn.json
mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
free = np.zeros(nm, bool); free[:mets.index("E_1") if "E_1" in mets else nm] = True
skip = zero_complexes(net)
cx = [c for c in range(nc) if c not in skip]                 # the complexes that enter the modules


def w_entries(flux):
    """W of Eq. 3 as {(metabolite, column): value} and its number of columns, filled reaction by reaction as
    kimonu.steady_state_map did up to commit f948ec4, so the float sums below add up in the order of the recorded runs"""
    column, ncol = rate_columns(net, True)
    by_reaction = {}
    for i, j, v in entries(net.N)[0]:
        by_reaction.setdefault(j, []).append((i, v))
    W = {}
    for j, col in column.items():
        for i, v in by_reaction.get(j, ()):
            W[(i, col)] = W.get((i, col), Fraction(0)) + flux[j] * v
    return W, ncol


def yt_rows(B):
    """Rows of Yᵀ B for the complexes that enter the modules"""
    rows = complex_directions(net, B)
    return [rows[c] for c in cx]


Yf = sp.csc_matrix(Y, dtype=float)
iu = np.triu_indices(nm, 1)
out = dict(model=M, m=nm, n=nc, taus=TAUS)
log = lambda *a: print(f"[{time.strftime('%H:%M:%S')}] {M}:", *a, flush=True)


if DO_TANGENT:
    # ---------------------------------------------------------------- tangent route
    t0 = time.time()
    fluxT = interior_flux(net, 0)
    W, ncol = w_entries(fluxT)
    J, Jf = {}, np.zeros((nm, nm + ncol - nc))
    for (i, c), w in W.items():
        if c >= nc:
            J[(i, nm + c - nc)] = w; Jf[i, nm + c - nc] += float(w)
            continue
        for k, y in net.Yt[c]:
            J[(i, k)] = J.get((i, k), Fraction(0)) + w * y
            Jf[i, k] += float(w) * float(y)
    rows = null_space_rows(exact(((i, k, v) for (i, k), v in J.items()), (nm, nm + ncol - nc)))
    ex = readout(net, [as_tuple(r) for r in rows[:nm]])
    ref = kimonu.tangent_route(net, detach_uptakes=True)
    assert sorted(ex.acr) == sorted(ref.acr), "tangent: own readout differs from kimonu.tangent_route"
    acr_x = np.zeros(nm, bool); acr_x[list(ex.acr)] = True
    E = pair_matrix(nm, ex.acrr_groups)
    out["tangent_exact_relations"] = dict(acr=int(acr_x.sum()), acrr_pairs=int(E.sum()), acr_free=int((acr_x & free).sum()),
                                          acrr_pairs_free=int(E[np.ix_(free, free)].sum()), seconds=time.time() - t0)
    out["tangent_exact_relations"].update(exact_hashes(acr_x, E, ex.modules))
    log("tangent exact", out["tangent_exact_relations"])
    rows_m = [as_tuple(r) for r in rows[:nm]]


    def tangent_exact():
        eb = ExactBasis(rows)
        return eb.T[:nm], eb.apply(yt_rows(rows_m))


    def tangent_svd():
        K = svd_kernel(Jf)[:nm]
        return K, np.asarray(Yf.T @ K)[cx]


    for arm, make in (("tangent_exact", tangent_exact), ("tangent_svd", tangent_svd)):
        t0 = time.time()
        T, YtT = make()
        D, nrm = rel_dist(T)
        DY, _ = rel_dist(YtT); eY = mst(DY); del DY

        def at(t):
            r = score(D, nrm, t, acr_x, E, free, iu)
            r["h_mod"] = hmods(groups(len(cx), eY, t, cx))
            return r

        out[arm] = dict(sweep=[at(t) for t in TAUS], check=at(EXACT_TAU), seconds=time.time() - t0)
        del D, T, YtT
        log(arm, "done", round(out[arm]["seconds"], 1), "s; at exact equality:", out[arm]["check"])
    del Jf, rows

if ITER:
    # ---------------------------------------------------------------- iterated route
    t0 = time.time()
    km = kimonu.iterated_route(net, detach_uptakes=True, seed=0)
    fluxI = km.flux
    Rx = readout(net, free_directions(net, km.modules))
    acr_xi = np.zeros(nm, bool); acr_xi[list(Rx.acr)] = True
    Ei = pair_matrix(nm, Rx.acrr_groups)
    mods_x = sorted(sorted(s) for s in km.modules if len(s) > 1)
    out["iterated_exact_relations"] = dict(acr=int(acr_xi.sum()), acrr_pairs=int(Ei.sum()), acr_free=int((acr_xi & free).sum()),
                                           acrr_pairs_free=int(Ei[np.ix_(free, free)].sum()), modules=len(mods_x),
                                           rounds=km.rounds, seconds=time.time() - t0)
    out["iterated_exact_relations"].update(exact_hashes(acr_xi, Ei, km.modules))
    log("iterated exact", out["iterated_exact_relations"])
    Wx, ncolI = w_entries(fluxI); mW = nm
    Wf = np.zeros((mW, ncolI))
    for (i, c), w in Wx.items():
        Wf[i, c] = float(w)


    def make_arm(kind):
        """Memoized kernels and their comparison structures for one arm of the iterated route."""
        memo_W, memo_L = {}, {}
        K1 = svd_kernel(Wf) if kind == "svd" else None

        def kerW(equal):
            if equal not in memo_W:
                if kind == "svd":
                    if equal:
                        Ek = np.array([K1[c] - K1[d] for c, d in equal])
                        X = K1 @ svd_kernel(Ek)
                    else:
                        X = K1
                else:
                    Wq = dict(Wx)
                    for k, (c, d) in enumerate(equal):
                        Wq[(mW + k, c)], Wq[(mW + k, d)] = Fraction(1), Fraction(-1)
                    X = ExactBasis(null_space_rows(exact(((i, j, v) for (i, j), v in Wq.items()),
                                                         (mW + len(equal), ncolI)))).T
                D, _ = rel_dist(X[cx])
                memo_W[equal] = mst(D)
                del D
            return memo_W[equal]

        def kerL(sets):
            key = frozenset(frozenset(s) for s in sets if len(s) > 1)
            if key not in memo_L:
                if kind == "svd":
                    cols = []
                    for s in key:
                        ms = sorted(s)
                        for o in ms[1:]:
                            cols.append((Yf[:, ms[0]] - Yf[:, o]).toarray().ravel())
                    L = svd_left_kernel(np.array(cols).T if cols else np.zeros((nm, 0)))
                    YtL = np.asarray(Yf.T @ L)[cx]
                else:
                    Lrows = free_directions(net, [set(s) for s in key])
                    eb = ExactBasis(Lrows)
                    L = eb.T
                    YtL = eb.apply(yt_rows(Lrows))
                DL, nrmL = rel_dist(L)
                DY, _ = rel_dist(YtL)
                memo_L[key] = (DL.astype(np.float32), nrmL, mst(DY))
                del DL, DY
            return memo_L[key]

        def run(tau):
            equal, rnd = (), 0
            while True:
                rnd += 1
                coupled = groups(len(cx), kerW(equal), tau, cx)
                DL, nrmL, eY = kerL(coupled)
                merged = groups(len(cx), eY, tau, cx)
                if len(merged) == len(coupled) or rnd >= 20:
                    break
                equal = tuple(sorted((min(s), c) for s in merged for c in s if c != min(s)))
            r = score(DL, nrmL, tau, acr_xi, Ei, free, iu)
            r.update(rounds=rnd, modules=sum(len(s) > 1 for s in merged), h_mod=hmods(merged))
            return r, merged

        return run


    for kind, chunk in ITER:
        arm = f"iterated_{kind}"
        t0 = time.time()
        run = make_arm(kind)
        taus = TAUS if chunk is None else [t for k, t in enumerate(TAUS) if k * NCHUNK // len(TAUS) == chunk]
        out[arm] = dict(sweep=[], taus=taus)
        if chunk in (None, 0):
            chk, merged = run(EXACT_TAU)
            chk["modules_match"] = sorted(sorted(s) for s in merged if len(s) > 1) == mods_x
            out[arm]["check"] = chk
            log(arm, "at exact equality:", chk)
        for t in taus:
            out[arm]["sweep"].append(run(t)[0])
            log(arm, f"tau {t:.1e}", out[arm]["sweep"][-1])
        out[arm]["seconds"] = time.time() - t0

out["svd_fallbacks"] = FALLBACKS
dest = f"{OUT}/{M}.json" if PART == "all" else f"{OUT}/parts/{M}__{PART.replace(':', '_')}.json"
os.makedirs(os.path.dirname(dest), exist_ok=True)
json.dump(out, open(dest, "w"))
log("wrote", dest)
