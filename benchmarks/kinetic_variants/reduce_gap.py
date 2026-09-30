"""Shrink a model to a minimal network that still shows a tangent-only relation (delta debugging).

Removal units are enzyme mechanisms: every reaction that involves one enzyme (E_k or a bound form Ek_...)
forms one unit, and every other reaction is its own unit. Removing a mechanism removes one enzyme's whole
contribution; nothing is lumped, so the kinetics of the kept reactions are unchanged.

The predicate: after removing a set of reactions (and whatever becomes blocked), the tangent space still
reports metabolite TARGET as fixed and the iterated route still does not. Reactions are removed in chunks
of decreasing size (ddmin); every accepted removal is logged and the current network saved, so a run can
be stopped at any time.

Units are tried far-to-near: ordered by their distance from the target in the metabolite-reaction
graph (the 20 most connected metabolites left out, so cofactors do not make everything adjacent), and
removed as contiguous blocks, so large distant parts go first. All checks are exact: a floating-point
screen (least squares for e_target in the row space of J) called the resumed iHN637 network "not fixed"
where the exact test says fixed.

    python kinetic_variants/reduce_gap.py <model> <target metabolite> <out.json> [resume.json]
"""
import json, os, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import kimonu
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.tangent_one_shot import one_shot_exact

M, target, out = sys.argv[1:4]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A0 = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64).tocsc(); Y0 = triplets_to_sparse(crn["Y"], nm, nc).tocsc()
N0 = triplets_to_sparse(crn["S"], nm, nr).tocsc(); mets = crn["metabolites"]; t_idx = mets.index(target)

def unblocked(keep):
    """support of the flux cone of the kept reactions: one LP, max sum t with t <= v, t <= 1"""
    Nk = N0[:, keep]; q = len(keep)
    lp = linprog(np.concatenate([np.zeros(q), -np.ones(q)]),
                 A_ub=sp.hstack([-sp.identity(q), sp.identity(q)]).tocsr(), b_ub=np.zeros(q),
                 A_eq=sp.hstack([Nk, sp.csr_matrix((nm, q))]).tocsr(), b_eq=np.zeros(nm),
                 bounds=[(0, None)] * q + [(0, 1)] * q, method="highs")
    return [keep[j] for j in range(q) if lp.status == 0 and lp.x[q + j] > 0.5]

def network(keep):
    A = A0[:, keep]; used_c = np.flatnonzero(np.abs(A).sum(axis=1).A.ravel())
    A = A[used_c, :]; Y = Y0[:, used_c]; used_m = np.flatnonzero(np.abs(Y).sum(axis=1).A.ravel())
    return A.tocsr().astype(np.int64), Y[used_m, :].tocsr(), to_exact(N0[used_m, :][:, keep]), [mets[i] for i in used_m]

from scipy.sparse.linalg import lsqr
from kimonu.kinetic import positive_flux

def tangent_fixed_float(A, Y, N, ti):
    """e_ti in the row space of J = N diag(r*) G, in floating point"""
    r = np.array([float(v) for v in positive_flux(N)])
    Nf = N0_float(N); G = ((A < 0).astype(float).T @ Y.T).tocsr()
    J = (Nf @ sp.diags(r) @ G).tocsr()
    e = np.zeros(J.shape[1]); e[ti] = 1
    w = lsqr(J.T.tocsr(), e, atol=1e-14, btol=1e-14, iter_lim=20000)[0]
    return np.linalg.norm(J.T @ w - e) < 1e-6

def N0_float(N):
    coo = N.to_coo_exact(); return sp.csr_matrix((np.array(coo.data, float) / float(coo.denom), (coo.rows, coo.cols)), shape=coo.shape)

cache = {}
def gap(keep, exact=False):
    key = frozenset(keep)
    if (key, exact) in cache: return cache[(key, exact)]
    ok = False
    try:
        kept = unblocked(sorted(keep))
        A, Y, N, ms = network(kept)
        if target in ms:
            ti = ms.index(target)
            it = kimonu.kinetic_modules(A, Y, N)
            acr_it = set(kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr).acr)
            if ti not in acr_it:
                if exact:
                    _, acr_t, _, _ = one_shot_exact(A, Y, N, False)
                    ok = ti in set(acr_t)
                else:
                    ok = tangent_fixed_float(A, Y, N, ti)
    except Exception as e:
        print("   predicate error:", type(e).__name__, str(e)[:80], flush=True)
    cache[(key, exact)] = ok if exact else cache.get((key, exact), ok)
    return ok

def save(keep, note):
    kept = unblocked(sorted(keep)); A, Y, N, ms = network(kept)
    Acsc, Ycsc = A.tocsc(), Y.tocsc()
    def cname(c): return "+".join(f"{'' if v == 1 else v}{ms[i]}" for i, v in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]])) or "0"
    rx = []
    for j in range(Acsc.shape[1]):
        col = Acsc[:, j].toarray().ravel(); rx.append([cname(int(np.flatnonzero(col < 0)[0])), cname(int(np.flatnonzero(col > 0)[0]))])
    json.dump(dict(model=M, target=target, note=note, reactions=kept, readable=rx, n_metabolites=len(ms)), open(out, "w"), indent=1)
    return len(kept), len(ms)

keep = unblocked(list(range(nr)))
if len(sys.argv) > 4:
    keep = unblocked(json.load(open(sys.argv[4]))["reactions"])
import re
Ncsc = N0.tocsc(); by_enzyme = {}
for j in keep:
    names = [mets[i] for i in Ncsc.indices[Ncsc.indptr[j]:Ncsc.indptr[j + 1]]]
    ids = {m.group(1) for n in names for m in [re.match(r"E_?(\d+)(?:_|$)", n)] if m}
    by_enzyme.setdefault(min(ids, key=int) if ids else f"r{j}", []).append(j)
units = list(by_enzyme.values())
# distance of each unit from the target, hubs excluded
Nr = N0.tocsr(); deg = np.diff(Nr.indptr); hubs = set(np.argsort(-deg)[:20]) - {t_idx}
dist = {t_idx: 0}; frontier = [t_idx]; rx_d = {}
while frontier:
    nxt = []
    for i in frontier:
        for j in Nr.indices[Nr.indptr[i]:Nr.indptr[i + 1]]:
            if j in rx_d: continue
            rx_d[j] = dist[i]
            for i2 in Ncsc.indices[Ncsc.indptr[j]:Ncsc.indptr[j + 1]]:
                if i2 not in dist and i2 not in hubs: dist[i2] = dist[i] + 1; nxt.append(i2)
    frontier = nxt
units.sort(key=lambda u: -min(rx_d.get(j, 10**6) for j in u))
t0 = time.time()
print(f"{len(keep)} reactions in {len(units)} units ({sum(1 for k in by_enzyme if not k.startswith('r'))} enzyme mechanisms)", flush=True)
print(f"{M} {target}: {len(keep)} unblocked reactions; gap on full model: {gap(keep, exact=True)} ({time.time() - t0:.0f}s)", flush=True)
n = 2
while len(units) >= 2:
    size = -(-len(units) // n); chunks = [units[i:i + size] for i in range(0, len(units), size)]   # far first
    removed = False
    for ch in chunks:
        drop = {r for u in ch for r in u}
        trial = [r for r in keep if r not in drop]
        if trial and gap(trial, exact=True):
            keep = unblocked(sorted(trial)); kept = set(keep)
            units = [u for u in units if kept.issuperset(u)]
            n = max(n - 1, 2); removed = True
            r_, m_ = save(keep, "in progress")
            print(f"[{time.time() - t0:7.0f}s] removed {len(ch)} groups: {r_} reactions in {len(units)} groups, {m_} metabolites left", flush=True)
            break
    if not removed:
        if n >= len(units): break
        n = min(len(units), 2 * n)
        print(f"[{time.time() - t0:7.0f}s] no chunk removable, splitting finer: n = {n} of {len(units)} groups", flush=True)
r_, m_ = save(keep, "1-minimal: removing any single reaction loses the gap")
print(f"done: {r_} reactions, {m_} metabolites ({time.time() - t0:.0f}s)", flush=True)
