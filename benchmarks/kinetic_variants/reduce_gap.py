"""Shrink a model to a minimal network that still shows a tangent-only relation (delta debugging).

The predicate: after removing a set of reactions (and whatever becomes blocked), the tangent space still
reports metabolite TARGET as fixed and the iterated route still does not. Reactions are removed in chunks
of decreasing size (ddmin); every accepted removal is logged and the current network saved, so a run can
be stopped at any time.

    python kinetic_variants/reduce_gap.py <model> <target metabolite> <out.json>
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

cache = {}
def gap(keep):
    key = frozenset(keep)
    if key in cache: return cache[key]
    ok = False
    try:
        kept = unblocked(sorted(keep))
        A, Y, N, ms = network(kept)
        if target in ms:
            ti = ms.index(target)
            it = kimonu.kinetic_modules(A, Y, N)
            acr_it = set(kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr).acr)
            if ti not in acr_it:
                _, acr_t, _, _ = one_shot_exact(A, Y, N, False)
                ok = ti in set(acr_t)
    except Exception as e:
        print("   predicate error:", type(e).__name__, str(e)[:80], flush=True)
    cache[key] = ok
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
t0 = time.time()
print(f"{M} {target}: {len(keep)} unblocked reactions; gap on full model: {gap(keep)} ({time.time() - t0:.0f}s)", flush=True)
n = 2
while len(keep) >= 2:
    chunks = [keep[i::n] for i in range(n)]            # interleaved, so each chunk samples the whole network
    removed = False
    for ch in chunks:
        trial = [r for r in keep if r not in set(ch)]
        if trial and gap(trial):
            keep = unblocked(sorted(trial)); n = max(n - 1, 2); removed = True
            r_, m_ = save(keep, "in progress")
            print(f"[{time.time() - t0:7.0f}s] removed a chunk of {len(ch)}: {r_} reactions, {m_} metabolites left", flush=True)
            break
    if not removed:
        if n >= len(keep): break
        n = min(len(keep), 2 * n)
        print(f"[{time.time() - t0:7.0f}s] no chunk removable, splitting finer: n = {n}", flush=True)
r_, m_ = save(keep, "1-minimal: removing any single reaction loses the gap")
print(f"done: {r_} reactions, {m_} metabolites ({time.time() - t0:.0f}s)", flush=True)
