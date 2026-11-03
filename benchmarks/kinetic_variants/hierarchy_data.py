"""Per model: the iterated route (proven) against the tangent route (local ACR), both attached uptakes.

The iterated route is kimonu.kinetic_modules + acr_acrr from src (terminal complexes join where coupled).
The tangent route reads T = ker(N diag(r*) S^T Y^T) at an integer interior flux (tangent_fast, scale 10^6);
modules are equal rows of Y^T T over all complexes, as in src, ACR the zero rows of T, ACRR groups equal
nonzero rows. Every coupled pair, ACR metabolite and ACRR pair of the iterated route must also be in T
(proven <= true <= local); the script records any violation instead of asserting, so a bug shows up in
the table.

    python kinetic_variants/hierarchy_data.py <model>   -> kinetic_variants/hierarchy/<model>.json
"""
import json, os, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import straindesign as sd
import kimonu
from kimonu.linalg import entries
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.tangent_fast import integer_positive_flux


def tangent(A, Y, N, seed=0):
    """T at an integer interior flux; module groups over all complexes, ACR, ACRR groups (indices)."""
    n = A.shape[0]; Ac = sp.csc_matrix(A); Yc = sp.csc_matrix(Y)
    flux = integer_positive_flux(N, seed)
    sub = {}
    for j in range(Ac.shape[1]):
        for i, a in zip(Ac.indices[Ac.indptr[j]:Ac.indptr[j + 1]], Ac.data[Ac.indptr[j]:Ac.indptr[j + 1]]):
            if a < 0: sub[j] = int(i)
    ycol = {c: [(int(i), Fraction(str(v))) for i, v in zip(Yc.indices[Yc.indptr[c]:Yc.indptr[c + 1]],
                                                          Yc.data[Yc.indptr[c]:Yc.indptr[c + 1]])] for c in range(n)}
    triples, (m, q) = entries(N)
    byr = {}
    for i, j, v in triples: byr.setdefault(j, []).append((i, v))
    J = {}
    for j, s in sub.items():                      # J = N diag(r*) S^T Y^T, column k = metabolite
        for i, v in byr.get(j, ()):
            for k, w in ycol[s]:
                J[(i, k)] = J.get((i, k), Fraction(0)) + flux[j] * v * w
    Jx = sd.RationalMatrix.from_fractions(((i, k, x) for (i, k), x in J.items() if x != 0), (m, m))
    K = sd.sparse_nullspace(Jx)
    rows = {}
    for i, j, v in entries(K)[0]: rows.setdefault(i, {})[j] = v
    key = lambda d: tuple(sorted((c, v) for c, v in d.items() if v != 0))
    groups = {}
    for c in range(n):                            # every complex, terminal ones included
        acc = {}
        for i, w in ycol[c]:
            for j, v in rows.get(i, {}).items(): acc[j] = acc.get(j, 0) + w * v
        groups.setdefault(key(acc), set()).add(c)
    acr = [i for i in range(m) if not rows.get(i)]
    mets = {}
    for i in range(m):
        if rows.get(i): mets.setdefault(key(rows[i]), []).append(i)
    return list(groups.values()), acr, [g for g in mets.values() if len(g) > 1], K.shape[1]


def pairs(groups):
    return {(a, b) for g in groups for a in g for b in g if a < b}


M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
mets = crn["metabolites"]; mets = [m if isinstance(m, str) else m.get("id") for m in mets]
free = set(range(mets.index("E_1") if "E_1" in mets else nm))

t0 = time.time()
it = kimonu.kinetic_modules(A, Y, S, detach_uptakes=False)
ar = kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr)
t_it = time.time() - t0
t0 = time.time(); mods, acr_t, acrr_g, dim = tangent(A, Y, S); t_tan = time.time() - t0

cp_it, cp_t = pairs([sorted(m) for m in it.modules]), pairs([sorted(m) for m in mods])
acr_it = set(ar.acr); acr_t = set(acr_t)
acrr_it = {tuple(sorted(p)) for p in ar.acrr}; acrr_t = pairs([sorted(g) for g in acrr_g])
nontriv = lambda gs: [len(g) for g in gs if len(g) > 1]
out = dict(
    model=M, n_complexes=nc, n_metabolites=nm, tangent_dim=dim, seconds_iterated=t_it, seconds_tangent=t_tan,
    it=dict(modules=len(nontriv(it.modules)), giant=max(nontriv(it.modules), default=1), pairs=len(cp_it),
            acr=len(acr_it), acrr_pairs=len(acrr_it),
            acr_free=len(acr_it & free), acrr_pairs_free=sum(a in free and b in free for a, b in acrr_it)),
    tan=dict(modules=len(nontriv(mods)), giant=max(nontriv(mods), default=1), pairs=len(cp_t),
             acr=len(acr_t), acrr_pairs=len(acrr_t),
             acr_free=len(acr_t & free), acrr_pairs_free=sum(a in free and b in free for a, b in acrr_t)),
    violations=dict(pairs=len(cp_it - cp_t), acr=len(acr_it - acr_t), acrr_pairs=len(acrr_it - acrr_t)))
out["gap"] = any(out["it"][k] != out["tan"][k] for k in ("pairs", "acr", "acrr_pairs"))
os.makedirs(f"{B}/kinetic_variants/hierarchy", exist_ok=True)
json.dump(out, open(f"{B}/kinetic_variants/hierarchy/{M}.json", "w"), indent=1)
print(json.dumps(out), flush=True)
