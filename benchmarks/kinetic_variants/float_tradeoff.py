"""Float64 against exact, as a trade-off: relations missed and relations made up, over the tolerance.

The tangent route (detached uptakes) is computed exactly (kimonu.tangent_route) and in float64 at the same
interior flux r* (integer, scale 10^6, seed 0): J = W Y^T assembled in floats, its kernel from a dense SVD.
Float needs two decisions the exact route does not: which singular values count as zero (rank tolerance,
relative to sigma_max) and when two rows of the kernel basis count as equal (row tolerance). A metabolite is
float-ACR when its row has norm below the row tolerance (the basis is orthonormal, so rows have norm <= 1);
two metabolites are float-ACRR when their nonzero rows differ by less than the row tolerance times the larger
norm. Every float ACR metabolite and ACRR pair is compared with the exact ones: a missed exact relation is a
false negative, a reported relation the exact route does not have is a false positive.

    python float_tradeoff.py <model>   -> $FT_OUT/<model>.json
"""
import json, os, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.kinetic import steady_state_map
from kimonu.linalg import positive_flux
from run_kimonu_crn import triplets_to_sparse, to_exact

OUT = os.environ.get("FT_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_float_tradeoff")
ROW_TOLS = [10.0 ** e for e in np.arange(-14, -0.5 + 1e-9, 0.5)]
RANK_TOLS = ["numpy", 1e-12, 1e-10, 1e-8, 1e-6]

M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
free = np.zeros(nm, bool); free[:mets.index("E_1") if "E_1" in mets else nm] = True

t0 = time.time()
ex = kimonu.tangent_route(A, Y, N, detach_uptakes=True)
t_exact = time.time() - t0
acr_x = set(int(i) for i in ex.acr)
E = np.zeros((nm, nm), bool)
for a, b in ex.acrr:
    E[min(a, b), max(a, b)] = True

# the same J in float64: column of complex c spread over its metabolites, detached sources after them
flux = positive_flux(N, 0, scale=10**6)
W, (_, ncol) = steady_state_map(A, Y, N, flux, True)
Yc = sp.csc_matrix(Y)
J = np.zeros((nm, nm + ncol - nc))
for (i, c), w in W.items():
    if c >= nc:
        J[i, nm + c - nc] += float(w)
        continue
    for k, y in zip(Yc.indices[Yc.indptr[c]:Yc.indptr[c + 1]], Yc.data[Yc.indptr[c]:Yc.indptr[c + 1]]):
        J[i, k] += float(w) * float(y)
t0 = time.time()
U, s, Vt = np.linalg.svd(J, full_matrices=True)
t_svd = time.time() - t0
eps_tol = max(J.shape) * np.finfo(float).eps

out = dict(model=M, m=nm, exact=dict(seconds=t_exact, dim_T=ex.dimension, acr=len(acr_x), acrr_pairs=int(E.sum()),
                                     acr_free=int(sum(free[i] for i in acr_x)),
                                     acrr_pairs_free=int(E[np.ix_(free, free)].sum())),
           svd_seconds=t_svd, sigma_ratio=float(s[s > eps_tol * s[0]][-1] / s[0]), row_tols=ROW_TOLS, rank={})
iu = np.triu_indices(nm, 1)
for rt in RANK_TOLS:
    tol = eps_tol if rt == "numpy" else rt
    rank = int((s > tol * s[0]).sum())
    K = Vt[rank:].T[:nm]                               # kernel basis, metabolite rows
    norm = np.linalg.norm(K, axis=1)
    G = K @ K.T                                        # pairwise distances from the Gram matrix
    D = np.sqrt(np.maximum(norm[:, None] ** 2 + norm[None, :] ** 2 - 2 * G, 0))
    R = (D / np.maximum(np.maximum(norm[:, None], norm[None, :]), 1e-300))[iu]
    Ex = E[iu]; nz = norm
    res = []
    for tau in ROW_TOLS:
        acr_f = set(np.flatnonzero(nz < tau).tolist())
        keep = (nz[iu[0]] >= tau) & (nz[iu[1]] >= tau)
        found = keep & (R < tau)
        res.append(dict(acr_fp=len(acr_f - acr_x), acr_fn=len(acr_x - acr_f),
                        acrr_fp=int((found & ~Ex).sum()), acrr_fn=int((Ex & ~found).sum()),
                        acr_fp_free=int(sum(free[i] for i in acr_f - acr_x)), acr_fn_free=int(sum(free[i] for i in acr_x - acr_f)),
                        acrr_fp_free=int((found & ~Ex & free[iu[0]] & free[iu[1]]).sum()),
                        acrr_fn_free=int((Ex & ~found & free[iu[0]] & free[iu[1]]).sum())))
    out["rank"][str(rt)] = dict(tol=tol, dim_T=nm + ncol - nc - rank, sweep=res)
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"))
print(json.dumps({k: v for k, v in out.items() if k != "rank"}), flush=True)
