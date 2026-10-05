"""The tangent route in floating point against the exact one, on one genome-scale model.

Same interior flux r* (integer, seed 0) and the same J = N diag(r*) S^T Y^T; the exact kernel comes from
straindesign's sparse rational nullspace, the float kernel from a dense SVD of J in float64. Rank is decided
at three relative tolerances (singular values below tol * sigma_max count as zero; the first is numpy's
default, max(m, n) * eps). ACR is a zero row of the float kernel basis (norm below 1e-8; the basis is
orthonormal, so rows have norm <= 1), ACRR two rows that agree to 1e-8 relative to their norm. A first
version compared rows to an absolute 1e-8 and called rows of norm ~1e-4 that differ in the fifth digit
equal; the relative test removes those, so what remains are the errors of the float rank decision.
Float rank and float row equality each need a tolerance; to be fair to floating point, three rank tolerances
are crossed with three row tolerances (1e-8, 1e-6, 1e-4), and the best of the nine per model is reported as
well ("oracle": a tolerance one could only pick knowing the exact answer). Attached uptakes only.

    python tangent_float.py <model>   -> json line, also written to $TF_OUT/<model>.json
"""
import json, os, sys, time
from fractions import Fraction
from itertools import combinations
import numpy as np, scipy.sparse as sp
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import straindesign as sd
from run_kimonu_crn import triplets_to_sparse, to_exact
from kimonu.linalg import entries
from kinetic_variants.tangent_fast import integer_positive_flux

OUT = os.environ.get("TF_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_tangent_float")
ROW_TOL = 1e-8


def jacobian(A, Y, N, flux):
    """J = N diag(r*) G with G row j = composition of reaction j's substrate complex (attached uptakes: zero row)"""
    Ac = sp.csc_matrix(A); Yc = sp.csc_matrix(Y)
    sub = {}
    for j in range(Ac.shape[1]):
        for i, a in zip(Ac.indices[Ac.indptr[j]:Ac.indptr[j + 1]], Ac.data[Ac.indptr[j]:Ac.indptr[j + 1]]):
            if a < 0: sub[j] = int(i)
    triples, (m, q) = entries(N)
    J = {}
    for i, j, v in triples:
        if j not in sub: continue
        c = sub[j]
        for k, w in zip(Yc.indices[Yc.indptr[c]:Yc.indptr[c + 1]], Yc.data[Yc.indptr[c]:Yc.indptr[c + 1]]):
            J[(i, int(k))] = J.get((i, int(k)), Fraction(0)) + flux[j] * v * Fraction(str(w))
    return {key: x for key, x in J.items() if x != 0}, m


def readout_exact(J, m):
    K = sd.sparse_nullspace(sd.RationalMatrix.from_fractions(((i, k, x) for (i, k), x in J.items()), (m, m)))
    rows = {}
    for i, j, v in entries(K)[0]:
        rows.setdefault(i, []).append((j, v))
    acr = {i for i in range(m) if i not in rows}
    groups = {}
    for i, r in rows.items():
        groups.setdefault(tuple(sorted(r)), []).append(i)
    pairs = {(a, b) for g in groups.values() for a, b in combinations(sorted(g), 2)}
    return acr, pairs, K.shape[1]


def readout_float(s, Vt, tol_rel, row_tol=ROW_TOL):
    rank = int((s > tol_rel * s[0]).sum()) if s.size and s[0] > 0 else 0
    T = Vt[rank:].T                                      # orthonormal basis of the float kernel, m x (m - rank)
    m = T.shape[0]
    norms = np.linalg.norm(T, axis=1)
    acr = {i for i in range(m) if norms[i] < row_tol}
    rest = [i for i in range(m) if i not in acr]
    # equal rows: sort by a random projection, then compare neighbours within the tolerance window
    proj = T[rest] @ np.random.default_rng(0).standard_normal(T.shape[1]) if T.shape[1] else np.zeros(len(rest))
    order = np.argsort(proj)
    pairs = set()
    for a in range(len(order)):
        b = a + 1
        while b < len(order) and proj[order[b]] - proj[order[a]] < row_tol * 10 * norms[rest[order[a]]] * np.sqrt(T.shape[1]):
            i, j = rest[order[a]], rest[order[b]]
            if np.linalg.norm(T[i] - T[j]) < row_tol * max(norms[i], norms[j]):
                pairs.add((min(i, j), max(i, j)))
            b += 1
    return acr, pairs, m - rank, s


if __name__ == "__main__":
    M = sys.argv[1]
    os.makedirs(OUT, exist_ok=True)
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    flux = integer_positive_flux(N, seed=0)
    J, m = jacobian(A, Y, N, flux)
    t0 = time.time(); acr_x, pairs_x, dim_x = readout_exact(J, m); t_exact = time.time() - t0
    Jd = np.zeros((m, m))
    for (i, k), x in J.items(): Jd[i, k] = float(x)
    out = dict(model=M, m=m, exact=dict(seconds=t_exact, dim_T=dim_x, acr=len(acr_x), acrr_pairs=len(pairs_x)), float={})
    t0 = time.time(); U, s, Vt = np.linalg.svd(Jd); t_svd = time.time() - t0
    out["float_svd_seconds"] = t_svd
    for name, tol in (("numpy_default", m * np.finfo(float).eps), ("1e-10", 1e-10), ("1e-6", 1e-6)):
        for row_tol in (1e-8, 1e-6, 1e-4):
            acr_f, pairs_f, dim_f, _ = readout_float(s, Vt, tol, row_tol)
            out["float"][f"{name}|{row_tol:g}"] = dict(tol=tol, row_tol=row_tol, dim_T=dim_f,
                                      acr_fp=len(acr_f - acr_x), acr_fn=len(acr_x - acr_f),
                                      acrr_fp=len(pairs_f - pairs_x), acrr_fn=len(pairs_x - pairs_f))
    wrong = lambda v: v["acr_fp"] + v["acr_fn"] + v["acrr_fp"] + v["acrr_fn"]
    best = min(out["float"], key=lambda k: wrong(out["float"][k]))
    out["oracle"] = dict(setting=best, **out["float"][best])
    out["sigma_min_nonzero_over_max"] = float(s[s > 0][-1] / s[0]) if (s > 0).any() else None
    out["rate_span_log10"] = float(np.log10(float(max(flux) / min(flux))))
    json.dump(out, open(f"{OUT}/{M}.json", "w"), indent=1)
    print(json.dumps(out), flush=True)
