"""Does the basis matter for the tolerance comparison? Tangent route, exact kernel, on small models.

Two exact comparisons of the same exact kernel T = ker(W Y^T), each for the half-decade tolerance grid:
    rref    Kimonu's RREF basis, absolute per element: zero row if max|r_i| <= tau, equal if max|r_i - r_j| <= tau
    ortho   an orthonormal basis, as the SVD arm uses, computed exactly through the projector P = T (T^T T)^-1 T^T:
            |r_i|^2 = P_ii and |r_i - r_j|^2 = P_ii + P_jj - 2 P_ij in every orthonormal basis of the same space;
            zero row if P_ii <= tau^2, equal if |r_i - r_j|^2 <= tau^2 max(P_ii, P_jj)
Both are exact (python-flint fmpq_mat). Reports, per tolerance, made-up and missed ACR/ACRR against the exact
tangent route for each comparison.

    python ortho_vs_rref.py <model> [...]
"""
import json, os, sys
from fractions import Fraction
import numpy as np, scipy.sparse as sp
import flint
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src"), os.path.join(B, "kinetic_variants")]
from kimonu.acr import readout
from kimonu.io import zero_complexes
from kimonu.kinetic import steady_state_map
from kimonu.linalg import exact, null_space_rows, positive_flux
from run_kimonu_crn import triplets_to_sparse, to_exact

TAUS = [float(10.0 ** e) for e in np.arange(-14, -0.5 + 1e-9, 0.5)]
TAU2 = [Fraction(10) ** int(round(2 * e)) for e in np.arange(-14, -0.5 + 1e-9, 0.5)]
fq = lambda x: flint.fmpq(x.numerator, x.denominator)
tofrac = lambda q: Fraction(int(q.p), int(q.q))

for M in sys.argv[1:]:
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    N = to_exact(triplets_to_sparse(crn["S"], nm, nr)); Ycsc = sp.csc_matrix(Y)
    flux = positive_flux(N, 0, scale=10**6)
    W, (_, ncol) = steady_state_map(A, Y, N, flux, True)
    J = {}
    for (i, c), w in W.items():
        if c >= nc: J[(i, nm + c - nc)] = w; continue
        for k, y in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
            J[(i, int(k))] = J.get((i, int(k)), Fraction(0)) + w * Fraction(str(y))
    rows = null_space_rows(exact(((i, k, v) for (i, k), v in J.items()), (nm, nm + ncol - nc)))
    ex = readout(Y, [tuple(sorted(r.items())) for r in rows[:nm]], zero_complexes(Y))
    acr_x = set(ex.acr); E = {(a, b) for g in ex.acrr for a in g for b in g if a < b}
    d = max((j for r in rows for j in r), default=-1) + 1
    T = flint.fmpq_mat(len(rows), d)
    for i, r in enumerate(rows):
        for j, v in r.items(): T[i, j] = fq(v)
    G = (T.transpose() * T).inv()
    Tm = flint.fmpq_mat(nm, d)
    for i in range(nm):
        for j in range(d): Tm[i, j] = T[i, j]
    P = Tm * G * Tm.transpose()                                  # metabolite block of the projector
    Pd = [tofrac(P[i, i]) for i in range(nm)]
    # candidate pairs for "equal": rows with equal exact RREF row, or close in float
    Rf = np.array([[float(rows[i].get(j, 0)) for j in range(d)] for i in range(nm)])
    res = {"rref": [], "ortho": []}
    pairs_rref, pairs_ortho = {}, {}
    for a in range(nm):
        for b in range(a + 1, nm):
            if np.abs(Rf[a] - Rf[b]).max() <= TAUS[-1] * 1.01:
                ra, rb = rows[a], rows[b]
                pairs_rref[(a, b)] = max(((ra.get(k, 0) - rb.get(k, 0)) ** 2 for k in set(ra) | set(rb)), default=Fraction(0))
            if Pd[a] or Pd[b]:
                d2 = Pd[a] + Pd[b] - 2 * tofrac(P[a, b])
                pairs_ortho[(a, b)] = (d2, max(Pd[a], Pd[b]))
    for k, (t, t2) in enumerate(zip(TAUS, TAU2)):
        for name in ("rref", "ortho"):
            if name == "rref":
                acr = {i for i in range(nm) if max((v * v for v in rows[i].values()), default=Fraction(0)) <= t2}
                found = {p for p, d2 in pairs_rref.items() if d2 <= t2}
            else:
                acr = {i for i in range(nm) if Pd[i] <= t2}
                found = {p for p, (d2, mx) in pairs_ortho.items() if d2 <= t2 * mx}
            found = {p for p in found if p[0] not in acr and p[1] not in acr}
            miss = len(acr_x - acr) + len({p for p in E - found if not (p[0] in acr and p[1] in acr)})
            fake = len(acr - acr_x) + len(found - E)
            res[name].append("right" if not (miss or fake) else "made up" if not miss else "missed" if not fake else "both")
    first_wrong = {n: next((f"{TAUS[k]:.0e}" for k, x in enumerate(v) if x != "right"), "none") for n, v in res.items()}
    print(M, "m", nm, "dim", d, "| first tolerance with an error:", first_wrong,
          "| missed ever:", {n: any(x in ("missed", "both") for x in v) for n, v in res.items()},
          "| same verdict at every tolerance:", res["rref"] == res["ortho"], flush=True)
