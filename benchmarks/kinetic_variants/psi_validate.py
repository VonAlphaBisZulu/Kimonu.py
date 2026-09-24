"""Numerical check of ψ-nullspace couplings: two positive mass-action steady states, same rate constants.

Rate constants are built from a positive flux v* and random concentrations x*, so x* is a steady
state by construction. A second steady state is found by Newton's method in log-concentrations after
shifting the conserved totals. Every pair the ψ-nullspace calls coupled must keep its ψ ratio across
the two states; pairs it does not couple should in general change.
"""
import json, sys, os
import numpy as np, scipy.sparse as sp, scipy.linalg as la
from fractions import Fraction
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
exec(open(os.path.join(B, "kinetic_variants/psi_nullspace.py")).read().split("DETACH = ")[0])
M = sys.argv[1]; rng = np.random.default_rng(7)
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc).tocsc(); Sx = to_exact(triplets_to_sparse(crn["S"], nm, nr))
Nf = triplets_to_sparse(crn["S"], nm, nr).toarray()
Kx = sd.sparse_nullspace(Sx)
if not isinstance(Kx, sd.ExactCOO):
    Kc = Kx.tocoo(); Kx = type("K", (), dict(rows=Kc.row, cols=Kc.col, data=Kc.data, denom=1))()
sub = np.zeros(nr, int)
for i, j, v in zip(*sp.find(A)):
    if v < 0: sub[j] = i
vstar = np.array([float(x) for x in exact_positive_flux(crn, Sx, Kx, 1)])
Yd = Y.toarray(); Ys = Yd[:, sub]                     # metabolites x reactions: substrate compositions
zstar = rng.normal(0, 0.5, nm)                          # log x*
k = vstar / np.exp(Ys.T @ zstar)
flux = lambda z: k * np.exp(Ys.T @ z)
print("residual at x*", np.abs(Nf @ flux(zstar)).max())
Lc = la.null_space(Nf.T).T                              # conservation laws  Lc x = const
U, s, Vt = la.svd(Nf, full_matrices=False); rk = int((s > 1e-9 * s[0]).sum()); P = U[:, :rk].T   # independent steady-state equations
ctot = Lc @ np.exp(zstar); target = ctot * (1 + 0.2 * rng.random(len(ctot)))
vs = np.maximum(vstar, 1e-12)
def resid(z, c):
    return np.concatenate([P @ (Nf @ flux(z)) / np.abs(P @ Nf @ vs + 1).max(), (Lc @ np.exp(z) - c) / (np.abs(c) + 1)])
z = zstar.copy(); ok = True
for lam in np.linspace(0.02, 1, 50):
    c = ctot + lam * (target - ctot)
    for it in range(40):
        v = flux(z); F = resid(z, c); nF = np.linalg.norm(F)
        if nF < 1e-11: break
        J = np.vstack([P @ Nf @ (v[:, None] * Ys.T) / np.abs(P @ Nf @ vs + 1).max(), (Lc * np.exp(z)[None, :]) / (np.abs(c) + 1)[:, None]])
        dz = la.lstsq(J, -F, lapack_driver="gelsy")[0]
        t = 1.0
        while t > 1e-4 and not np.linalg.norm(resid(z + t * dz, c)) < nF: t /= 2
        z = z + t * dz
    if not np.isfinite(z).all() or np.linalg.norm(resid(z, c)) > 1e-8:
        ok = False; print("continuation failed at", lam, np.linalg.norm(resid(z, c))); break
if not ok: sys.exit(1)
print("second state: residual", np.abs(Nf @ flux(z)).max(), "conservation error", np.abs(Lc @ np.exp(z) - target).max(), "moved", np.abs(z - zstar).max())
lpsi1 = Yd.T @ zstar; lpsi2 = Yd.T @ z                   # log ψ at both states
# ψ-nullspace blocks (detached), same construction as psi_nullspace.py
Nent = {}
for i, j, v in crn["S"]: Nent.setdefault(int(j) - 1, []).append((int(i) - 1, Fraction(str(v))))
zero = [c for c in range(nc) if Y[:, c].nnz == 0]
vx = exact_positive_flux(crn, Sx, Kx, 2); W = {}; extra = {}
for r_ in range(nr):
    s = int(sub[r_]); col = extra.setdefault(r_, nc + len(extra)) if s in zero else s
    for i, a in Nent.get(r_, []): W[(i, col)] = W.get((i, col), Fraction(0)) + vx[r_] * a
ncol = nc + len(extra)
psi = kimonu.concordance_from_matrices(sp.identity(ncol, format="csr", dtype=np.int64), sd.RationalMatrix.from_fractions(((i, j, x) for (i, j), x in W.items() if x != 0), (nm, ncol)))
res = kimonu.concordance_from_matrices(A, Sx); ko = kimonu.kinetic_modules(A, res.modules, res.balanced, Y=None)
of = {c: frozenset(m) for m in ko.modules for c in m}
shift = lpsi2 - lpsi1                                   # a coupled block must shift by one common amount
worst_coupled = 0; worst_extra = 0; n_extra = 0
for m in psi.modules:
    m = [c for c in m if c < nc and c not in zero]
    if len(m) < 2: continue
    spread = shift[m].max() - shift[m].min(); worst_coupled = max(worst_coupled, spread)
    core = max((of[c] for c in m), key=len)
    ex = [c for c in m if c not in core]
    if ex:
        n_extra += len(ex); worst_extra = max(worst_extra, np.abs(shift[ex] - shift[list(core)[0]]).max())
ns = [c for c in range(nc) if c not in zero]
rnd = [abs(shift[a] - shift[b]) for a, b in rng.integers(0, len(ns), (2000, 2))]
print(f"{M}: coupled blocks worst log-ratio change {worst_coupled:.2e}; extra couplings ({n_extra} complexes) worst {worst_extra:.2e}; random pairs median change {np.median(rnd):.2e}")
