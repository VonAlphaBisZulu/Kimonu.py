"""Search for two positive steady states (same rate constants) that differ in one metabolite.

Unknowns a = ln r (rates at x) and z = ln x' - ln x with z_target fixed to +-delta. Residuals
N exp(a) = 0 and N exp(a + G z) = 0, scaled per metabolite, minimised by sparse least squares from
the interior flux r* (a = ln r*, z = delta e_target). A zero residual is a candidate second steady
state; it is then checked by the exact criterion: with z fixed, an LP must find r >= 1 in ker N with
exp(G z) o r in ker N.

    python branch_search.py <model> <metabolite> [delta ...]
"""
import json, os, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import least_squares, linprog
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
from run_kimonu_crn import triplets_to_sparse, to_exact
from kimonu.kinetic import positive_flux

M, target = sys.argv[1], sys.argv[2]
deltas = [float(d) for d in sys.argv[3:]] or [0.5, -0.5, 2.0, -2.0]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = triplets_to_sparse(crn["S"], nm, nr).tocsr()
G = ((A < 0).astype(float).T @ Y.T).tocsr()
mets = crn["metabolites"]; t = mets.index(target)
r0 = np.array([float(v) for v in positive_flux(to_exact(N), 0)])
scale = 1 / np.maximum(np.abs(N) @ r0, 1e-12)                       # per-metabolite turnover
free = np.array([i for i in range(nm) if i != t])

def unpack(p, delta):
    a = p[:nr]; z = np.zeros(nm); z[free] = p[nr:]; z[t] = delta
    return a, z
def res(p, delta):
    a, z = unpack(p, delta)
    return np.concatenate([scale * (N @ np.exp(a)), scale * (N @ np.exp(a + G @ z)), [a[0] - np.log(r0[0])]])
def jac(p, delta):
    a, z = unpack(p, delta); e1, e2 = np.exp(a), np.exp(a + G @ z)
    S = sp.diags(scale)
    top = sp.hstack([S @ N @ sp.diags(e1), sp.csr_matrix((nm, len(free)))])
    bot = sp.hstack([S @ N @ sp.diags(e2), S @ N @ sp.diags(e2) @ G[:, free]])
    last = sp.csr_matrix(([1.0], ([0], [0])), shape=(1, nr + len(free)))
    return sp.vstack([top, bot, last]).tocsr()

def exact_check(z):
    """LP: r >= 1, N r = 0, N (exp(G z) o r) = 0"""
    e = np.exp(G @ z)
    lp = linprog(np.zeros(nr), A_eq=sp.vstack([N, N @ sp.diags(e)]).tocsr(), b_eq=np.zeros(2 * nm),
                 bounds=[(1, None)] * nr, method="highs")
    return lp.status == 0

print(f"{M}: {nm} metabolites, {nr} reactions; target {target}", flush=True)
for delta in deltas:
    p0 = np.concatenate([np.log(r0), np.zeros(len(free))])
    t0 = time.time()
    sol = least_squares(res, p0, jac=jac, args=(delta,), method="trf", tr_solver="lsmr", x_scale="jac",
                        max_nfev=300, ftol=1e-14, xtol=1e-14, gtol=1e-14)
    a, z = unpack(sol.x, delta)
    rr = np.abs(res(sol.x, delta)).max()
    ok = exact_check(z) if rr < 1e-8 else False
    print(f"delta {delta:+.2f}: max scaled residual {rr:.2e}, |z| max {np.abs(z).max():.2f}, "
          f"LP check {'FEASIBLE: second steady state found' if ok else 'not feasible'} ({time.time() - t0:.0f}s)", flush=True)
