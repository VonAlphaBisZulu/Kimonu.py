"""Blocked reactions from one LP, against the FVA, on one network.

    max sum z   subject to   N r = 0,  r >= 0,  0 <= z <= 1,  z <= r

All rates are nonnegative and all lower bounds zero, so the steady states contain r = 0 and finite upper bounds
only scale: a reaction can carry flux under the bounds exactly when it can in the cone. Then z_j = 1 for every
reaction that can carry flux, and z_j = 0 for the blocked ones. Compared here with the plain capped version,
max sum r with r <= ub, which stops at a vertex and can leave a reaction that can carry flux at zero.

    python one_lp.py <sweep crn.json> <unblocked crn.json>
Prints one JSON line.
"""
import json, sys, time
t0 = time.perf_counter()
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
t_import = time.perf_counter() - t0

crn, ref = (json.load(open(f)) for f in sys.argv[1:3])
q, m = crn["n_reactions"], crn["n_metabolites"]
S = np.array(crn["S"], dtype=float)
N = sp.csr_matrix((S[:, 2], (S[:, 0].astype(int) - 1, S[:, 1].astype(int) - 1)), shape=(m, q))
lb, ub = np.array(crn["lb"], float), np.array(crn["ub"], float)
fva_kept = {tuple(r) if isinstance(r, list) else r for r in map(lambda x: x, ref["reactions"])}
ids = [tuple(r) if isinstance(r, list) else r for r in crn["reactions"]]
fva = np.array([i in fva_kept for i in ids])
assert fva.sum() == ref["n_reactions"], (fva.sum(), ref["n_reactions"])

t = time.perf_counter()
Aeq = sp.hstack([N, sp.csr_matrix((m, q))]).tocsr()
Aub = sp.hstack([-sp.identity(q), sp.identity(q)]).tocsr()
bounds = [(0, None if u > 0 else 0) for u in ub] + [(0, 1)] * q
c = np.r_[np.zeros(q), -np.ones(q)]
lp = linprog(c, A_ub=Aub, b_ub=np.zeros(q), A_eq=Aeq, b_eq=np.zeros(m), bounds=bounds, method="highs")
t_lp = time.perf_counter() - t
assert lp.status == 0, lp.message
z = lp.x[q:]; r = lp.x[:q]
one = z > 0.5

t = time.perf_counter()
cap = linprog(-np.ones(q), A_eq=N, b_eq=np.zeros(m), bounds=list(zip(np.zeros(q), ub)), method="highs")
t_cap = time.perf_counter() - t
assert cap.status == 0, cap.message
capped = cap.x > 1e-9

print(json.dumps({
    "model": sys.argv[1].split("/")[-1].split("_fixed")[0], "q": q, "min_lb": float(lb.min()),
    "fva_blocked": int((~fva).sum()),
    "lp_seconds": round(t_lp, 3), "import_seconds": round(t_import, 3),
    "lp_blocked": int((~one).sum()), "lp_wrong": int((one != fva).sum()),
    "z_between": int(((z > 1e-6) & (z < 1 - 1e-6)).sum()), "r_max": float(r.max()), "r_min_kept": float(r[one].min()),
    "capped_seconds": round(t_cap, 3), "capped_blocked": int((~capped).sum()), "capped_wrong": int((capped != fva).sum()),
}), flush=True)
