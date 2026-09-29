"""Two positive steady states that differ in c . ln x, with the magnitudes kept (global MINLP, SCIP).

Variables a = ln r (rates at x) and z = ln x' - ln x; the rates at x' are exp(a + G z). Both are steady
states for the same rate constants iff N exp(a) = 0 and N exp(a + G z) = 0. Rate constants are free,
so a ranges over all positive steady-state fluxes; a_0 = 0 fixes the scale. Asking for c . z >= delta
and proving infeasibility within the box |a|, |z| <= L certifies that c . ln x cannot change by delta
or more between two steady states whose concentrations differ by at most e^L per metabolite.

    python branch_global.py <cubic|envz|model> <target> [delta] [L]
"""
import json, os, sys, time
import numpy as np, scipy.sparse as sp
from pyscipopt import Model, exp, quicksum

def build(N, G, c, delta, L, time_limit=3600):
    N, G = sp.csr_matrix(N), sp.csr_matrix(G); m, q = N.shape
    M = Model(); M.hideOutput(); M.setParam("limits/time", time_limit)
    a = [M.addVar(lb=-L, ub=L, name=f"a{j}") for j in range(q)]
    z = [M.addVar(lb=-L, ub=L, name=f"z{i}") for i in range(m)]
    r = [M.addVar(lb=np.exp(-L), ub=np.exp(L), name=f"r{j}") for j in range(q)]
    rp = [M.addVar(lb=np.exp(-2 * L * 10), ub=np.exp(2 * L * 10), name=f"rp{j}") for j in range(q)]
    for j in range(q):
        M.addCons(r[j] == exp(a[j]))
        gz = quicksum(G.data[k] * z[G.indices[k]] for k in range(G.indptr[j], G.indptr[j + 1]))
        M.addCons(rp[j] == exp(a[j] + gz))
    for i in range(m):
        idx = range(N.indptr[i], N.indptr[i + 1])
        M.addCons(quicksum(N.data[k] * r[N.indices[k]] for k in idx) == 0)
        M.addCons(quicksum(N.data[k] * rp[N.indices[k]] for k in idx) == 0)
    M.addCons(a[0] == 0)
    M.addCons(quicksum(float(c[i]) * z[i] for i in np.flatnonzero(c)) >= delta)
    M.setObjective(0)
    return M, z

def build_poly(N, G, c_index, ratio, L, time_limit=3600):
    """Polynomial form: w = x'/x (bounded in [e^-L, e^L]), r > 0 with N r = 0, N (r o w^G) = 0, and
    w_target >= ratio. G is integer here (molecularities), so w^G is a monomial."""
    N, G = sp.csr_matrix(N), sp.csr_matrix(G); m, q = N.shape
    M = Model(); M.hideOutput(); M.setParam("limits/time", time_limit)
    w = [M.addVar(lb=np.exp(-L), ub=np.exp(L), name=f"w{i}") for i in range(m)]
    r = [M.addVar(lb=1e-3, ub=1e3, name=f"r{j}") for j in range(q)]
    rp = []
    for j in range(q):
        mono = r[j]
        for k in range(G.indptr[j], G.indptr[j + 1]):
            mono = mono * w[G.indices[k]] ** int(round(G.data[k]))
        v = M.addVar(lb=0, ub=None, name=f"rp{j}"); M.addCons(v == mono); rp.append(v)
    for i in range(m):
        idx = range(N.indptr[i], N.indptr[i + 1])
        M.addCons(quicksum(N.data[k] * r[N.indices[k]] for k in idx) == 0)
        M.addCons(quicksum(N.data[k] * rp[N.indices[k]] for k in idx) == 0)
    M.addCons(r[0] == 1)
    M.addCons(w[c_index] >= ratio)
    M.setObjective(0)
    return M, w

if __name__ == "__main__" and len(sys.argv) > 5 and sys.argv[5] == "poly":
    B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "tests"))
    from test_kinetic import _crn, ENVZ, ENVZ_SPECIES
    name, target, ratio, L = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    R, mets = ([("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")], ["A"]) if name == "cubic" else (ENVZ, ENVZ_SPECIES)
    A, Y, N, cx = _crn(R, mets)
    G = (A.toarray() < 0).astype(float).T @ Y.toarray().T
    M, w = build_poly(N.toarray().astype(float), G, mets.index(target), ratio, L, time_limit=900)
    t0 = time.time(); M.optimize(); st = M.getStatus()
    ws = {mets[i]: round(M.getVal(w[i]), 3) for i in range(len(mets))} if M.getNSols() else None
    print(f"[poly] {name} {target} x'/x >= {ratio} (box e^+-{L}): {st} in {time.time() - t0:.1f}s" + (f"  w={ws}" if ws else ""))
    sys.exit()

if __name__ == "__main__":
    B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "tests"))
    name, target = sys.argv[1], sys.argv[2]
    delta = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1; L = float(sys.argv[4]) if len(sys.argv) > 4 else 5
    from test_kinetic import _crn, ENVZ, ENVZ_SPECIES
    R, mets = ([("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")], ["A"]) if name == "cubic" else (ENVZ, ENVZ_SPECIES)
    A, Y, N, cx = _crn(R, mets)
    G = (A.toarray() < 0).astype(float).T @ Y.toarray().T
    c = np.zeros(len(mets)); a_, _, b_ = target.partition("/"); c[mets.index(a_)] = 1
    if b_: c[mets.index(b_)] = -1
    M, z = build(N.toarray().astype(float), G, c, delta, L)
    t0 = time.time(); M.optimize(); st = M.getStatus()
    zs = {mets[i]: round(M.getVal(z[i]), 3) for i in range(len(mets))} if M.getNSols() else None
    print(f"{name} {target} delta={delta} L={L}: {st} in {time.time() - t0:.1f}s" + (f"  z={zs}" if zs else ""))
