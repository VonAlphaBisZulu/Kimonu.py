"""Can two positive steady states with the same rate constants differ in a given log-linear quantity?

Let x, x' be positive steady states for rate constants k, z = ln x' - ln x and D = r' - r. Mass action
gives ln r' - ln r = G z (G_j = composition of reaction j's substrate complex), so sign(D) = sign(G z),
and both steady states give N D = 0. A quantity c . ln x (c = e_i for ACR of metabolite i, e_i - e_j
for a ratio, y_c - y_d for two complexes) differs between them iff c . z != 0. Hence

    exists z, D:  N D = 0,  sign(D_j) = sign((G z)_j) for every reaction j,  c . z >= 1

is necessary for two steady states that differ in c . ln x. If it is infeasible, c . ln x takes the
same value at every positive steady state, for every choice of rate constants: a global certificate.
If it is feasible the test is inconclusive (the sign pattern exists; the magnitudes may not).

Too weak as it stands: on EnvZ-OmpR, where Yp has ACR for all rate constants (Shinar-Feinberg), a
sign-compatible pair exists, so the test cannot certify it. Signs drop the magnitude coupling
D = r o (exp(G z) - 1) with r > 0 in ker N; certifying ACR needs that coupling (a nonconvex, exp-type
constraint), not just the sign pattern.

Signs are modelled with two binaries per reaction and indicator constraints (CPLEX).

    python branch_certificate.py <model|toy> <target> [...]      target: metabolite, or a/b for a ratio
"""
import json, os, sys, time
import numpy as np, scipy.sparse as sp
import cplex

EPS, BIG = 1e-3, 1e4

def certify(N, G, c, time_limit=3600, threads=2):
    """N: metabolites x reactions, G: reactions x metabolites (sparse), c: metabolite weights (dense)"""
    m, q = N.shape
    P = cplex.Cplex(); P.set_log_stream(None); P.set_results_stream(None); P.set_warning_stream(None)
    P.parameters.timelimit.set(time_limit); P.parameters.threads.set(threads)
    P.parameters.emphasis.mip.set(1)                                          # feasibility first
    z = [f"z{i}" for i in range(m)]; d = [f"d{j}" for j in range(q)]; v = [f"v{j}" for j in range(q)]
    pos = [f"p{j}" for j in range(q)]; neg = [f"n{j}" for j in range(q)]
    P.variables.add(names=z, lb=[-BIG] * m, ub=[BIG] * m)
    P.variables.add(names=d, lb=[-cplex.infinity] * q, ub=[cplex.infinity] * q)
    P.variables.add(names=v, lb=[-cplex.infinity] * q, ub=[cplex.infinity] * q)
    P.variables.add(names=pos + neg, types="B" * (2 * q))
    Nc, Gr = N.tocsr(), G.tocsr()
    rows = [cplex.SparsePair([d[j] for j in Nc.indices[Nc.indptr[i]:Nc.indptr[i + 1]]],
                             list(Nc.data[Nc.indptr[i]:Nc.indptr[i + 1]])) for i in range(m)]
    P.linear_constraints.add(lin_expr=rows, senses="E" * m, rhs=[0] * m)                     # N D = 0
    rows = [cplex.SparsePair([v[j]] + [z[i] for i in Gr.indices[Gr.indptr[j]:Gr.indptr[j + 1]]],
                             [1.0] + list(-Gr.data[Gr.indptr[j]:Gr.indptr[j + 1]])) for j in range(q)]
    P.linear_constraints.add(lin_expr=rows, senses="E" * q, rhs=[0] * q)                     # v = G z
    P.linear_constraints.add(lin_expr=[cplex.SparsePair([pos[j], neg[j]], [1, 1]) for j in range(q)],
                             senses="L" * q, rhs=[1] * q)
    nz = np.flatnonzero(c)
    P.linear_constraints.add(lin_expr=[cplex.SparsePair([z[i] for i in nz], list(c[nz]))], senses="G", rhs=[1])
    ind = P.indicator_constraints
    for j in range(q):
        ind.add(indvar=pos[j], complemented=0, lin_expr=cplex.SparsePair([d[j]], [1]), sense="G", rhs=1)
        ind.add(indvar=pos[j], complemented=0, lin_expr=cplex.SparsePair([v[j]], [1]), sense="G", rhs=EPS)
        ind.add(indvar=neg[j], complemented=0, lin_expr=cplex.SparsePair([d[j]], [1]), sense="L", rhs=-1)
        ind.add(indvar=neg[j], complemented=0, lin_expr=cplex.SparsePair([v[j]], [1]), sense="L", rhs=-EPS)
    # neither: D_j = 0 and (G z)_j = 0, via a third binary that is on exactly when both p and n are off
    zer = [f"o{j}" for j in range(q)]
    P.variables.add(names=zer, types="B" * q)
    P.linear_constraints.add(lin_expr=[cplex.SparsePair([pos[j], neg[j], zer[j]], [1, 1, 1]) for j in range(q)],
                             senses="E" * q, rhs=[1] * q)
    for j in range(q):
        for var in (d[j], v[j]):
            ind.add(indvar=zer[j], complemented=0, lin_expr=cplex.SparsePair([var], [1]), sense="E", rhs=0)
    t0 = time.time(); P.solve(); st = P.solution.get_status_string()
    return st, time.time() - t0

if __name__ == "__main__":
    B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src")); sys.path.insert(0, os.path.join(os.path.dirname(B), "tests"))
    from run_kimonu_crn import triplets_to_sparse
    name, targets = sys.argv[1], sys.argv[2:]
    if name in ("cubic", "envz"):
        from test_kinetic import _crn, ENVZ, ENVZ_SPECIES
        R, mets = ([("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")], ["A"]) if name == "cubic" else (ENVZ, ENVZ_SPECIES)
        A, Y, N, cx = _crn(R, mets); A, Y, N = sp.csr_matrix(A), sp.csr_matrix(Y, dtype=float), sp.csr_matrix(N, dtype=float)
    else:
        crn = json.load(open(f"{B}/unblocked/{name}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
        A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc); N = triplets_to_sparse(crn["S"], nm, nr)
        mets = crn["metabolites"]
    G = ((A < 0).astype(float).T @ Y.T).tocsr()
    for t in targets:
        c = np.zeros(len(mets))
        a, _, b = t.partition("/"); c[mets.index(a)] = 1
        if b: c[mets.index(b)] = -1
        st, sec = certify(N, G, c)
        verdict = "CERTIFIED: same value at every positive steady state" if "infeasible" in st else "inconclusive: sign pattern for two different values exists" if "feasible" in st or "optimal" in st else "undecided"
        print(f"{name:14s} {t:28s} {st:28s} {sec:7.1f}s  -> {verdict}", flush=True)
