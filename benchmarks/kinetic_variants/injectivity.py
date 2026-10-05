"""Injectivity: at most one positive steady state per compatibility class, for every choice of rate constants.

The mass-action map f_k(x) = N (k o psi(x)) is injective on every compatibility class (x + S) and for every
k > 0 iff there are no

    sigma in S = im N,  D in ker N,  dxi in R^m, dxi != 0   with   sign(dxi) = sign(sigma),  sign(D) = sign(G dxi)

(G = S^T Y^T, row j = composition of reaction j's substrate complex; zero for uptakes, whose rate is k_j).
Necessity: x != x' in one class with f_k(x) = f_k(x') give sigma = x - x', dxi = ln x - ln x' (same signs,
ln is monotone) and D = k o (psi(x) - psi(x')) in ker N with sign(D_j) = sign((G dxi)_j). Sufficiency: from
such a triple, x'_i = sigma_i / (exp(dxi_i) - 1), x_i = x'_i exp(dxi_i) where sigma_i != 0 (positive, since
the signs agree), x_i = x'_i = 1 elsewhere, and k_j = D_j / (psi_j(x) - psi_j(x')) where (G dxi)_j != 0, any
k_j > 0 elsewhere. This is the sign condition of Mueller, Feliu, Regensburger, Conradi, Shiu & Dickenstein
(Found. Comput. Math. 2016) in the form of Shinar & Feinberg's concordance test (Math. Biosci. 2012); for
mass action it is equivalent to the Jacobian on S being nonsingular at every x and k (Craciun & Feinberg
2005, Feliu & Wiuf 2012), so injective networks also have only nondegenerate steady states.

Compared with the two-point sign MILP (two_point.certify), which asks for two steady states with the same k
anywhere, the difference is confined to one class (sigma in S with matching signs) and need not be a pair of
steady states (D in ker N is f(x) = f(x'), not f = 0); so neither the iterated route's equalities nor its
module classes are valid here.

A detached uptake (source species with free level) adds a conserved source per uptake; within one class the
source levels are equal, so injectivity with detached uptakes is the same question as with attached ones.

Signs are modelled with three binaries per metabolite and per reaction and indicator constraints (CPLEX);
"infeasible" proves injectivity up to the dynamic range big/eps of the sign thresholds; a witness is only
trusted after realise_witness rebuilds explicit x, x', k from it (tolerance artefacts fail that check).

    python injectivity.py <model> [time_limit_s] [threads]       -> json line on stdout
"""
import json, os, sys, time
import numpy as np, scipy.sparse as sp


def certify_injective(N, G, time_limit=3600, threads=4, eps=1.0, big=1e7):
    """'infeasible' = injective for every k > 0; otherwise returns a sign-compatible (sigma, D, dxi) if found.

    Returns (status, seconds, witness dict or None).
    """
    import cplex
    N = sp.csr_matrix(N, dtype=float); G = sp.csr_matrix(G, dtype=float); m, q = N.shape
    P = cplex.Cplex(); P.set_log_stream(None); P.set_results_stream(None); P.set_warning_stream(None)
    P.parameters.timelimit.set(time_limit); P.parameters.threads.set(threads); P.parameters.emphasis.mip.set(1)
    # eps = 1, not 1e-3: with substrate coefficients as small as 1e-4 (biomass-like complexes), G dxi at the
    # 1e-3 scale falls below CPLEX's feasibility tolerance and "zero" signs are accepted for nonzero products,
    # which produced spurious witnesses; numerical emphasis on for the same reason
    P.parameters.emphasis.numerical.set(1)
    z = [f"z{i}" for i in range(m)]; s = [f"s{i}" for i in range(m)]; w = [f"w{j}" for j in range(q)]
    d = [f"d{j}" for j in range(q)]; v = [f"v{j}" for j in range(q)]
    inf = cplex.infinity
    P.variables.add(names=z, lb=[-big] * m, ub=[big] * m)
    P.variables.add(names=s + w + d + v, lb=[-inf] * (m + 3 * q), ub=[inf] * (m + 3 * q))
    lin = P.linear_constraints
    # sigma = N w  (sigma in S),  N D = 0,  v = G dxi
    lin.add(lin_expr=[cplex.SparsePair([s[i]] + [w[j] for j in N.indices[N.indptr[i]:N.indptr[i + 1]]],
                                       [1.0] + list(-N.data[N.indptr[i]:N.indptr[i + 1]])) for i in range(m)],
            senses="E" * m, rhs=[0] * m)
    lin.add(lin_expr=[cplex.SparsePair([d[j] for j in N.indices[N.indptr[i]:N.indptr[i + 1]]],
                                       list(N.data[N.indptr[i]:N.indptr[i + 1]])) for i in range(m)],
            senses="E" * m, rhs=[0] * m)
    lin.add(lin_expr=[cplex.SparsePair([v[j]] + [z[i] for i in G.indices[G.indptr[j]:G.indptr[j + 1]]],
                                       [1.0] + list(-G.data[G.indptr[j]:G.indptr[j + 1]])) for j in range(q)],
            senses="E" * q, rhs=[0] * q)
    ind = P.indicator_constraints

    def signs(tag, n, small, large):
        """three binaries per index: small[i] and large[i] share a sign (or are both zero)"""
        p, ng, o = [f"{tag}p{i}" for i in range(n)], [f"{tag}n{i}" for i in range(n)], [f"{tag}o{i}" for i in range(n)]
        P.variables.add(names=p + ng + o, types="B" * (3 * n))
        lin.add(lin_expr=[cplex.SparsePair([p[i], ng[i], o[i]], [1, 1, 1]) for i in range(n)], senses="E" * n, rhs=[1] * n)
        for i in range(n):
            ind.add(indvar=p[i], lin_expr=cplex.SparsePair([small[i]], [1]), sense="G", rhs=eps)
            ind.add(indvar=p[i], lin_expr=cplex.SparsePair([large[i]], [1]), sense="G", rhs=1)
            ind.add(indvar=ng[i], lin_expr=cplex.SparsePair([small[i]], [1]), sense="L", rhs=-eps)
            ind.add(indvar=ng[i], lin_expr=cplex.SparsePair([large[i]], [1]), sense="L", rhs=-1)
            ind.add(indvar=o[i], lin_expr=cplex.SparsePair([small[i]], [1]), sense="E", rhs=0)
            ind.add(indvar=o[i], lin_expr=cplex.SparsePair([large[i]], [1]), sense="E", rhs=0)
        return p, ng
    mp, mn = signs("m", m, z, s)          # sign(dxi_i) = sign(sigma_i)
    signs("r", q, v, d)                   # sign((G dxi)_j) = sign(D_j)
    lin.add(lin_expr=[cplex.SparsePair(mp + mn, [1] * (2 * m))], senses="G", rhs=[1])     # dxi != 0
    t0 = time.time(); P.solve(); st = P.solution.get_status_string(); sec = time.time() - t0
    wit = None
    if P.solution.is_primal_feasible():
        wit = dict(dxi=np.array(P.solution.get_values(z)), sigma=np.array(P.solution.get_values(s)),
                   D=np.array(P.solution.get_values(d)), mp=np.array(P.solution.get_values(mp)), mn=np.array(P.solution.get_values(mn)))
    return st, sec, wit


def realise_witness(N, G, wit, tol=1e-6):
    """Turn a sign-compatible triple into explicit x != x' in one class and k > 0 with f_k(x) = f_k(x')
    (the sufficiency construction); returns the largest relative mismatch, or None if the signs do not hold."""
    N = sp.csr_matrix(N, dtype=float); G = sp.csr_matrix(G, dtype=float)
    z, sg, D = wit["dxi"], wit["sigma"], wit["D"]
    z = np.where(np.abs(z) < 1e-9, 0, z); sg = np.where(np.abs(sg) < 1e-9, 0, sg)
    z = z * (2 / np.abs(z).max())         # signs of dxi and G dxi are scale-free; keep exp() in range
    if np.any(np.sign(z) != np.sign(sg)):
        return None
    x2 = np.ones(len(z)); nz = z != 0
    x2[nz] = sg[nz] / np.expm1(z[nz]); x1 = x2 * np.exp(z)
    lp1, lp2 = G @ np.log(x1), G @ np.log(x2)
    dpsi = np.exp(lp1) - np.exp(lp2)
    k = np.ones(G.shape[0]); act = np.abs(G @ z) > 1e-9
    if np.any(np.sign(D[act]) != np.sign(dpsi[act])) or np.any(np.abs(D[~act]) > 1e-6):
        return None
    k[act] = D[act] / dpsi[act]
    f1, f2 = N @ (k * np.exp(lp1)), N @ (k * np.exp(lp2))
    return float(np.abs(f1 - f2).max() / max(np.abs(f1).max(), 1e-12)), float(np.abs(x1 - x2 - sg).max())


def jacobian_sign_test(A, Y, N):
    """Independent check for small networks (Feliu & Wiuf, Appl. Math. Comput. 2012): f_k is injective on every
    class for every k > 0 iff det(M) has nonzero coefficients of one sign only, as a polynomial in x and k,
    where M is the Jacobian of f_k with the pivot rows replaced by the reduced row echelon basis of S^perp
    (the conservation laws). Symbolic; small networks only."""
    import sympy as sy
    A, Y, N = (sy.Matrix(np.asarray(M.todense() if sp.issparse(M) else M).astype(int).tolist()) for M in (A, Y, N))
    m, q = N.shape
    xs = sy.symbols(f"x0:{m}", positive=True); ks = sy.symbols(f"k0:{q}", positive=True)
    psi = []
    for j in range(q):
        c = [i for i in range(A.shape[0]) if A[i, j] < 0]
        psi.append(sy.Mul(*[xs[s] ** Y[s, c[0]] for s in range(m)]) if c else sy.Integer(1))
    f = N * sy.Matrix([ks[j] * psi[j] for j in range(q)])
    J = f.jacobian(xs)
    Z = sy.Matrix.hstack(*N.T.nullspace()).T if N.rank() < m else sy.zeros(0, m)   # rows: conservation laws
    M = J.copy()
    if Z.shape[0]:
        # replace one row per conservation law: pick pivot columns of Z's rref
        R, piv = Z.rref()
        for r, p in enumerate(piv):
            M[p, :] = R[r, :]
    det = sy.Poly(sy.expand(M.det()), *xs, *ks)
    signs = {sy.sign(c) for c in det.coeffs() if c != 0}
    return len(signs) == 1 and 0 not in signs


if __name__ == "__main__":
    B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
    from run_kimonu_crn import triplets_to_sparse
    from kinetic_variants.two_point import substrate_composition
    M = sys.argv[1]; tl = float(sys.argv[2]) if len(sys.argv) > 2 else 3600; th = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    Nf = triplets_to_sparse(crn["S"], nm, nr)
    st, sec, wit = certify_injective(Nf, substrate_composition(A, Y), time_limit=tl, threads=th)
    out = dict(model=M, status=st, seconds=sec, m=nm, q=nr)
    if wit is not None:
        chk = realise_witness(Nf, substrate_composition(A, Y), wit)
        out["witness_realised"] = chk is not None and chk[0] < 1e-6
        out["witness_check"] = chk
    print(json.dumps(out), flush=True)
