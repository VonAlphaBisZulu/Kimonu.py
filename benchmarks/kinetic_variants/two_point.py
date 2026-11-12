"""Two positive steady states with the same rate constants: the exact criterion behind ACR, ACRR and coupling.

Let x, x' be positive steady states for rate constants k > 0, r = k o psi(x), r' = k o psi(x') their rates,
and dxi = ln x' - ln x. Mass action gives r'_j / r_j = exp((G dxi)_j), G = S^T Y^T (row j = composition of
reaction j's substrate complex; zero for uptakes, whose psi = 1). Both are steady states: N r = 0, N r' = 0.
Conversely, any r > 0 in ker N and dxi with N (r o exp(G dxi)) = 0 are two such steady states: take x = 1,
k = r, x' = exp(dxi). Hence, for a log-linear quantity c . ln x (c = e_s for ACR of s, e_i - e_j for a fixed
ratio, y_c - y_d for kinetic coupling of two complexes),

    c . ln x is the same at every positive steady state, for every k > 0
        <=>  no r > 0 and dxi with  N r = 0,  N (r o exp(G dxi)) = 0,  c . dxi != 0.

Nothing ties x and x' to one conservation class, as ACR requires. A detached uptake has its own source,
whose level is free, so its rate ratio is a free positive number rather than exp(0) = 1.

Three tools act on this criterion:
  realise   for a fixed dxi the criterion is an LP in r; feasible = an explicit counterexample.
  search    looks for such a dxi: nonlinear least squares in (ln r, dxi) from random starts, then realise.
  certify   the sign relaxation sign(r' - r) = sign(G dxi), as a MILP; infeasible = a proof. Valid
            equalities from the iterated route (V dxi = 0 for its fixed log-combinations) and one sign per
            kinetic module (equal G dxi inside a module) shrink it.
"""
import time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog, least_squares


def substrate_composition(A, Y):
    """G: reactions x metabolites, row j = composition of the substrate complex of reaction j."""
    return sp.csr_matrix(((sp.csr_matrix(A) < 0).astype(float).T @ sp.csr_matrix(Y, dtype=float).T))


def uptakes(A, Y):
    """Reactions whose substrate is the zero complex."""
    G = substrate_composition(A, Y)
    return np.flatnonzero(np.diff(G.indptr) == 0)


def realise(N, G, dxi, free=(), tol=1e-7):
    """An r > 0 with N r = 0 and N r' = 0, r'_j = r_j exp((G dxi)_j) (r'_j free for j in `free`), or None.

    Returns (r, r', relative residual). The ratios are scaled by their maximum, which N r' = 0 allows.
    """
    N = sp.csr_matrix(N, dtype=float); m, q = N.shape
    v = G @ dxi; lam = np.exp(v - v.max())
    free = np.asarray(sorted(free), dtype=int)
    # variables: r (q), then r'_free; r' on the other reactions is lam o r
    lam_c = lam.copy(); lam_c[free] = 0
    A_eq = sp.vstack([sp.hstack([N, sp.csr_matrix((m, len(free)))]),
                      sp.hstack([N @ sp.diags(lam_c), N[:, free]])]).tocsr()
    lp = linprog(np.zeros(q + len(free)), A_eq=A_eq, b_eq=np.zeros(2 * m), bounds=[(1, None)] * (q + len(free)),
                 method="highs")
    if lp.status != 0:
        return None
    r = lp.x[:q]; r2 = lam * r; r2[free] = lp.x[q:]
    res = max(np.abs(N @ r).max() / np.abs(r).max(), np.abs(N @ r2).max() / np.abs(r2).max())
    return (r, r2, res) if res < tol else None


def search(N, G, c, V=None, free=(), starts=20, scales=(0.5, 1, 2, 4), seed=0, time_limit=600, tol=1e-9, rho0=None):
    """Look for two steady states that differ in c . ln x. Returns (dxi, r, r', residual) or None.

    Unknowns rho = ln r (q) and dxi (m); residuals N e^rho, N e^(rho + G dxi), c . dxi - 1, mean(rho), and
    V dxi when V (rows = log-combinations already proven fixed) is given, since every pair must satisfy them.
    rho0: start ln r near this (e.g. ln of an interior flux) instead of a random point.
    Free (detached-uptake) reactions get their own log-ratio variable. A converged point is re-checked with
    `realise` at its dxi, independently of the least-squares residual.
    """
    N = sp.csr_matrix(N, dtype=float); G = sp.csr_matrix(G); m, q = N.shape
    free = np.asarray(sorted(free), dtype=int); nf = len(free)
    Fm = sp.csr_matrix((np.ones(nf), (free, np.arange(nf))), shape=(q, nf))   # free log-ratios enter like G
    Gx = sp.hstack([G, Fm]).tocsr()
    cx = np.concatenate([c, np.zeros(nf)])
    Vx = sp.hstack([V, sp.csr_matrix((V.shape[0], nf))]).tocsr() if V is not None and V.shape[0] else None
    rng = np.random.default_rng(seed); t0 = time.time()

    def unpack(u): return u[:q], u[q:]
    def fun(u):
        rho, w = unpack(u); r = np.exp(np.clip(rho, -50, 50)); r2 = np.exp(np.clip(rho + Gx @ w, -50, 50))
        out = [N @ r / r.max(), N @ r2 / r2.max(), [cx @ w - 1], [rho.mean()]]
        if Vx is not None: out.append(Vx @ w)
        return np.concatenate(out)
    def jac(u):
        rho, w = unpack(u); r = np.exp(np.clip(rho, -50, 50)); r2 = np.exp(np.clip(rho + Gx @ w, -50, 50))
        J1 = sp.hstack([N @ sp.diags(r / r.max()), sp.csr_matrix((m, m + nf))])
        Nr2 = N @ sp.diags(r2 / r2.max())
        J2 = sp.hstack([Nr2, Nr2 @ Gx])
        J3 = sp.hstack([sp.csr_matrix((1, q)), sp.csr_matrix(cx[None, :])])
        J4 = sp.hstack([sp.csr_matrix(np.full((1, q), 1 / q)), sp.csr_matrix((1, m + nf))])
        Js = [J1, J2, J3, J4]
        if Vx is not None: Js.append(sp.hstack([sp.csr_matrix((Vx.shape[0], q)), Vx]))
        return sp.vstack(Js).tocsr()

    best = None
    for k in range(starts):
        for s in scales:
            if time.time() - t0 > time_limit:
                return best
            w0 = rng.normal(size=m + nf); w0 *= s / max(abs(cx @ w0), 1e-9) * np.sign(cx @ w0 or 1)
            rstart = rng.normal(size=q) if rho0 is None else rho0 - rho0.mean() + 0.1 * rng.normal(size=q)
            u0 = np.concatenate([rstart, w0])
            sol = least_squares(fun, u0, jac=jac, method="trf", tr_solver="lsmr", max_nfev=200, xtol=1e-12, ftol=1e-14, gtol=1e-12)
            if np.abs(sol.fun).max() > 1e-6:
                continue
            w = sol.x[q:]; dxi = w[:m]
            if nf:   # free log-ratios: hand realise their reactions as free
                hit = realise(N, G, dxi, free=free)
            else:
                hit = realise(N, G, dxi)
            if hit is not None and abs(c @ dxi) > 1e-6:
                return dxi, hit[0], hit[1], hit[2]
    return best


def certify(N, G, c, V=None, classes=None, free=(), time_limit=3600, threads=4, eps=1e-3, big=1e4):
    """Sign relaxation as a MILP. 'infeasible' proves c . ln x fixed at every positive steady state.

    classes: list of reaction-index lists whose G dxi are equal for every admissible dxi (reactions whose
    substrates lie in one kinetic module); they share one sign. Only valid together with the V that makes
    them equal. Returns (status, seconds, dxi or None).
    """
    import cplex
    N = sp.csr_matrix(N, dtype=float); G = sp.csr_matrix(G); m, q = N.shape
    free = set(int(j) for j in free)
    if classes is None:
        classes = [[j] for j in range(q)]
    classes = [[j for j in cl if j not in free] for cl in classes]
    classes = [cl for cl in classes if cl]
    P = cplex.Cplex(); P.set_log_stream(None); P.set_results_stream(None); P.set_warning_stream(None)
    P.parameters.timelimit.set(time_limit); P.parameters.threads.set(threads); P.parameters.emphasis.mip.set(1)
    z = [f"z{i}" for i in range(m)]; d = [f"d{j}" for j in range(q)]
    P.variables.add(names=z, lb=[-big] * m, ub=[big] * m)
    P.variables.add(names=d, lb=[-cplex.infinity] * q, ub=[cplex.infinity] * q)
    P.linear_constraints.add(lin_expr=[cplex.SparsePair([d[j] for j in N.indices[N.indptr[i]:N.indptr[i + 1]]],
                                                        list(N.data[N.indptr[i]:N.indptr[i + 1]])) for i in range(m)],
                             senses="E" * m, rhs=[0] * m)                                    # N (r' - r) = 0
    if V is not None and V.shape[0]:
        V = sp.csr_matrix(V, dtype=float)
        P.linear_constraints.add(lin_expr=[cplex.SparsePair([z[i] for i in V.indices[V.indptr[k]:V.indptr[k + 1]]],
                                                            list(V.data[V.indptr[k]:V.indptr[k + 1]])) for k in range(V.shape[0])],
                                 senses="E" * V.shape[0], rhs=[0] * V.shape[0])
    nz = np.flatnonzero(c)
    P.linear_constraints.add(lin_expr=[cplex.SparsePair([z[i] for i in nz], list(map(float, c[nz])))], senses="G", rhs=[1])
    ind = P.indicator_constraints
    for k, cl in enumerate(classes):
        g = G[cl[0]]
        vk, pk, nk, ok = f"v{k}", f"p{k}", f"n{k}", f"o{k}"
        P.variables.add(names=[vk], lb=[-cplex.infinity], ub=[cplex.infinity])
        P.variables.add(names=[pk, nk, ok], types="BBB")
        P.linear_constraints.add(lin_expr=[cplex.SparsePair([vk] + [z[i] for i in g.indices], [1.0] + list(-g.data)),
                                           cplex.SparsePair([pk, nk, ok], [1, 1, 1])], senses="EE", rhs=[0, 1])
        ind.add(indvar=pk, lin_expr=cplex.SparsePair([vk], [1]), sense="G", rhs=eps)
        ind.add(indvar=nk, lin_expr=cplex.SparsePair([vk], [1]), sense="L", rhs=-eps)
        ind.add(indvar=ok, lin_expr=cplex.SparsePair([vk], [1]), sense="E", rhs=0)
        for j in cl:
            ind.add(indvar=pk, lin_expr=cplex.SparsePair([d[j]], [1]), sense="G", rhs=1)
            ind.add(indvar=nk, lin_expr=cplex.SparsePair([d[j]], [1]), sense="L", rhs=-1)
            ind.add(indvar=ok, lin_expr=cplex.SparsePair([d[j]], [1]), sense="E", rhs=0)
    t0 = time.time(); P.solve(); st = P.solution.get_status_string(); sec = time.time() - t0
    dxi = None
    if P.solution.is_primal_feasible():
        dxi = np.array(P.solution.get_values(z))
    return st, sec, dxi


def iterated_constraints(A, Y, N, detach=False):
    """V (rows: log-combinations the iterated route proves fixed) and reaction classes (one per module)."""
    import kimonu
    from kimonu.readout import composition_differences
    it = kimonu.iterated_route(A, Y, N, detach_uptakes=detach)
    D = composition_differences(sp.csc_matrix(Y), it.modules, it.acr)
    V = sp.csr_matrix(D.T, dtype=float)
    Ac = sp.csc_matrix(A)
    sub = {j: int(Ac.indices[Ac.indptr[j]:Ac.indptr[j + 1]][Ac.data[Ac.indptr[j]:Ac.indptr[j + 1]] < 0][0]) for j in range(A.shape[1])}
    cls = {}
    for j, s in sub.items():
        cls.setdefault(int(it.module_of[s]), []).append(j)
    return V, list(cls.values()), it


def exact_realise(N, G, xprime):
    """Exact counterexample for integer G and rational x' (x = 1): an r > 0 in ker N with
    N (r o psi(x')) = 0, from the exact nullspace of the stacked matrix. Returns r (Fractions) or None."""
    from fractions import Fraction
    from kimonu import interior_flux
    from kimonu.linalg import entries, exact
    triples, (m, q) = entries(N)
    Gc = sp.csr_matrix(G)
    psi = []
    for j in range(q):
        p = Fraction(1)
        for i, g in zip(Gc.indices[Gc.indptr[j]:Gc.indptr[j + 1]], Gc.data[Gc.indptr[j]:Gc.indptr[j + 1]]):
            assert float(g).is_integer(), "exact check needs integer substrate compositions"
            p *= Fraction(xprime[i]) ** int(g)
        psi.append(p)
    stacked = triples + [(m + i, j, v * psi[j]) for i, j, v in triples]
    try:
        return interior_flux(exact(stacked, (2 * m, q)), draw="lp")
    except ValueError:
        return None
