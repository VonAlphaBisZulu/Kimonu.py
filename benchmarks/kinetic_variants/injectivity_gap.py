"""The gap between the iterated and the tangent route as a subspace, and the two certificates that could close it.

Every log-linear quantity c . ln x that is fixed at all positive steady states (ACR: c = e_s, ACRR: e_i - e_j)
lies in one subspace C_true, and C_iter <= C_true <= C_tan: C_iter is spanned by the iterated route's fixed
log-combinations (rows of V: composition differences inside its modules, unit vectors of its ACR metabolites),
C_tan = T^perp is the row space of J = N diag(r*) S^T Y^T at a generic interior point. ACRR links are vectors
of the same subspaces, so they are proven or open exactly when the ACR question in their direction is.
T(r) depends on the interior flux r, so the tangent bound is C_tan* = (sum of T(r) over generic r)^perp, sampled
over seeds until the sum stops growing; one point alone (C_tan = T(r*)^perp) overstates it. The gap is
C_tan* / C_iter, of dimension g* = dim ker V - dim sum T(r). For the ACR/ACRR question only the span of the
tangent-only relation vectors matters: g_rel = the dimension they add to C_iter.

    python injectivity_gap.py <model> dims                  -> g, containment check, quotient functionals
    python injectivity_gap.py <model> inject  [s] [threads] -> injectivity MILP on the full network
    python injectivity_gap.py <model> quotient [s] [threads]-> one two-point sign MILP for the whole quotient

quotient: one two-point sign MILP for all tangent-only ACR/ACRR relations at once (their disjunction), with the
iterated equalities; infeasible proves all of them for generic k, a sign pattern is inconclusive.
"""
import json, os, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import straindesign as sd
from run_kimonu_crn import triplets_to_sparse, to_exact
from kimonu.linalg import entries, exact
from kinetic_variants.two_point import substrate_composition, iterated_constraints, realise
from kinetic_variants.tangent_fast import integer_positive_flux

OUT = os.environ.get("INJ_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_injectivity")


def load(M):
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    Nf = triplets_to_sparse(crn["S"], nm, nr)
    return crn, A, Y, Nf


def tangent_matrix(A, Y, N, seed=0):
    """Exact J = N diag(r*) G (m x m, attached uptakes), r* an integer interior flux; T = ker J."""
    flux = integer_positive_flux(N, seed=seed)
    G = substrate_composition(A, Y)
    Gx = to_exact(sp.csr_matrix(G))
    gtrip, _ = entries(Gx)
    byj = {}
    for j, i, v in gtrip: byj.setdefault(j, []).append((i, v))
    ntrip, (m, q) = entries(N)
    J = {}
    for i, j, v in ntrip:
        for k, w in byj.get(j, ()):
            J[(i, k)] = J.get((i, k), Fraction(0)) + flux[j] * v * w
    return exact(((i, k, x) for (i, k), x in J.items()), (m, m))


def rows_of(Mx):
    triples, shape = entries(Mx)
    out = {}
    for i, j, v in triples: out.setdefault(i, {})[j] = v
    return out, shape


def exact_V(Y_triplets, nm, modules, known_acr):
    """V, the iterated route's fixed log-combinations, built in exact arithmetic: composition differences inside
    each module (first member against the others) and unit vectors of its ACR metabolites. The float version
    (kimonu.acr.composition_differences on a float Y) leaves rounding residues such as 0.0003 - 0.0001 =
    0.00019999999999999998 in biomass-like rows, which an exact rank then counts as real."""
    col = {}
    for i, c, v in Y_triplets:                      # crn.json triplets are 1-based
        col.setdefault(int(c) - 1, {})[int(i) - 1] = Fraction(str(v))
    rows = []
    for mod in modules:
        mem = sorted(mod)
        for o in mem[1:]:
            d = dict(col.get(mem[0], {}))
            for i, v in col.get(o, {}).items(): d[i] = d.get(i, 0) - v
            d = {i: v for i, v in d.items() if v != 0}
            if d: rows.append(d)
    for s_ in sorted(known_acr):
        rows.append({int(s_): Fraction(1)})
    trip = [(k, i, v) for k, d in enumerate(rows) for i, v in d.items()]
    return sd.RationalMatrix.from_fractions(trip, (len(rows), nm)), rows


def kernel_dim(rowsets, m):
    """dim of the common kernel of several exact matrices given as lists of {col: value} rows"""
    trip, off = [], 0
    for rows in rowsets:
        trip += [(off + k, j, v) for k, d in enumerate(rows) for j, v in d.items()]; off += len(rows)
    if off == 0:
        return m
    return sd.sparse_nullspace(sd.RationalMatrix.from_fractions(trip, (off, m))).shape[1]


def matrix_rows(Mx):
    rows = {}
    for i, j, v in entries(Mx)[0]: rows.setdefault(i, {})[j] = v
    return [rows[i] for i in sorted(rows)]


def column_rows(Mx):
    """basis vectors (columns of a nullspace matrix) as {row index: value} dicts"""
    cols = {}
    for i, j, v in entries(Mx)[0]: cols.setdefault(j, {})[i] = v
    return [cols[j] for j in sorted(cols)]


def gap_dims(crn, A, Y, N, Nf, relations=(), seeds=4, log=print):
    """C_iter <= C_true <= C_tan*, with C_tan* = (sum over generic r of T(r))^perp, sampled over seeds until the
    sum stops growing. g* = dim C_tan* - dim C_iter; g_rel = dimension that the given relation vectors (ACR e_s,
    ACRR e_i - e_j) add to C_iter. All ranks exact."""
    import kimonu
    m = Nf.shape[0]; t0 = time.time()
    it = kimonu.kinetic_modules(A, Y, N)
    Vx, Vrows = exact_V(crn["Y"], m, it.modules, it.known_acr)
    dKV = kernel_dim([Vrows], m); log(f"iterated route {time.time() - t0:.0f}s, V {len(Vrows)} rows, dim ker V {dKV}")
    # T(r) for several generic r. If every later T lies in the first one (J_0 T_k = 0, an exact product), the
    # tangent space does not depend on r and C_tan* = T(r_0)^perp. Otherwise dim(T_0 + T_k) = dT_0 + dT_k -
    # dim ker [J_0; J_k] (sparse) gives a lower bound on dim sum T(r), i.e. an upper bound on g*.
    J0 = tangent_matrix(A, Y, N, seed=0); J0rows = matrix_rows(J0); T0 = sd.sparse_nullspace(J0); dT0 = T0.shape[1]
    J0cols = {}
    for i, row in enumerate(J0rows):
        for j, v in row.items(): J0cols.setdefault(j, []).append((i, v))
    dims_T, contained, in_T0, grow = [dT0], [kernel_dim([Vrows, J0rows], m) == dT0], [True], [dT0]
    log(f"seed 0: dim T {dT0}, T <= ker V {contained[-1]} ({time.time() - t0:.0f}s)")
    for seed in range(1, seeds):
        J = tangent_matrix(A, Y, N, seed=seed); Jrows = matrix_rows(J); T = sd.sparse_nullspace(J); dT = T.shape[1]
        dims_T.append(dT); contained.append(kernel_dim([Vrows, Jrows], m) == dT)
        sub = True
        for b in column_rows(T):
            acc = {}
            for j, v in b.items():
                for i, w in J0cols.get(j, ()): acc[i] = acc.get(i, 0) + w * v
            if any(x != 0 for x in acc.values()):
                sub = False; break
        in_T0.append(sub)
        grow.append(dT0 if sub else dT0 + dT - kernel_dim([J0rows, Jrows], m))
        log(f"seed {seed}: dim T {dT}, T <= ker V {contained[-1]}, T <= T_0 {sub}, dim(T_0 + T) {grow[-1]} ({time.time() - t0:.0f}s)")
    out = dict(m=m, dim_kerV=dKV, dim_C_iter=m - dKV, dims_T=dims_T, T_in_kerV=contained, T_in_T0=in_T0, dim_T0_plus_Tk=grow,
               dim_C_tan_single=m - dims_T[0], dim_C_tan_star_upper=m - max(grow), g_single=dKV - dims_T[0], g_star_upper=dKV - max(grow),
               g_star_exact=all(in_T0))
    if relations:
        mets = crn["metabolites"]
        C = [{mets.index(k): Fraction(v) for k, v in r["c"].items()} for r in relations]
        out["n_relations"] = len(C)
        out["g_rel"] = dKV - kernel_dim([Vrows, C], m)
        out["relations_in_C_iter"] = [kernel_dim([Vrows, [c]], m) == dKV for c in C]
    out["seconds"] = time.time() - t0
    return out, Vrows, it


def tangent_sum(M, max_seeds=20, primes=(4611686018427387847, 2305843009213693951)):
    """dim of sum T(r) over generic interior fluxes, by modular rank (python-flint nmod_mat) of the stacked exact
    bases, adding seeds until the rank stops growing or reaches its ceiling dim ker [V; relations] (every
    tangent-only relation is orthogonal to every generic T(r)). A rank mod p never exceeds the rank over Q, so
    g*_upper = dim ker V - rank is an upper bound on g*; at the ceiling it equals g_rel and is exact."""
    import flint
    crn, A, Y, Nf = load(M); N = to_exact(Nf); m = Nf.shape[0]; mets = crn["metabolites"]; t0 = time.time()
    import kimonu
    it = kimonu.kinetic_modules(A, Y, N)
    _, Vrows = exact_V(crn["Y"], m, it.modules, it.known_acr)
    relf = f"/scratch/sadreyev_lab/pschneider/kimonu_twopoint/{M}_relations.json"
    C = [{mets.index(k): Fraction(v) for k, v in r["c"].items()} for r in json.load(open(relf))["relations"]]
    dKV = kernel_dim([Vrows], m); ceiling = kernel_dim([Vrows, C], m)
    print(f"[{M}] dim ker V {dKV}, ceiling dim ker [V; relations] {ceiling}", flush=True)
    rows = []; ranks = []
    for seed in range(max_seeds):
        T = sd.sparse_nullspace(tangent_matrix(A, Y, N, seed=seed))
        rows += column_rows(T)
        rk = []
        for p in primes:
            Mp = flint.nmod_mat(len(rows), m, p)
            for k, d in enumerate(rows):
                for j, v in d.items():
                    Mp[k, j] = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
            rk.append(Mp.rank())
        ranks.append(max(rk))
        print(f"[{M}] seed {seed}: rank of stacked tangent bases {ranks[-1]} (primes {rk}) ({time.time() - t0:.0f}s)", flush=True)
        if ranks[-1] == ceiling or (len(ranks) >= 3 and ranks[-1] == ranks[-2] == ranks[-3]):
            break
    out = dict(model=M, m=m, dim_kerV=dKV, ceiling=ceiling, g_rel=dKV - ceiling, tsum_ranks=ranks,
               g_star_upper=dKV - ranks[-1], g_star_exact=ranks[-1] == ceiling, seconds=time.time() - t0)
    json.dump(out, open(f"{OUT}/{M}_tsum.json", "w"))
    print(json.dumps(out), flush=True)


def dims(M):
    crn, A, Y, Nf = load(M); N = to_exact(Nf)
    relf = f"/scratch/sadreyev_lab/pschneider/kimonu_twopoint/{M}_relations.json"
    rel = json.load(open(relf))["relations"] if os.path.exists(relf) else []
    out, Vrows, it = gap_dims(crn, A, Y, N, Nf, relations=rel, log=lambda s: print(f"[{M}] {s}", flush=True))
    out["model"] = M
    os.makedirs(OUT, exist_ok=True)
    m = Nf.shape[0]
    sp.save_npz(f"{OUT}/{M}_V.npz", sp.csr_matrix(([float(v) for d in Vrows for v in d.values()],
                ([k for k, d in enumerate(Vrows) for _ in d], [i for d in Vrows for i in d])), shape=(len(Vrows), m)))
    json.dump(out, open(f"{OUT}/{M}_dims.json", "w"))
    print(json.dumps(out), flush=True)


def inject(M, tl=7200, threads=8):
    from kinetic_variants.injectivity import certify_injective, realise_witness
    crn, A, Y, Nf = load(M); G = substrate_composition(A, Y)
    st, sec, wit = certify_injective(Nf, G, time_limit=tl, threads=threads)
    out = dict(model=M, status=st, seconds=sec)
    if wit is not None:
        chk = realise_witness(Nf, G, wit)
        out.update(witness_found=True, witness_check=chk)
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(f"{OUT}/{M}_inject.json", "w"))
    print(json.dumps(out), flush=True)


def quotient(M, tl=10800, threads=8, eps=1e-3, big=1e4):
    """One two-point sign MILP for all tangent-only relations at once: a pair of steady states with equal k whose
    log-difference satisfies the iterated equalities (V dxi = 0) and moves at least one relation (c_i . dxi >= 1
    or <= -1 for some i). Infeasible proves every relation (generic k); a sign pattern is inconclusive and is
    tried as a counterexample along its ray."""
    import cplex
    crn, A, Y, Nf = load(M); G = substrate_composition(A, Y); mets = crn["metabolites"]
    rel = json.load(open(f"/scratch/sadreyev_lab/pschneider/kimonu_twopoint/{M}_relations.json"))
    classes = rel["classes"]; C = [{mets.index(k): float(v) for k, v in r["c"].items()} for r in rel["relations"]]
    V = sp.load_npz(f"{OUT}/{M}_V.npz")
    N = sp.csr_matrix(Nf, dtype=float); m, q = N.shape; g = len(C)
    P = cplex.Cplex(); P.set_log_stream(None); P.set_results_stream(None); P.set_warning_stream(None)
    P.parameters.timelimit.set(tl); P.parameters.threads.set(threads); P.parameters.emphasis.mip.set(1)
    z = [f"z{i}" for i in range(m)]; d = [f"d{j}" for j in range(q)]
    P.variables.add(names=z, lb=[-big] * m, ub=[big] * m)
    P.variables.add(names=d, lb=[-cplex.infinity] * q, ub=[cplex.infinity] * q)
    lin = P.linear_constraints; ind = P.indicator_constraints
    lin.add(lin_expr=[cplex.SparsePair([d[j] for j in N.indices[N.indptr[i]:N.indptr[i + 1]]], list(N.data[N.indptr[i]:N.indptr[i + 1]])) for i in range(m)],
            senses="E" * m, rhs=[0] * m)
    V = sp.csr_matrix(V, dtype=float)
    lin.add(lin_expr=[cplex.SparsePair([z[i] for i in V.indices[V.indptr[k]:V.indptr[k + 1]]], list(V.data[V.indptr[k]:V.indptr[k + 1]])) for k in range(V.shape[0])],
            senses="E" * V.shape[0], rhs=[0] * V.shape[0])
    hp, hn = [f"hp{k}" for k in range(g)], [f"hn{k}" for k in range(g)]
    P.variables.add(names=hp + hn, types="B" * (2 * g))
    for k, c in enumerate(C):
        e = cplex.SparsePair([z[i] for i in c], list(c.values()))
        ind.add(indvar=hp[k], lin_expr=e, sense="G", rhs=1)
        ind.add(indvar=hn[k], lin_expr=e, sense="L", rhs=-1)
    lin.add(lin_expr=[cplex.SparsePair(hp + hn, [1] * (2 * g))], senses="G", rhs=[1])          # some relation moves
    Gc = sp.csr_matrix(G)
    for k, cl in enumerate(classes):
        gr = Gc[cl[0]]
        vk, pk, nk, ok = f"v{k}", f"p{k}", f"n{k}", f"o{k}"
        P.variables.add(names=[vk], lb=[-cplex.infinity], ub=[cplex.infinity]); P.variables.add(names=[pk, nk, ok], types="BBB")
        lin.add(lin_expr=[cplex.SparsePair([vk] + [z[i] for i in gr.indices], [1.0] + list(-gr.data)),
                          cplex.SparsePair([pk, nk, ok], [1, 1, 1])], senses="EE", rhs=[0, 1])
        ind.add(indvar=pk, lin_expr=cplex.SparsePair([vk], [1]), sense="G", rhs=eps)
        ind.add(indvar=nk, lin_expr=cplex.SparsePair([vk], [1]), sense="L", rhs=-eps)
        ind.add(indvar=ok, lin_expr=cplex.SparsePair([vk], [1]), sense="E", rhs=0)
        for j in cl:
            ind.add(indvar=pk, lin_expr=cplex.SparsePair([d[j]], [1]), sense="G", rhs=1)
            ind.add(indvar=nk, lin_expr=cplex.SparsePair([d[j]], [1]), sense="L", rhs=-1)
            ind.add(indvar=ok, lin_expr=cplex.SparsePair([d[j]], [1]), sense="E", rhs=0)
    t0 = time.time(); P.solve(); st = P.solution.get_status_string(); sec = time.time() - t0
    out = dict(model=M, relations=g, status=st, seconds=sec)
    if P.solution.is_primal_feasible():
        dxi = np.array(P.solution.get_values(z))
        out["moved"] = [k for k, c in enumerate(C) if abs(sum(v * dxi[i] for i, v in c.items())) > 0.5]
        out["sign_pattern_realised"] = any(realise(N, G, t * dxi) is not None for t in (0.001, 0.01, 0.1, 1))
    json.dump(out, open(f"{OUT}/{M}_quotient.json", "w"))
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    M, cmd = sys.argv[1], sys.argv[2]; rest = [float(x) for x in sys.argv[3:]]
    if cmd == "dims": dims(M)
    elif cmd == "tsum": tangent_sum(M)
    elif cmd == "inject": inject(M, *(rest[:1] or [7200]), *([int(rest[1])] if len(rest) > 1 else []))
    elif cmd == "quotient": quotient(M, *(rest[:1] or [10800]), *([int(rest[1])] if len(rest) > 1 else []))
