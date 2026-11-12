"""The sign certificate on a local subsystem, grown by lazy row generation.

certify (two_point.py) asks for a sign pattern of the whole network: d = r' - r in ker N with sign(d_j) =
sign((G dxi)_j), V dxi = 0 and c . dxi >= 1. Infeasible proves c . ln x fixed (generic k, up to the dynamic
range of the sign thresholds eps and big). On the gap models this MILP times out after 3 h per relation.

Dropping balance rows only relaxes it: a subset R of the rows of N, with d and its sign binaries kept only
for the reactions that touch R, is a weaker system, so if it is already infeasible, the full one is too and
the relation is proven. All of dxi, all of V and the c . dxi >= 1 constraint stay. Reactions of one kinetic
module share one sign, as in certify; a module is cut down to its reactions that touch R.

Per relation, R starts at the metabolites of c. Each round solves the local MILP. Infeasible: proven.
Feasible: the class sign pattern is checked by an LP on the full N (d_j >= 1, <= -1 or = 0 on the signed
reactions, free elsewhere), recorded for information only, and R grows by the metabolites of every reaction
with d_j != 0 in the solution; the 20 metabolites of highest degree (hubs) enter only when nothing else is
left, and when the solution touches nothing new, R grows by the neighbours of the reactions it touches. The
loop stops after 8 rounds or when R exceeds 1500 rows.

    python kinetic_variants/two_point_local.py <model> <relation index> [time limit per round, s]
        -> $TWO_POINT_LOCAL_OUT/<model>_rel<i>.json    (relations, V and classes from two_point_gap.py prepare)
    python kinetic_variants/two_point_local.py toys     soundness check on networks with a known answer
"""
import json, os, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]

PREP = os.environ.get("TWO_POINT_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_twopoint")
OUT = os.environ.get("TWO_POINT_LOCAL_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_twopoint_local")
HUBS, MAX_ROUNDS, MAX_ROWS = 20, 8, 1500


def local_certify(N, G, c, rows, V=None, classes=None, time_limit=900, threads=8, eps=1e-3, big=1e4):
    """certify with only the balance rows in `rows`. Returns (status, seconds, signs, d, stats): signs maps each
    signed reaction to -1, 0 or +1 and d holds its value in the solution (both None unless feasible)."""
    import cplex
    N = sp.csr_matrix(N, dtype=float); G = sp.csr_matrix(G); m, q = N.shape
    NR = N[sorted(rows)]
    touch = sorted(set(NR.indices.tolist()))                          # reactions that touch R
    if classes is None:
        classes = [[j] for j in range(q)]
    tset = set(touch)
    local = [[j for j in cl if j in tset] for cl in classes]
    local = [cl for cl in local if cl]
    P = cplex.Cplex(); P.set_log_stream(None); P.set_results_stream(None); P.set_warning_stream(None)
    P.parameters.timelimit.set(time_limit); P.parameters.threads.set(threads); P.parameters.emphasis.mip.set(1)
    z = [f"z{i}" for i in range(m)]; d = {j: f"d{j}" for j in touch}
    P.variables.add(names=z, lb=[-big] * m, ub=[big] * m)
    P.variables.add(names=[d[j] for j in touch], lb=[-cplex.infinity] * len(touch), ub=[cplex.infinity] * len(touch))
    P.linear_constraints.add(lin_expr=[cplex.SparsePair([d[j] for j in NR.indices[NR.indptr[i]:NR.indptr[i + 1]]],
                                                        list(NR.data[NR.indptr[i]:NR.indptr[i + 1]])) for i in range(NR.shape[0])],
                             senses="E" * NR.shape[0], rhs=[0] * NR.shape[0])            # N_R (r' - r) = 0
    if V is not None and V.shape[0]:
        V = sp.csr_matrix(V, dtype=float)
        P.linear_constraints.add(lin_expr=[cplex.SparsePair([z[i] for i in V.indices[V.indptr[k]:V.indptr[k + 1]]],
                                                            list(V.data[V.indptr[k]:V.indptr[k + 1]])) for k in range(V.shape[0])],
                                 senses="E" * V.shape[0], rhs=[0] * V.shape[0])
    nz = np.flatnonzero(c)
    P.linear_constraints.add(lin_expr=[cplex.SparsePair([z[i] for i in nz], list(map(float, c[nz])))], senses="G", rhs=[1])
    ind = P.indicator_constraints
    for k, cl in enumerate(local):
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
    stats = dict(rows=len(rows), reactions=len(touch), classes=len(local))
    t0 = time.time(); P.solve(); st = P.solution.get_status_string(); sec = time.time() - t0
    signs = dv = None
    if P.solution.is_primal_feasible():
        vals = dict(zip([d[j] for j in touch], P.solution.get_values([d[j] for j in touch])))
        dv = {j: vals[d[j]] for j in touch}
        signs = {j: (1 if dv[j] > 0.5 else -1 if dv[j] < -0.5 else 0) for j in touch}
    return st, sec, signs, dv, stats


def full_lp(N, signs):
    """Does the local sign pattern extend to N d = 0 on the whole network (other reactions free)? For information."""
    N = sp.csr_matrix(N, dtype=float); m, q = N.shape
    bounds = [(None, None)] * q
    for j, s in signs.items():
        bounds[j] = (1, None) if s > 0 else (None, -1) if s < 0 else (0, 0)
    lp = linprog(np.zeros(q), A_eq=N, b_eq=np.zeros(m), bounds=bounds, method="highs")
    return lp.status == 0


def grow(N, rows, signs, hubs):
    """New rows: metabolites of the reactions that move (d != 0), non-hubs first; else neighbours of R."""
    Nc = sp.csc_matrix(N)
    moving = [j for j, s in signs.items() if s != 0]
    cand = set()
    for j in moving:
        cand |= set(Nc.indices[Nc.indptr[j]:Nc.indptr[j + 1]].tolist())
    cand -= set(rows)
    if not cand:                                                      # nothing new moves: widen by neighbours
        for j in signs:
            cand |= set(Nc.indices[Nc.indptr[j]:Nc.indptr[j + 1]].tolist())
        cand -= set(rows)
    plain = cand - hubs
    return plain if plain else cand


def decide(N, G, c, V=None, classes=None, time_limit=900, threads=8, log=print):
    """the local loop for one relation: ('proven' | 'feasible' | 'timeout' | 'closed', rounds)"""
    N = sp.csr_matrix(N, dtype=float)
    deg = np.diff(N.indptr)
    hubs = set(np.argsort(-deg)[:HUBS].tolist())
    rows = set(np.flatnonzero(c).tolist())
    rounds = []
    for rnd in range(MAX_ROUNDS):
        st, sec, signs, dv, stats = local_certify(N, G, c, rows, V=V, classes=classes, time_limit=time_limit, threads=threads)
        rec = dict(round=rnd, status=st, seconds=sec, **stats)
        if "infeasible" in st:
            rounds.append(rec); log(json.dumps(rec))
            return "proven", rounds
        if signs is None:
            rounds.append(rec); log(json.dumps(rec))
            return "timeout", rounds
        rec["moving"] = sum(s != 0 for s in signs.values())
        rec["extends_to_full_N"] = full_lp(N, signs)
        new = grow(N, rows, signs, hubs)
        rec["new_rows"] = len(new)
        rounds.append(rec); log(json.dumps(rec))
        if not new:
            return "closed", rounds                                  # R holds every row the pattern reaches
        rows |= new
        if len(rows) > MAX_ROWS:
            return "feasible", rounds
    return "feasible", rounds


def run_relation(M, i, time_limit=900):
    from kinetic_variants.two_point_gap import load
    from kinetic_variants.two_point import substrate_composition
    crn, A, Y, Nf = load(M); mets = crn["metabolites"]
    info = json.load(open(f"{PREP}/{M}_relations.json")); r = info["relations"][i]
    V = sp.load_npz(f"{PREP}/{M}_V.npz")
    G = substrate_composition(A, Y)
    c = np.zeros(len(mets))
    for k, v in r["c"].items():
        c[mets.index(k)] = v
    t0 = time.time()
    verdict, rounds = decide(Nf, G, c, V=V, classes=info["classes"], time_limit=time_limit,
                             log=lambda s: print(f"[{M} rel {i}] {s}", flush=True))
    out = dict(model=M, index=i, relation=r, verdict=verdict, rounds=rounds, seconds=time.time() - t0,
               max_rounds=MAX_ROUNDS, max_rows=MAX_ROWS, time_limit_per_round=time_limit)
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(f"{OUT}/{M}_rel{i}.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "rounds"}), flush=True)


def toys():
    """Known answers: Figure 2's kinase (X proven; C and P have a second steady state), the completeness example
    (A, B proven; C, P refuted) and Figure 1 (G proven with the iterated equalities only)."""
    sys.path.insert(0, os.path.join(os.path.dirname(B), "tests"))
    from networks import matrices
    from kinetic_variants.two_point import substrate_composition, iterated_constraints
    cases = [
        ("kinase", [("0", "X"), ("X", "0"), ("X+X", "C"), ("C", "0"), ("C+P+P", "P+P+P"), ("P", "0")],
         ["X", "C", "P"], {"X": "proven", "C": "not", "P": "not"}, False),
        ("completeness", [("0", "A"), ("A", "0"), ("A+A", "0"), ("0", "B"), ("B", "0"), ("A+B", "0"),
                          ("0", "C"), ("C", "0"), ("P", "0"), ("C+P+P", "P+P+P")],
         ["A", "B", "C", "P"], {"A": "proven", "B": "proven", "C": "not", "P": "not"}, False),
        ("figure1", [("A", "B"), ("B", "A"), ("B", "C"), ("C", "B"), ("C", "D"), ("D+E", "F"), ("F", "D+E"), ("F", "B+G"),
                     ("C+G", "H"), ("H", "C+G"), ("H", "C+E"), ("A+G", "I"), ("I", "A+G"), ("I", "A+E")],
         list("ABCDEFGHI"), {"G": "proven", "D": "not", "E": "not"}, True),
    ]
    ok = True
    for name, R, S, want, use_V in cases:
        A, Y, N, _ = matrices(R, S)
        G = substrate_composition(A, Y); Nf = sp.csr_matrix(N, dtype=float)
        V = classes = None
        if use_V:
            V, classes, _ = iterated_constraints(A, Y, N)
        for s, w in want.items():
            c = np.eye(len(S))[S.index(s)]
            verdict, rounds = decide(Nf, G, c, V=V, classes=classes, time_limit=60, threads=2, log=lambda _: None)
            good = (verdict == "proven") == (w == "proven")
            ok &= good
            print(f"{name:13s} {s}: {verdict:8s} rows {rounds[-1]['rows']} (want {w}) {'ok' if good else 'WRONG'}")
    print("all toys as expected" if ok else "TOY MISMATCH")
    return ok


if __name__ == "__main__":
    if sys.argv[1] == "toys":
        sys.exit(0 if toys() else 1)
    run_relation(sys.argv[1], int(sys.argv[2]), *(float(x) for x in sys.argv[3:4]))
