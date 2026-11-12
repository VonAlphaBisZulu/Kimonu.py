"""Local ACR as published (symbolic) against the tangent route (one exact nullspace), with ground truth.

Pascual-Escudero & Feliu (2022, Math. Methods Appl. Sci. 45:359-382, Theorem 5.3) decide local ACR in x_i
for every k: write the steady-state system as g_k(x) = N diag(k) x^B with N of full row rank s (linear
dependencies removed) and B the exponent (substrate composition) matrix; then x_i has local ACR iff all
s x s minors of N diag(v) B^T without column i vanish for every v in ker N. With v = W a for a kernel basis
W, the minors are polynomials in a, and the test is whether they vanish identically. Two symbolic versions
of this criterion are implemented here; both are ours, written from the paper, which gives no code:

  minors  the criterion as stated: enumerate the s x s minors of (N diag(v) B^T)_{-i} and expand each as a
          polynomial in a (early exit at the first nonzero one). ACR only.
  rref    the same rank condition by symbolic Gaussian elimination over the field Q(a): e_i (or e_i - e_j for
          a ratio) lies in the row space of N diag(v) B^T over Q(a). ACR and ACRR.
  kimonu  the tangent route: one exact sparse nullspace T = ker(N diag(r*) S^T Y^T) at a random integer
          interior flux r* (tangent_fast); zero rows of T are ACR, equal rows ACRR. A polynomial that is not
          identically zero is nonzero at a random point with high probability (Schwartz-Zippel), so this
          decides the same identity.

With N diag(v) B^T of rank < s for all v (degenerate steady states) every s x s minor vanishes and the
published criterion reports every species; such networks are flagged ("degenerate") instead of compared.
The iterated route (kimonu.iterated_route) is run alongside as the proven lower bound.

    python local_acr_symbolic.py small              -> ground-truth suite, all methods
    python local_acr_symbolic.py ladder <name> [s]  -> one network of the size ladder, time limit s per method
"""
import json, os, signal, sys, time
from fractions import Fraction
from itertools import combinations
import numpy as np, scipy.sparse as sp
B_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B_DIR, os.path.join(os.path.dirname(B_DIR), "src"), os.path.join(os.path.dirname(B_DIR), "tests")]
import sympy
from sympy.polys.matrices import DomainMatrix
import straindesign as sd
import kimonu
from kimonu.linalg import entries
from networks import matrices, ENVZ, ENVZ_SPECIES
from kinetic_variants.tangent_fast import one_shot_fast

OUT = os.environ.get("LACR_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_local_acr")


# ---- networks: (reactions, species, truth) -- truth: species -> True/False for unconditional ACR, where known
def _rev(a, b): return [(a, b), (b, a)]
NETWORKS = {
    # García Puente et al. 2025, Example 2.4/2.18: ACR in A, not in B, for every k
    "GP2.18": ([("B", "A"), ("A+A+B", "A+B+B")], ["A", "B"], {"A": True, "B": False}),
    # Example 2.7/2.19: capacity for ACR, but two positive steady states for k = (2, 1.5, 1): no unconditional ACR
    "GP2.19": ([("A+A", "A"), ("A+A+A", "A+A+A+A+A"), ("A+A+A+A", "A+A+A")], ["A"], {"A": False}),
    # Example 3.9 (Joshi & Nguyen bifunctional enzyme): ACR in S3 for every k
    "GP3.9": (_rev("S1+E", "C1") + [("C1", "S2+E")] + _rev("S2+E", "C2") + [("C2", "S3+E")]
              + _rev("S2+C1", "C3") + [("C3", "S3+C1")] + _rev("S3+C1", "C4") + [("C4", "S1+C1")],
              ["S1", "S2", "S3", "E", "C1", "C2", "C3", "C4"], {"S3": True}),
    # Example 3.12: for its special k = (1, 1, 1, 1, 1/2) the steady states form the line x_A = x_B (no ACR, so no
    # unconditional ACR); for generic k admitting a positive steady state it is unique, x_B = c x_A with
    # c = sqrt(k1 / 2k5) and x_A = (k4 - k3 c) / (k1 - k2 c) (derived here, checked on sampled k): generic ACR in both
    "GP3.12": ([("A+A", "A+A+A+B")] + _rev("A+B", "B") + [("A", "0"), ("B+B", "0")], ["A", "B"], {"A": True, "B": True}),
    # Example 3.14: two components x_A = 1 and x_A = 2 for its k: no ACR in A (nor B)
    "GP3.14": ([("A+A+B", "A+A+A"), ("A+B", "B+B"), ("B", "A")], ["A", "B"], {"A": False, "B": False}),
    # Example 3.15 (Shinar-Feinberg EnvZ-OmpR, 7 species): ACR in Yp for every k
    "GP3.15": (_rev("X", "XT") + [("XT", "Xp")] + _rev("Xp+Y", "XpY") + [("XpY", "X+Yp")]
               + _rev("XT+Yp", "XTYp") + [("XTYp", "XT+Y")], ["X", "XT", "Xp", "Y", "Yp", "XpY", "XTYp"], {"Yp": True}),
    # Example 3.17: stated for one k; for every k x_A^2 = k1/(2k2) and x_B = (k4 x_A - k5)/k3 are unique (derived here)
    "GP3.17": ([("0", "A"), ("A+A", "0"), ("B+B", "B+B+B"), ("A+B", "A"), ("B", "B+B")], ["A", "B"],
               {"A": True, "B": True}),
    # Example B.3: ACR in A (x_A = k2/k1), up to three steady-state values of B
    "GPB.3": (_rev("A+B", "B") + _rev("B", "0") + _rev("B+B", "B+B+B"), ["A", "B"], {"A": True, "B": False}),
    # Pascual-Escudero & Feliu 2022, Example 4.4: ACR in x1 for every k
    "PEF4.4": ([("X1+X2", "X2+X2"), ("X2", "X1")], ["X1", "X2"], {"X1": True, "X2": False}),
    # Example 4.6: two half-lines for D(k) > 0: local ACR in x1, but not ACR
    "PEF4.6": ([("X1+X1+X1+X2", "X1+X1+X1+X1"), ("X1+X1+X2", "X2+X2+X2"), ("X1+X2", "X1+X1")], ["X1", "X2"],
               {"X1": False, "X2": False}),
    # Example 5.5 (IDHKP-IDH core, Shinar & Feinberg 2010): ACR in x4
    "PEF5.5": (_rev("X1+X2", "X3") + [("X3", "X1+X4")] + _rev("X3+X4", "X5") + [("X5", "X3+X2")],
               ["X1", "X2", "X3", "X4", "X5"], {"X4": True}),
    # ours: Figure 1 (Langary's EnvZ-OmpR); truth from the iterated route (proof) and two-point counterexamples
    "Fig1": (ENVZ, ENVZ_SPECIES, {s: s == "Yp" for s in ENVZ_SPECIES}),
    # former completeness toy (dimers plus a switch): A, B have ACR (unique steady state of their subsystem); C, P take two values
    "Fig3": ([("0", "A"), ("A", "0"), ("A+A", "0"), ("0", "B"), ("B", "0"), ("A+B", "0"), ("0", "C"), ("C", "0"),
              ("P", "0"), ("C+P+P", "P+P+P")], ["A", "B", "C", "P"], {"A": True, "B": True, "C": False, "P": False}),
    "dimer": ([("0", "X"), ("X", "0"), ("X+X", "0")], ["X"], {"X": True}),
    "heterodimer": ([("0", "X"), ("X", "0"), ("X+X", "0"), ("0", "Y"), ("Y", "0"), ("X+Y", "0")], ["X", "Y"],
                    {"X": True, "Y": True}),
    "cubic": ([("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")], ["A"], {"A": False}),
}


# the truth above is for generic rate constants, the question both routes answer; where unconditional ACR
# (every k, García Puente et al. Definition 2.17) differs, it is recorded here
UNCONDITIONAL = {"GP3.12": {"A": False, "B": False}}


def envz_copies(k):
    """k disjoint copies of the Figure 1 network: a size ladder with a known answer (Yp_i in every copy)"""
    R = [(_ren(a, i), _ren(b, i)) for i in range(k) for a, b in ENVZ]
    return R, [f"{s}_{i}" for i in range(k) for s in ENVZ_SPECIES]


def _ren(c, i):
    return "+".join(f"{t}_{i}" for t in c.split("+")) if c != "0" else c


def ecoli_core(name="e_coli_core"):
    """e_coli_core with reversible reactions split and blocked reactions removed (FVA), as (A, Y, N, species)"""
    import cobra
    cobra.Configuration().solver = "glpk"                 # Gurobi is node-locked on the cluster
    from kimonu.io import build_crn
    # the .mat copy: the SBML reader needs libSBML, which this environment lacks
    m = cobra.io.load_matlab_model(os.path.join(os.path.dirname(os.path.dirname(B_DIR)), "shared_models", f"{name}.mat"))
    for r in list(m.reactions):
        if r.lower_bound < 0:
            rv = cobra.Reaction(r.id + "_rev", lower_bound=0, upper_bound=-r.lower_bound if r.lower_bound > -1e9 else 1000)
            m.add_reactions([rv]); rv.add_metabolites({k: -v for k, v in r.metabolites.items()})
            r.lower_bound = 0
    ranges = sd.fva(m, solver="glpk")              # Gurobi is node-locked on the cluster
    m.remove_reactions([r for r in ranges.index if abs(ranges.at[r, "maximum"]) <= 1e-9], remove_orphans=True)
    crn = build_crn(m)                              # the zero complex has its own row of A
    return sp.csr_matrix(crn.A, dtype=np.int64), sp.csr_matrix(crn.Y), crn.N, crn.metabolite_ids


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def timed(f, limit):
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(int(limit)); t0 = time.time()
    try:
        return f(), time.time() - t0
    except Timeout:
        return None, time.time() - t0
    finally:
        signal.alarm(0)


# ---- the published criterion, symbolically
def symbolic_setup(A, Y, N):
    """J(a) = N_red diag(W a) B^T over Q[a]: N_red a row basis of N, B = substrate compositions, W a basis of ker N"""
    triples, (m, q) = entries(N)
    Nr = sympy.zeros(m, q)
    for i, j, v in triples: Nr[i, j] = sympy.Rational(v.numerator, v.denominator)
    rows = Nr.T.rref()[1]                                   # pivot rows of N = a row basis
    Nred = Nr.extract(list(rows), list(range(q)))
    W = sympy.Matrix.hstack(*Nr.nullspace())
    d = W.shape[1]
    a = sympy.symbols(f"a0:{d}")
    v = W * sympy.Matrix(a)
    Ac = sp.csc_matrix(A); Yc = sp.csc_matrix(Y)
    Bt = sympy.zeros(q, m)                                  # row j = composition of reaction j's substrate complex
    for j in range(q):
        for c, x in zip(Ac.indices[Ac.indptr[j]:Ac.indptr[j + 1]], Ac.data[Ac.indptr[j]:Ac.indptr[j + 1]]):
            if x < 0:
                for s_, y in zip(Yc.indices[Yc.indptr[c]:Yc.indptr[c + 1]], Yc.data[Yc.indptr[c]:Yc.indptr[c + 1]]):
                    Bt[j, s_] = sympy.Rational(str(y))
    J = Nred * sympy.diag(*v) * Bt
    return J, a, len(rows), m


def minors_acr(J, a, s, m, deadline):
    """Theorem 5.3(3) as stated: species i is local ACR iff every s x s minor of J without column i is the zero
    polynomial; returns (acr list, minors evaluated)"""
    dom = sympy.QQ[a] if a else sympy.QQ
    Jd = DomainMatrix.from_Matrix(J).convert_to(dom)
    acr, count = [], 0
    for i in range(m):
        cols = [c for c in range(m) if c != i]
        zero = True
        for comb in combinations(cols, s):
            if time.time() > deadline: raise Timeout()
            count += 1
            if not Jd.extract(list(range(s)), list(comb)).det().is_zero:
                zero = False; break
        if zero: acr.append(i)
    return acr, count


def rref_readout(J, a, m):
    """the same rank condition over Q(a) by Gaussian elimination: e_i in rowspace -> ACR, e_i - e_j -> ACRR"""
    dom = sympy.QQ.frac_field(*a) if a else sympy.QQ
    R, piv = DomainMatrix.from_Matrix(J).convert_to(dom).rref()
    R = R.to_list(); rank = len(piv); piv = list(piv)
    def in_rowspace(c):
        # RREF: c is in the row space iff c minus its pivot-coordinate combination of rows is zero
        for col in range(m):
            if col in piv: continue
            x = dom.convert(c.get(col, 0)) - sum((dom.convert(c.get(p, 0)) * R[k][col] for k, p in enumerate(piv)), dom.zero)
            if x != dom.zero: return False
        return True
    acr = [i for i in range(m) if in_rowspace({i: 1})]
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m) if i not in acr and j not in acr and in_rowspace({i: 1, j: -1})]
    return acr, pairs, rank


def groups_to_pairs(groups):
    return sorted((min(i, j), max(i, j)) for g in groups for i, j in combinations(sorted(g), 2))


def run_network(name, A, Y, N, species, truth=None, limit=3600):
    m = len(species)
    rec = dict(network=name, species=m, reactions=int(A.shape[1]), complexes=int(A.shape[0]))
    # tangent route (ours)
    res, t = timed(lambda: one_shot_fast(A, Y, N, False, want_groups=True), limit)
    _, acr_k, groups_k, dimT = res
    rec["kimonu"] = dict(seconds=t, acr=sorted(acr_k), acrr=groups_to_pairs(groups_k), dim_T=dimT)
    # iterated route (proven)
    it = kimonu.iterated_route(A, Y, N, detach_uptakes=False)
    rec["iterated"] = dict(acr=sorted(int(x) for x in it.acr), acrr=sorted((min(i, j), max(i, j)) for i, j in it.acrr))
    # symbolic criterion
    t0 = time.time()
    J, a, s, _ = symbolic_setup(A, Y, N)
    rec["setup_seconds"] = time.time() - t0; rec["s"] = s; rec["kernel_params"] = len(a)
    res, t = timed(lambda: rref_readout(J, a, m), limit)
    rec["rref"] = None if res is None else dict(seconds=t, acr=res[0], acrr=res[1], rank=res[2])
    if res is None: rec["rref_timeout"] = t
    rec["degenerate"] = res is not None and res[2] < s
    if os.environ.get("LACR_SKIP_MINORS"):
        res, t = None, float("nan")
    else:
        res, t = timed(lambda: minors_acr(J, a, s, m, time.time() + limit), limit)
    rec["minors"] = None if res is None else dict(seconds=t, acr=res[0], n_minors=res[1])
    if res is None and not os.environ.get("LACR_SKIP_MINORS"): rec["minors_timeout"] = t
    if os.environ.get("LACR_SKIP_MINORS"): rec["minors_skipped"] = "timed out at a smaller ladder step"
    # agreement
    k = rec["kimonu"]
    rec["agree_rref"] = rec["rref"] is not None and (rec["rref"]["acr"], sorted(map(tuple, rec["rref"]["acrr"]))) == (k["acr"], k["acrr"])
    rec["agree_minors"] = rec["minors"] is not None and (rec["minors"]["acr"] == k["acr"] or rec["degenerate"])
    if truth:
        rec["truth"] = {species.index(s_): v for s_, v in truth.items()}
        rec["chain_ok"] = all((i in rec["iterated"]["acr"]) <= v and v <= (i in k["acr"]) for i, v in rec["truth"].items())
    rec["species_names"] = species
    return rec


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "small":
        out = []
        for name, (R, S, truth) in NETWORKS.items():
            A, Y, N, cx = matrices(R, S)
            r = run_network(name, A, Y, N, S, truth, limit=600)
            if name in UNCONDITIONAL: r["truth_unconditional"] = UNCONDITIONAL[name]
            out.append(r)
            print(f"{name:12s} m={r['species']:2d} s={r['s']:2d} d={r['kernel_params']:2d} | kimonu {r['kimonu']['seconds']:6.2f}s "
                  f"rref {('%.2fs' % r['rref']['seconds']) if r['rref'] else 'timeout'} minors "
                  f"{('%.2fs' % r['minors']['seconds']) if r['minors'] else 'timeout'} | agree {r['agree_rref']} {r['agree_minors']} "
                  f"deg {r['degenerate']} chain {r.get('chain_ok')}", flush=True)
        json.dump(out, open(f"{OUT}/small.json", "w"), indent=1, default=str)
    else:
        name = sys.argv[2]; limit = float(sys.argv[3]) if len(sys.argv) > 3 else 3600
        if name.startswith("envz"):
            R, S = envz_copies(int(name[4:])); A, Y, N, cx = matrices(R, S)
        else:
            A, Y, N, S = ecoli_core(name)
        r = run_network(name, A, Y, N, S, limit=limit)
        json.dump(r, open(f"{OUT}/ladder_{name}.json", "w"), indent=1, default=str)
        print(json.dumps({k: v for k, v in r.items() if k not in ("species_names",)}, default=str)[:2000], flush=True)
