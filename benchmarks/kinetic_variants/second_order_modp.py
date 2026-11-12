"""Second-order test of the tangent-only relations, exact modulo a prime (attached uptakes).

second_order.py did this in float64 and could not tell anything apart: at an integer r* of scale 10^6 the
exact kernel spans too many orders of magnitude, lsqr stalls at a relative residual of 1e-3, and relations the
iterated route proves came out larger than the candidates. Here every step is exact over GF(p).

In log coordinates z = ln x - ln x*, the steady states are F(z) = N (r o exp(G z)) = W u(z) = 0, with
G = S^T Y^T. A curve z(t) = t h + t^2 q / 2 + ... of steady states has h in ker(J), J = N diag(r) G, and

    J q = -N (r o (G h)^2).

A relation v orthogonal to ker(J) is fixed to second order along h exactly when v^T q = 0 (v^T q does not
depend on the choice of q). The curve exists for every h when z* is regular, rank J = rank W, since F maps
into im W. All of this is polynomial in r on ker N, so it can be tested at a random r in ker N mod p: a
nonzero v^T q mod p is nonzero over Q at that r (and so at generic real r), and zero for random r and h means
zero with high probability (Schwartz-Zippel).

Candidates are the tangent-only relations of the exact routes at the real r* (tangent route, scale 10^6, and
the iterated route from src): ACR e_i, ACRR e_a - e_b, coupled complex pairs y_c - y_d. Checks reported:
each candidate is still orthogonal to ker(J) at the random r mod p (the first-order claim at a second,
independent r); proven relations give v^T q = 0; random v = J^T w give v^T q != 0 (the test has power).

    python kinetic_variants/second_order_modp.py <model>   -> json line, also $SO_OUT/<model>_modp.json
"""
import json, os, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import straindesign as sd
import kimonu
from kimonu.linalg import entries
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.tangent_fast import integer_positive_flux

OUT = os.environ.get("SO_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_second_order")
P = 2147483647                                   # 2^31 - 1: products of two residues fit in int64
N_H = 4
rng = np.random.default_rng(7)


def modp(x):
    x = Fraction(x)
    return x.numerator % P * pow(x.denominator % P, P - 2, P) % P


def rref(A):
    """reduced row echelon form mod P, in place; returns (A, pivot columns)"""
    A %= P
    nr, ncol = A.shape
    piv, r = [], 0
    for c in range(ncol):
        if r == nr: break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0: continue
        k = r + nz[0]
        if k != r: A[[r, k]] = A[[k, r]]
        A[r, c:] = A[r, c:] * pow(int(A[r, c]), P - 2, P) % P
        rows = np.nonzero(A[:, c])[0]
        rows = rows[rows != r]
        if rows.size:
            A[rows, c:] = (A[rows, c:] - (A[rows, c][:, None] * A[r, c:][None, :]) % P) % P
        piv.append(c); r += 1
    return A, piv


def mm(A, B):
    """A @ B mod P for int64 residues: B split into 16-bit halves so no partial sum overflows"""
    lo, hi = B & 0xFFFF, B >> 16
    return ((A @ lo) % P + ((A @ hi) % P) * 65536 % P) % P


def kernel_sample(R, piv, ncol, k):
    """k random kernel vectors from an rref (ncol x k)"""
    free = np.array(sorted(set(range(ncol)) - set(piv)), dtype=np.int64)
    X = np.zeros((ncol, k), dtype=np.int64)
    X[free] = rng.integers(0, P, size=(free.size, k))
    X[np.array(piv, dtype=np.int64)] = (-mm(R[:len(piv)][:, free], X[free])) % P
    return X


def kernel_basis(R, piv, ncol):
    """kernel basis mod P from an rref, ncol x (#free)"""
    free = sorted(set(range(ncol)) - set(piv))
    K = np.zeros((ncol, len(free)), dtype=np.int64)
    for t, f in enumerate(free):
        K[f, t] = 1
    for i, c in enumerate(piv):
        K[c] = (-R[i, free]) % P
    return K


def pairs(groups):
    return {(a, b) for g in groups for a in g for b in g if a < b}


M = sys.argv[1]
t_start = time.time()
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
free_mets = set(range(mets.index("E_1") if "E_1" in mets else nm))

# ---- exact routes at the real r* (as in hierarchy_data.py)
it = kimonu.iterated_route(kimonu.crn_from_matrices(A, Y), detach_uptakes=False)
flux = integer_positive_flux(S, seed=0)
Ac = sp.csc_matrix(A); Yc = sp.csc_matrix(Y)
sub = {}
for j in range(nr):
    for i, a in zip(Ac.indices[Ac.indptr[j]:Ac.indptr[j + 1]], Ac.data[Ac.indptr[j]:Ac.indptr[j + 1]]):
        if a < 0: sub[j] = int(i)
ycol = {c: [(int(i), Fraction(str(v))) for i, v in zip(Yc.indices[Yc.indptr[c]:Yc.indptr[c + 1]],
                                                      Yc.data[Yc.indptr[c]:Yc.indptr[c + 1]])] for c in range(nc)}
triples, _ = entries(S)
byr = {}
for i, j, v in triples: byr.setdefault(j, []).append((i, v))
J = {}
for j, s in sub.items():
    for i, v in byr.get(j, ()):
        for k, w in ycol[s]:
            J[(i, k)] = J.get((i, k), Fraction(0)) + flux[j] * v * w
K = sd.sparse_nullspace(sd.RationalMatrix.from_fractions(((i, k, x) for (i, k), x in J.items() if x != 0), (nm, nm)))
rows = {}
for i, j, v in entries(K)[0]: rows.setdefault(i, {})[j] = v
key = lambda dd: tuple(sorted((c, v) for c, v in dd.items() if v != 0))
groups = {}
for c in range(nc):
    acc = {}
    for i, w in ycol[c]:
        for j, v in rows.get(i, {}).items(): acc[j] = acc.get(j, 0) + w * v
    groups.setdefault(key(acc), set()).add(c)
acr_t = {i for i in range(nm) if not rows.get(i)}
metg = {}
for i in range(nm):
    if rows.get(i): metg.setdefault(key(rows[i]), []).append(i)
acrr_t = pairs([sorted(g) for g in metg.values() if len(g) > 1])
cp_t = pairs([sorted(g) for g in groups.values()])
cp_it = pairs([sorted(m) for m in it.modules])
acr_it = set(it.acr); acrr_it = {tuple(sorted(p)) for p in it.acrr}
t_routes = time.time() - t_start

# ---- GF(p): random r in ker N, J = N diag(r) G, W = N diag(r) S^T
t0 = time.time()
Np = np.zeros((nm, nr), dtype=np.int64)
for i, j, v in triples: Np[i, j] = modp(v)
RN, pivN = rref(Np.copy())
rank_N = len(pivN)
r = kernel_sample(RN, pivN, nr, 1)[:, 0]
assert not mm(Np, r[:, None]).any()
Yp = {c: [(i, modp(w)) for i, w in lst] for c, lst in ycol.items()}
Jp = np.zeros((nm, nm), dtype=np.int64); Wp = np.zeros((nm, nc), dtype=np.int64)
for j, s in sub.items():
    for i, v in byr.get(j, ()):
        a = modp(v) * int(r[j]) % P
        Wp[i, s] = (Wp[i, s] + a) % P
        for k, w in Yp[s]:
            Jp[i, k] = (Jp[i, k] + a * w) % P
rank_W = len(rref(Wp.copy())[1])
RJ, pivJ = rref(Jp.copy())
rank_J = len(pivJ)
KJ = kernel_basis(RJ, pivJ, nm)                      # nm x d
d = KJ.shape[1]
# random h in ker J; G h per reaction (substrate composition . h); rhs = -N (r o (G h)^2)
H = mm(KJ, rng.integers(0, P, size=(d, N_H)))
assert not mm(Jp, H).any()
GH = np.zeros((nr, N_H), dtype=np.int64)
for j, s in sub.items():
    for k, w in Yp[s]:
        GH[j] = (GH[j] + w * H[k] % P) % P
coef = r[:, None] * (GH * GH % P) % P                     # nr x N_H
RHS = (-mm(Np, coef)) % P
aug = np.concatenate([Jp, RHS], axis=1)
Raug, pivA = rref(aug)
consistent = all(c < nm for c in pivA)
Q = np.zeros((nm, N_H), dtype=np.int64)
if consistent:
    for i, c in enumerate(pivA):
        Q[c] = Raug[i, nm:]
t_modp = time.time() - t0


def yvec(c):
    return {i: w for i, w in Yp[c]}


def evec(i):
    return {i: 1}


def diff(a, b):
    v = dict(a)
    for i, w in b.items(): v[i] = (v.get(i, 0) - w) % P
    return {i: w for i, w in v.items() if w}


def test(v):
    """(orthogonal to ker J mod p, v^T q nonzero for some h); v sparse {metabolite: residue}"""
    t, z = np.zeros(d, dtype=np.int64), np.zeros(N_H, dtype=np.int64)
    for i, w in v.items():
        t = (t + w * KJ[i] % P) % P
        z = (z + w * Q[i] % P) % P
    return (not t.any()), bool(z.any())


cand = {"acr": [(i, evec(i)) for i in sorted(acr_t - acr_it)],
        "acrr": [((a, b), diff(evec(a), evec(b))) for a, b in sorted(acrr_t - acrr_it)],
        "complex_pairs": [((c, e), diff(yvec(c), yvec(e))) for c, e in sorted(cp_t - cp_it)]}
proven = ([diff(evec(a), evec(b)) for a, b in sorted(acrr_it)] + [evec(i) for i in sorted(acr_it)]
          + [diff(yvec(c), yvec(e)) for c, e in sorted(cp_it)][:2000])
randrow = [{i: int(w) for i, w in enumerate(row) if w} for row in mm(rng.integers(0, P, size=(20, nm)), Jp)]

out = dict(model=M, p=P, n_metabolites=nm, dim_T_real=K.shape[1], dim_T_modp=d, rank_N=rank_N, rank_W=rank_W,
           rank_J=rank_J, regular=rank_J == rank_W, second_order_consistent=consistent, n_h=N_H,
           seconds_routes=t_routes, seconds_modp=t_modp)
res = [test(v) for v in proven]
out["control_proven"] = dict(n=len(res), not_orthogonal=sum(not o for o, _ in res), nonzero=sum(z for _, z in res))
res = [test(v) for v in randrow]
out["control_random_rowspace"] = dict(n=len(res), not_orthogonal=sum(not o for o, _ in res), nonzero=sum(z for _, z in res))
out["candidates"] = {}
for k, lst in cand.items():
    res = [test(v) for _, v in lst]
    isfree = [(x in free_mets) if k == "acr" else (k == "acrr" and x[0] in free_mets and x[1] in free_mets) for x, _ in lst]
    out["candidates"][k] = dict(n=len(lst), n_free=sum(isfree), not_orthogonal=sum(not o for o, _ in res),
                                nonzero=sum(z for _, z in res),
                                nonzero_items=[str(x) for (x, _), (_, z) in zip(lst, res) if z][:50])
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}_modp.json", "w"), indent=1)
print(json.dumps(out), flush=True)
