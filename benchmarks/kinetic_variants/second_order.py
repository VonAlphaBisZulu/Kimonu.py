"""Second-order test of the tangent-only relations on the gap models (attached uptakes).

In log coordinates z = ln x - ln x*, the steady states are F(z) = N (r* o exp(G z)) = 0, with G = S^T Y^T
(row j = composition of reaction j's substrate complex; attached uptakes give a zero row). A curve of steady
states through z* = 0 is z(t) = t h + t^2 q / 2 + ..., with h in T = ker(J), J = N diag(r*) G, and

    J q = -N (r* o (G h)^2).

A relation v with v ⊥ T is fixed to second order along h exactly when v^T q = 0; v^T q does not depend on
which q solves the system, since v lies in the row space of J. If v^T q != 0 for some h, and z* is a regular
point (rank J = rank N, so that the curve exists), v is not fixed.

Candidates: tangent-only ACR (e_i), ACRR pairs (e_a - e_b) and coupled complex pairs (y_c - y_d), i.e. those
in the tangent route and not in the iterated route. Controls: the same test on relations the iterated route
proves (must give 0), and on random v = J^T w (must give nonzero for the test to have power). Float64 on
the exact kernel and flux, several random h.

    python kinetic_variants/second_order.py <model>   -> json line, also $SO_OUT/<model>.json
"""
import json, os, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import lsqr
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import straindesign as sd
import kimonu
from kimonu.linalg import entries
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.tangent_fast import integer_positive_flux

OUT = os.environ.get("SO_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_second_order")
N_H = 8


def pairs(groups):
    return {(a, b) for g in groups for a in g for b in g if a < b}


M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
free = set(range(mets.index("E_1") if "E_1" in mets else nm))

# iterated route
t0 = time.time()
it = kimonu.kinetic_modules(A, Y, S, detach_uptakes=False)
ar = kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr)
t_it = time.time() - t0

# tangent route, exact, as in hierarchy_data.py
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
J = {k: x for k, x in J.items() if x != 0}
K = sd.sparse_nullspace(sd.RationalMatrix.from_fractions(((i, k, x) for (i, k), x in J.items()), (nm, nm)))
d = K.shape[1]
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
acr_it = set(ar.acr); acrr_it = {tuple(sorted(p)) for p in ar.acrr}

# float objects
Jf = sp.csr_matrix(([float(x) for x in J.values()], ([i for i, _ in J], [k for _, k in J])), shape=(nm, nm))
Kt, _ = entries(K)
Kf = sp.csr_matrix(([float(v) for *_, v in Kt], ([i for i, *_ in Kt], [j for _, j, _ in Kt])), shape=(nm, d)).toarray()
Nf = sp.csr_matrix(([float(v) for *_, v in triples], ([i for i, *_ in triples], [j for _, j, _ in triples])), shape=(nm, nr))
rf = np.array([float(x) for x in flux])
G = sp.lil_matrix((nr, nm))
for j, s in sub.items():
    for k, w in ycol[s]: G[j, k] = float(w)
G = G.tocsr()
Yf = Yc.astype(float)
assert abs(Jf - Nf @ sp.diags(rf) @ G).max() < 1e-6 * abs(Jf).max()
# regularity: rank J = rank N (float, dense SVD is fine at these sizes)
rankN = np.linalg.matrix_rank(Nf.toarray()); rankJ = nm - d

rng = np.random.default_rng(1)
H, Q, res = [], [], []
for _ in range(N_H):
    h = Kf @ rng.standard_normal(d)
    h /= np.linalg.norm(h)
    rhs = -(Nf @ (rf * (G @ h) ** 2))
    sol = lsqr(Jf, rhs, atol=1e-14, btol=1e-14, iter_lim=50 * nm)
    q = sol[0]
    res.append(float(np.linalg.norm(Jf @ q - rhs) / max(np.linalg.norm(rhs), 1e-300)))
    H.append(h); Q.append(q)
Q = np.array(Q)                                # N_H x nm
qscale = float(np.abs(Q).max())


def test(vecs):
    """max over h of |v^T q| / (||v|| * max|q|) per candidate"""
    out = []
    for v in vecs:
        v = np.asarray(v, float)
        out.append(float(np.abs(Q @ v).max() / (np.linalg.norm(v) * qscale)))
    return out


def e(i):
    v = np.zeros(nm); v[i] = 1; return v


def ydiff(c, dd):
    return np.asarray((Yf[:, c] - Yf[:, dd]).todense()).ravel()


cand = {"acr": [(i, e(i)) for i in sorted(acr_t - acr_it)],
        "acrr": [((a, b), e(a) - e(b)) for a, b in sorted(acrr_t - acrr_it)],
        "complex_pairs": [((c, dd), ydiff(c, dd)) for c, dd in sorted(cp_t - cp_it)]}
ctrl_rng = np.random.default_rng(2)
proven = [e(a) - e(b) for a, b in list(sorted(acrr_it))[:200]] + [e(i) for i in sorted(acr_it)]
randrow = [Jf.T @ ctrl_rng.standard_normal(nm) for _ in range(20)]

out = dict(model=M, n_metabolites=nm, dim_T=d, rank_J=int(rankJ), rank_N=int(rankN), regular=bool(rankJ == rankN),
           seconds_iterated=t_it, n_h=N_H, lsqr_rel_residual_max=max(res),
           first_order_check=float(max(np.abs(Jf @ h).max() for h in H)),
           control_proven_max=max(test(proven), default=None),
           control_random_rowspace_min=min(test(randrow)),
           candidates={})
for k, lst in cand.items():
    vals = test([v for _, v in lst])
    out["candidates"][k] = dict(n=len(lst), max=max(vals, default=None),
                                n_free=sum(1 for (x, _) in lst if (x in free if k == "acr" else
                                           (k == "acrr" and x[0] in free and x[1] in free))),
                                values=sorted(((str(x), val) for (x, _), val in zip(lst, vals)), key=lambda t: -t[1])[:40])
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "candidates"}), flush=True)
for k, v in out["candidates"].items():
    print(k, v["n"], v["n_free"], v["max"], flush=True)
