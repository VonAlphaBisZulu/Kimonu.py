"""Certificate of genericity: both routes again at a uniformly random point, modulo the prime p = 2^61 - 1.

Kimonu's relations are read off one interior flux r* (kimonu.interior_flux). They hold for
generic rate constants unless r* lies on the zero set of one of finitely many polynomials in r: the minors of
W(r) or J(r) = W(r) Y^T that decide its rank and whether a candidate vector (e_i, e_i - e_j, y_c - y_d) lies in
its row space. Their entries are linear in r, so these polynomials have degree at most rank + 1. The recorded
run used the LP point snapped exactly onto ker N (now draw="lp"), which is not a random draw, so no probability
bound applies to it. The opt-in uniform draw has a bound of its own.

Here the readout is computed once more at r' = K b, with b uniform over GF(p)^d, in arithmetic modulo p.
A polynomial that is nonzero over Q (and not identically zero mod p) vanishes at r' with probability at most
deg / p (Schwartz-Zippel). If the readout at r' equals the readout at r*, the result at r* can be non-generic
only if r' is non-generic in the same way, so

    P(result at r* non-generic and certificate passes) <= (number of tests) * (rank + 1) / p,

with one test per candidate relation (m ACR, m(m-1)/2 ACRR, n(n-1)/2 complex pairs) and one for the rank.
The iterated route computes ker W in modular arithmetic and its log-space closure exactly (that part does not
depend on r), round by round as kimonu.iterated_route. Detached uptakes.

    python cert_modp.py <model> [seed]   -> $CM_OUT/<model>.json
"""
import json, os, sys, time
from fractions import Fraction
import numpy as np
import flint
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu import free_directions, readout, steady_state_map
from kimonu.network import zero_complexes
from kimonu.linalg import nullspace
from run_kimonu_crn import triplets_to_sparse

OUT = os.environ.get("CM_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_cert_modp")
P = 2**61 - 1
M = sys.argv[1]; SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 0
log = lambda *a: print(f"[{time.strftime('%H:%M:%S')}] {M}:", *a, flush=True)


def mod(x):
    x = Fraction(x)
    assert x.denominator % P, f"p divides a denominator ({x})"
    return x.numerator % P * pow(x.denominator, -1, P) % P


crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
net = kimonu.crn_from_matrices(A, Y)                 # N = Y·A, the S of the crn.json
skip = zero_complexes(net)
Ycol = [[(k, mod(v)) for k, v in y] for y in net.Yt]

# ---- the exact results of the default pipeline
t0 = time.time()
tr = kimonu.tangent_route(net, detach_uptakes=True)
km = kimonu.iterated_route(net, detach_uptakes=True)
mods = lambda ms: sorted(sorted(int(c) for c in s) for s in ms if len(s) > 1)
exact_res = dict(tangent=dict(modules=mods(tr.modules), acr=sorted(map(int, tr.acr)),
                              acrr=sorted(tuple(sorted(map(int, a))) for a in tr.acrr), dim=tr.dimension),
                 iterated=dict(modules=mods(km.modules), acr=sorted(map(int, km.acr)),
                               acrr=sorted(tuple(sorted(map(int, a))) for a in km.acrr)))
log("exact", round(time.time() - t0, 1), "s")

# ---- a uniform point of ker N mod p
t0 = time.time()
K = nullspace(net.N)
d = max((j for r in K for j, _ in r), default=-1) + 1
rng = np.random.default_rng(SEED)
b = [int(rng.integers(0, P)) for _ in range(d)]
Kmod = [[(j, mod(v)) for j, v in r] for r in K]
rp = [sum(v * b[j] for j, v in r) % P for r in Kmod]
W, ncols = steady_state_map(net, rp, True)          # entries r'_j * N_ij, exact integers times Fractions
m = len(W)
Wm = {(i, c): mod(v) for i, row in enumerate(W) for c, v in row}
Wm = {key: v for key, v in Wm.items() if v}


def kernel_rows(entries_, shape):
    """Rows of a basis of ker(M) mod p, one tuple per column of M, and the nullity."""
    Mx = flint.nmod_mat(shape[0], shape[1], P)
    for (i, j), v in entries_.items():
        Mx[i, j] = v
    X, nul = Mx.nullspace()
    return [tuple(int(X[i, j]) for j in range(nul)) for i in range(shape[1])], nul


def normalized(row):
    """The row scaled so its first nonzero is 1 (proportional rows become equal); None for a zero row."""
    for v in row:
        if v:
            inv = pow(v, -1, P)
            return tuple(x * inv % P for x in row)
    return None


def groups(keys, items):
    g = {}
    for it, k in zip(items, keys):
        g.setdefault(k, []).append(it)
    return list(g.values())


# ---- tangent route at r' mod p: T' = ker(W Y^T)
J = {}
for (i, c), w in Wm.items():
    if c >= nc:
        J[(i, nm + c - nc)] = (J.get((i, nm + c - nc), 0) + w) % P
        continue
    for k, y in Ycol[c]:
        J[(i, k)] = (J.get((i, k), 0) + w * y) % P
J = {key: v for key, v in J.items() if v}
rowsT, nulT = kernel_rows(J, (nm, nm + ncols - nc))
Tm = rowsT[:nm]
acr_p = sorted(i for i, r in enumerate(Tm) if not any(r))
acrr_p = sorted(tuple(sorted((a, b2))) for g in groups(Tm, range(nm)) if any(Tm[g[0]]) and len(g) > 1
                for i, a in enumerate(g) for b2 in g[i + 1:])
YtT = {}
for c in range(nc):
    if c in skip: continue
    acc = [0] * nulT
    for k, y in Ycol[c]:
        for j, v in enumerate(Tm[k]):
            if v: acc[j] = (acc[j] + y * v) % P
    YtT[c] = tuple(acc)
mod_p = mods(groups(YtT.values(), YtT.keys()))
rankJ = nm + ncols - nc - nulT
tang = dict(modules=mod_p, acr=acr_p, acrr=acrr_p, dim=nulT)
tang_ok = {k: tang[k] == exact_res["tangent"][k] for k in tang}
log("tangent mod p", round(time.time() - t0, 1), "s", tang_ok)

# ---- iterated route at r' mod p: ker W in modular arithmetic, the log closure exactly
t0 = time.time()
equal, rnd, ranks = [], 0, []
drop = skip
while True:
    rnd += 1
    Wx = dict(Wm)
    for k, (c, e) in enumerate(equal):
        Wx[(m + k, c)], Wx[(m + k, e)] = 1, P - 1
    rowsW, nulW = kernel_rows(Wx, (m + len(equal), ncols))
    ranks.append(ncols - nulW)
    coupled = [set(g) for g in groups([normalized(r) for r in rowsW[:nc]], range(nc))]
    coupled = [{c for c in s if c not in drop} for s in coupled]
    coupled = [s for s in coupled if s]
    # the readout lists the detached zero complex last as a module of its own, which coupled leaves out
    sets = [s for s in readout(net, free_directions(net, coupled)).modules if not s & skip]
    if len(sets) == len(coupled):
        break
    equal = [(min(s), c) for s in sets for c in s if c != min(s)]
Rp = readout(net, free_directions(net, sets))
it = dict(modules=mods(sets), acr=sorted(map(int, Rp.acr)),
          acrr=sorted(tuple(sorted((a, b2))) for g in Rp.acrr_groups for i, a in enumerate(g) for b2 in g[i + 1:]))
it_ok = {k: it[k] == exact_res["iterated"][k] for k in it}
log("iterated mod p", round(time.time() - t0, 1), "s, rounds", rnd, it_ok)

# ---- the bound
ncx = nc - len(skip)
tests = nm + nm * (nm - 1) // 2 + ncx * (ncx - 1) // 2 + 1
bound_t = tests * (rankJ + 1) / P
bound_i = rnd * tests * (max(ranks) + 1) / P
out = dict(model=M, seed=SEED, p=P, m=nm, n=nc, tangent_match=all(tang_ok.values()), tangent_detail=tang_ok,
           iterated_match=all(it_ok.values()), iterated_detail=it_ok, rank_J=rankJ, ranks_W=ranks, rounds=rnd,
           tests=tests, bound_tangent=bound_t, bound_iterated=bound_i)
log("certificate", {k: out[k] for k in ("tangent_match", "iterated_match", "bound_tangent", "bound_iterated")})
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"))
