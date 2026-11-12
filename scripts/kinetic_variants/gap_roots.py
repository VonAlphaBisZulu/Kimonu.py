"""How much of the gap follows from its ACR relations: fix them, rerun the iterated route, recount.

If the metabolites in ROOTS have ACR, their concentrations are constants at every positive steady state of
given rate constants, and dropping them from every complex (their constant factor goes into the rate
constants, their balance is dropped) gives a network whose steady states contain the original ones. Every
relation the iterated route proves on that network then holds on the original network too, provided the
roots hold. The script reports which of the model's open relations (two_point_gap.py prepare, attached
uptakes) the iterated route proves after fixing the roots.

A root "a/b" is an ACRR relation instead: x_b = c x_a at every positive steady state, so b is replaced by a in
every complex (the constant goes into the rate constants) and the two balances are summed, which the
original steady states satisfy as well.

    python kinetic_variants/gap_roots.py <model> [root ...]   # root: metabolite (ACR) or a/b (ACRR); default: the open ACR relations

DETACH=1 runs the iterated route with detached uptakes (the open relations are the same in both modes).
"""
import json, os, sys
from fractions import Fraction
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import records                                         # paths of the records (KIMONU_RECORDS, KIMONU_OUT)
import kimonu
from kimonu.linalg import exact
from run_kimonu_crn import crn_of

M = sys.argv[1]
DETACH = os.environ.get("DETACH") == "1"
crn = json.load(open(records.network(M)))
mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
rel = json.load(open(records.read(records.TWO_POINT, f"{M}_relations.json")))["relations"]   # two_point_gap.py prepare
roots = sys.argv[2:] or [next(iter(r["c"])) for r in rel if r["kind"] == "ACR"]

# structured reactions (substrate, product) from the composition of each reaction's two complexes
ycol, sub, prod = {}, {}, {}
for i, c, v in crn["Y"]:
    ycol.setdefault(int(c) - 1, {})[mets[int(i) - 1]] = Fraction(str(v))
for c, j, v in crn["A"]:
    (sub if float(v) < 0 else prod)[int(j) - 1] = int(c) - 1
rx = [(dict(ycol.get(sub[j], {})), dict(ycol.get(prod[j], {}))) for j in range(crn["n_reactions"])]

# build() as in reduce_stage3.py: prune blocked reactions by one LP, merge complexes of equal composition
def build(rx):
    ms = sorted({n for s, p in rx for n in list(s) + list(p)}); mi = {n: i for i, n in enumerate(ms)}
    q = len(rx)
    Nf = sp.lil_matrix((len(ms), q))
    for j, (s, p) in enumerate(rx):
        for n, v in s.items(): Nf[mi[n], j] -= float(v)
        for n, v in p.items(): Nf[mi[n], j] += float(v)
    Nf = Nf.tocsr()
    lp = linprog(np.concatenate([np.zeros(q), -np.ones(q)]), A_ub=sp.hstack([-sp.identity(q), sp.identity(q)]).tocsr(),
                 b_ub=np.zeros(q), A_eq=sp.hstack([Nf, sp.csr_matrix((len(ms), q))]).tocsr(), b_eq=np.zeros(len(ms)),
                 bounds=[(0, None)] * q + [(0, 1)] * q, method="highs")
    rx = [rx[j] for j in range(q) if lp.status == 0 and lp.x[q + j] > 0.5]
    ms = sorted({n for s, p in rx for n in list(s) + list(p)}); mi = {n: i for i, n in enumerate(ms)}
    cx, cid, Ae, Ne = [], {}, [], {}
    def c(side):
        key = tuple(sorted(side.items()))
        if key not in cid: cid[key] = len(cx); cx.append(key)
        return cid[key]
    for j, (s, p) in enumerate(rx):
        Ae += [(c(s), j, -1), (c(p), j, 1)]
        for n, v in s.items(): Ne[(mi[n], j)] = Ne.get((mi[n], j), 0) - v
        for n, v in p.items(): Ne[(mi[n], j)] = Ne.get((mi[n], j), 0) + v
    A = sp.csr_matrix(([v for *_, v in Ae], ([i for i, *_ in Ae], [j for _, j, _ in Ae])), shape=(len(cx), len(rx)), dtype=np.int64)
    Ye = [(mi[n], k, float(v)) for k, key in enumerate(cx) for n, v in key]
    Y = sp.csr_matrix(([v for *_, v in Ye], ([i for i, *_ in Ye], [k for _, k, _ in Ye])), shape=(len(ms), len(cx)))
    N = exact(((i, j, v) for (i, j), v in Ne.items() if v != 0), (len(ms), len(rx)))
    return A, Y, N, ms, rx

def fixed(rx, roots):
    drop = {r for r in roots if "/" not in r}
    alias = {}                                    # b -> a for every ratio root a/b, chains resolved
    for r in roots:
        if "/" in r:
            a, b = r.split("/")
            while a in alias: a = alias[a]
            alias[b] = a
    def name(n):
        while n in alias: n = alias[n]
        return n
    def side(d):
        out = {}
        for n, v in d.items():
            if n not in drop:
                out[name(n)] = out.get(name(n), 0) + v
        return out
    out = []
    for s, p in rx:
        s2, p2 = side(s), side(p)
        if s2 != p2:
            out.append((s2, p2))
    return out


def proven(A, Y, N, ms):
    it = kimonu.iterated_route(crn_of(A, Y, N), detach_uptakes=DETACH)
    acr = {ms[i] for i in it.acr}
    parent = {}
    def find(a):
        while parent.setdefault(a, a) != a:
            a = parent[a]
        return a
    for a, b in it.acrr:
        parent[find(ms[a])] = find(ms[b])
    return acr, find


A, Y, N, ms, _ = build(rx)
acr0, find0 = proven(A, Y, N, ms)
A, Y, N, ms, kept = build(fixed(rx, set(roots)))
acr1, find1 = proven(A, Y, N, ms)
alias = {}
for r in roots:
    if "/" in r:
        a, b = r.split("/")
        while a in alias: a = alias[a]
        alias[b] = a
def name(n):
    while n in alias: n = alias[n]
    return n
fixed_acr = {r for r in roots if "/" not in r} | acr1
closed, still = [], []
for r in rel:
    names = list(r["c"])
    if r["kind"] == "ACR":
        ok = name(names[0]) in fixed_acr
    else:
        a, b = (name(n) for n in names)
        ok = a == b or (a in fixed_acr and b in fixed_acr) or (a in ms and b in ms and find1(a) == find1(b))
    (closed if ok else still).append(r)
print(json.dumps(dict(model=M, detached=DETACH, roots=roots, reactions_kept=len(kept), n_relations=len(rel), closed=len(closed),
                      open=[(r["kind"], r["c"]) for r in still]), indent=1))
