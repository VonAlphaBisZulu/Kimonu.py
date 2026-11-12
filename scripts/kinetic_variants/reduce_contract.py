"""Shrink a network that shows a tangent-only relation by contracting enzyme mechanisms.

Deleting reactions from a genome model soon blocks biomass and with it the flux couplings a relation rests
on. Contracting keeps them: an enzyme mechanism (elementary steps that move one enzyme between its free
form E_k and bound forms Ek_...) is replaced by its net reactions, one per cycle of the enzyme through its
forms that returns to E_k with a nonzero net conversion. The net stoichiometry is unchanged, so the flux
cone of the free metabolites is; only the kinetics inside the contracted mechanism change. A contraction
is kept if the tangent space still fixes the target and the iterated route still does not (both exact).
Mechanisms are tried far-to-near from the target in contiguous blocks, as in reduce_gap.py. A second
stage then deletes units from the contracted network (each remaining elementary mechanism as a whole,
every other reaction on its own), same order and predicate, re-pruning blocked reactions each time.

    python kinetic_variants/reduce_contract.py <model> <target> <start.json> <out.json>
"""
import json, os, re, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import records                                         # paths of the records (KIMONU_RECORDS, KIMONU_OUT)
import kimonu
from kimonu.linalg import exact
from run_kimonu_crn import crn_of, triplets_to_sparse
from kinetic_variants.tangent_one_shot import one_shot_exact

M, target, start, out = sys.argv[1:5]
crn = json.load(open(records.network(M))); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A0 = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64).tocsc(); Y0 = triplets_to_sparse(crn["Y"], nm, nc).tocsc(); mets = crn["metabolites"]

def comp(c):
    return {mets[i]: Fraction(str(v)) for i, v in zip(Y0.indices[Y0.indptr[c]:Y0.indptr[c + 1]], Y0.data[Y0.indptr[c]:Y0.indptr[c + 1]])}
def reaction(j):
    col = A0[:, j].toarray().ravel()
    return comp(int(np.flatnonzero(col < 0)[0])), comp(int(np.flatnonzero(col > 0)[0]))

ENZ = re.compile(r"E_?(\d+)(?:_|$)")
def enzyme_of(name):
    m = ENZ.match(name); return m.group(1) if m else None

base = [reaction(j) for j in json.load(open(start))["reactions"]]
mech = {}
for k, (s, p) in enumerate(base):
    ids = {enzyme_of(n) for n in list(s) + list(p)} - {None}
    if len(ids) == 1: mech.setdefault(ids.pop(), []).append(k)

def net_reactions(eid, idx):
    """cycles of enzyme eid through its forms, starting and ending at the free form, with nonzero net"""
    free = f"E_{eid}"; edges = []
    for k in idx:
        s, p = base[k]
        fs = [n for n in s if enzyme_of(n) == eid]; fp = [n for n in p if enzyme_of(n) == eid]
        if len(fs) != 1 or len(fp) != 1 or s[fs[0]] != 1 or p[fp[0]] != 1: return None
        rest = {n: -v for n, v in s.items() if n != fs[0]}
        for n, v in p.items():
            if n != fp[0]: rest[n] = rest.get(n, 0) + v
        edges.append((fs[0], fp[0], rest))
    succ = {}
    for a, b, d in edges: succ.setdefault(a, []).append((b, d))
    if free not in succ: return None
    nets, stack = set(), [(free, {}, {free})]
    while stack:
        node, acc, seen = stack.pop()
        for b, d in succ.get(node, []):
            tot = dict(acc)
            for n, v in d.items(): tot[n] = tot.get(n, 0) + v
            if b == free:
                tot = {n: v for n, v in tot.items() if v != 0}
                if tot: nets.add(tuple(sorted(tot.items())))
            elif b not in seen and len(seen) < 30:
                stack.append((b, tot, seen | {b}))
    if not nets: return None
    return [({n: -v for n, v in t if v < 0}, {n: v for n, v in t if v > 0}) for t in nets]

contractible = {e: nr_ for e, idx in mech.items() if (nr_ := net_reactions(e, idx)) is not None}

def build(contracted):
    drop = {k for e in contracted for k in mech[e]}
    rx = [base[k] for k in range(len(base)) if k not in drop] + [r for e in contracted for r in contractible[e]]
    return build_list(rx)

def build_list(rx):
    ms = sorted({n for s, p in rx for n in list(s) + list(p)}); mi = {n: i for i, n in enumerate(ms)}
    q = len(rx); Nf = sp.lil_matrix((len(ms), q))
    for j, (s, p) in enumerate(rx):
        for n, v in s.items(): Nf[mi[n], j] -= float(v)
        for n, v in p.items(): Nf[mi[n], j] += float(v)
    Nf = Nf.tocsr()
    lp = linprog(np.concatenate([np.zeros(q), -np.ones(q)]), A_ub=sp.hstack([-sp.identity(q), sp.identity(q)]).tocsr(),
                 b_ub=np.zeros(q), A_eq=sp.hstack([Nf, sp.csr_matrix((len(ms), q))]).tocsr(), b_eq=np.zeros(len(ms)),
                 bounds=[(0, None)] * q + [(0, 1)] * q, method="highs")
    rx = [rx[j] for j in range(q) if lp.status == 0 and lp.x[q + j] > 0.5]
    ms = sorted({n for s, p in rx for n in list(s) + list(p)}); mi = {n: i for i, n in enumerate(ms)}
    cx, cid = [], {}
    def c(side):
        key = tuple(sorted(side.items()))
        if key not in cid: cid[key] = len(cx); cx.append(key)
        return cid[key]
    Ae, Ne = [], {}
    for j, (s, p) in enumerate(rx):
        Ae += [(c(s), j, -1), (c(p), j, 1)]
        for n, v in s.items(): Ne[(mi[n], j)] = Ne.get((mi[n], j), 0) - v
        for n, v in p.items(): Ne[(mi[n], j)] = Ne.get((mi[n], j), 0) + v
    A = sp.csr_matrix(([v for *_, v in Ae], ([i for i, *_ in Ae], [j for _, j, _ in Ae])), shape=(len(cx), len(rx)), dtype=np.int64)
    Ye = [(mi[n], k, float(v)) for k, key in enumerate(cx) for n, v in key]
    Y = sp.csr_matrix(([v for *_, v in Ye], ([i for i, *_ in Ye], [k for _, k, _ in Ye])), shape=(len(ms), len(cx)))
    N = exact(((i, j, v) for (i, j), v in Ne.items() if v != 0), (len(ms), len(rx)))
    return A, Y, N, ms, rx

cache = {}
def gap(contracted):
    key = frozenset(contracted)
    if key in cache: return cache[key]
    ok = False
    try:
        A, Y, N, ms, _ = build(contracted)
        if target in ms:
            ti = ms.index(target)
            it = kimonu.iterated_route(crn_of(A, Y, N), detach_uptakes=False)
            if ti not in set(it.acr):
                ok = ti in set(one_shot_exact(A, Y, N, False)[1])
    except Exception as e:
        print("   predicate error:", type(e).__name__, str(e)[:100], flush=True)
    cache[key] = ok
    return ok

def save(contracted, note):
    A, Y, N, ms, rx = build(contracted)
    fmt = lambda d: "+".join(f"{'' if v == 1 else v}{n}" for n, v in sorted(d.items())) or "0"
    json.dump(dict(model=M, target=target, note=note, contracted=sorted(contracted), n_reactions=len(rx),
                   n_metabolites=len(ms), reactions=[[fmt(s), fmt(p)] for s, p in rx],
                   structured=[[[[n, str(v)] for n, v in sorted(s.items())], [[n, str(v)] for n, v in sorted(p.items())]] for s, p in rx]),
              open(out, "w"), indent=1)
    return len(rx), len(ms)

# far-to-near order of the contractible mechanisms, hubs excluded
rx0 = base; deg = {}
for s, p in rx0:
    for n in set(s) | set(p): deg[n] = deg.get(n, 0) + 1
hubs = set(sorted(deg, key=lambda n: -deg[n])[:20]) - {target}
touch = {}
for k, (s, p) in enumerate(rx0):
    for n in set(s) | set(p): touch.setdefault(n, []).append(k)
dist, rd, frontier = {target: 0}, {}, [target]
while frontier:
    nxt = []
    for n in frontier:
        for k in touch.get(n, []):
            if k in rd: continue
            rd[k] = dist[n]
            s, p = rx0[k]
            for n2 in set(s) | set(p):
                if n2 not in dist and n2 not in hubs: dist[n2] = dist[n] + 1; nxt.append(n2)
    frontier = nxt
units = sorted(contractible, key=lambda e: -min(rd.get(k, 10**6) for k in mech[e]))
t0 = time.time()
print(f"{M} {target}: start {len(base)} reactions, {len(mech)} mechanisms, {len(contractible)} contractible; "
      f"gap at start: {gap(set())} ({time.time() - t0:.0f}s)", flush=True)
done, n = set(), 2
while len(units) >= 1:
    size = -(-len(units) // n); chunks = [units[i:i + size] for i in range(0, len(units), size)]
    accepted = False
    for ch in chunks:
        trial = done | set(ch)
        if gap(trial):
            done = trial; units = [u for u in units if u not in set(ch)]; n = max(n - 1, 2); accepted = True
            r_, m_ = save(done, "in progress")
            print(f"[{time.time() - t0:7.0f}s] contracted {len(ch)} mechanisms ({len(done)} total): {r_} reactions, {m_} metabolites", flush=True)
            break
    if not accepted:
        if n >= len(units): break
        n = min(len(units), 2 * n)
        print(f"[{time.time() - t0:7.0f}s] no block contractible, splitting finer: n = {n} of {len(units)}", flush=True)
r_, m_ = save(done, "1-minimal: contracting any further single mechanism loses the gap")
print(f"done: {len(done)} mechanisms contracted, {r_} reactions, {m_} metabolites ({time.time() - t0:.0f}s)", flush=True)

# stage 2: delete units from the contracted network
_, _, _, _, rx2 = build(done)
def ok_list(rx):
    try:
        A, Y, N, ms, _ = build_list(rx)
        if target not in ms: return False
        ti = ms.index(target); it = kimonu.iterated_route(crn_of(A, Y, N), detach_uptakes=False)
        if ti in set(it.acr): return False
        return ti in set(one_shot_exact(A, Y, N, False)[1])
    except Exception as e:
        print("   predicate error:", type(e).__name__, str(e)[:100], flush=True); return False
groups = {}
for k, (s_, p_) in enumerate(rx2):
    ids = {enzyme_of(n) for n in list(s_) + list(p_)} - {None}
    groups.setdefault(min(ids, key=int) if ids else f"r{k}", []).append(k)
def far(u):
    names = {n for k in u for n in list(rx2[k][0]) + list(rx2[k][1])}
    return -min((dist.get(n, 10**6) for n in names), default=10**6)
units2 = sorted(groups.values(), key=far)
keep = set(range(len(rx2))); n = 2
print(f"stage 2: {len(rx2)} reactions in {len(units2)} units", flush=True)
while len(units2) >= 2:
    size = -(-len(units2) // n); chunks = [units2[i:i + size] for i in range(0, len(units2), size)]
    accepted = False
    for ch in chunks:
        drop = {k for u in ch for k in u}; trial = sorted(keep - drop)
        if trial and ok_list([rx2[k] for k in trial]):
            _, _, _, _, pruned = build_list([rx2[k] for k in trial])
            live = {tuple(tuple(sorted(d.items())) for d in r) for r in pruned}
            keep = {k for k in trial if tuple(tuple(sorted(d.items())) for d in rx2[k]) in live}
            units2 = [u for u in units2 if set(u) <= keep]; n = max(n - 1, 2); accepted = True
            A, Y, N, ms, rx_now = build_list([rx2[k] for k in sorted(keep)])
            fmt = lambda d: "+".join(f"{'' if v == 1 else v}{nm_}" for nm_, v in sorted(d.items())) or "0"
            json.dump(dict(model=M, target=target, note="stage 2 in progress", n_reactions=len(rx_now), n_metabolites=len(ms),
                           reactions=[[fmt(a), fmt(b)] for a, b in rx_now],
                           structured=[[[[x, str(v)] for x, v in sorted(a.items())], [[x, str(v)] for x, v in sorted(b.items())]] for a, b in rx_now]),
                      open(out, "w"), indent=1)
            print(f"[{time.time() - t0:7.0f}s] deleted {len(ch)} units: {len(rx_now)} reactions, {len(ms)} metabolites", flush=True)
            break
    if not accepted:
        if n >= len(units2): break
        n = min(len(units2), 2 * n)
        print(f"[{time.time() - t0:7.0f}s] stage 2: no block deletable, splitting finer: n = {n} of {len(units2)}", flush=True)
print(f"stage 2 done: {len(keep)} reactions kept ({time.time() - t0:.0f}s)", flush=True)

