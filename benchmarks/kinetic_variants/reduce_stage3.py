"""Stage 3: shrink a witness further by fixing metabolites and merging linear chains.

Starting from a witness saved by reduce_contract.py (structured reactions), three moves are tried, each kept
only if the tangent space still fixes the target and the iterated route still does not (exact):
  1. fix metabolites as constant (drop them from every complex), in far-to-near blocks (ddmin);
  2. merge linear chains: a metabolite made by exactly one reaction and used by exactly one is eliminated
     by adding the two reactions in the ratio that cancels it;
  3. delete single reactions.
Rounds repeat until none of the moves applies. The target and the protected enzyme's forms are never fixed.

    python kinetic_variants/reduce_stage3.py <witness.json> <out.json> [protected enzyme id ...]
"""
import json, os, re, sys, time
from fractions import Fraction
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import kimonu
from kimonu.io import exact_matrix
from kinetic_variants.tangent_one_shot import one_shot_exact

src, out = sys.argv[1:3]; protect_ids = set(sys.argv[3:])
W = json.load(open(src)); target = W["target"]
rx = [({n: Fraction(v) for n, v in s}, {n: Fraction(v) for n, v in p}) for s, p in W["structured"]]
ENZ = re.compile(r"E_?(\d+)(?:_|$)")
protected = lambda n: n == target or (ENZ.match(n) and ENZ.match(n).group(1) in protect_ids)

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
    N = exact_matrix(((i, j, v) for (i, j), v in Ne.items() if v != 0), (len(ms), len(rx)))
    return A, Y, N, ms, rx

def ok(rx):
    try:
        A, Y, N, ms, _ = build(rx)
        if target not in ms: return False
        ti = ms.index(target); it = kimonu.kinetic_modules(A, Y, N)
        if ti in set(kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr).acr): return False
        return ti in set(one_shot_exact(A, Y, N, False)[1])
    except Exception as e:
        print("   predicate error:", type(e).__name__, str(e)[:100], flush=True); return False

def fix(rx, names):
    names = set(names); new = []
    for s, p in rx:
        s2 = {n: v for n, v in s.items() if n not in names}; p2 = {n: v for n, v in p.items() if n not in names}
        if s2 != p2: new.append((s2, p2))
    return new

def net(r):
    d = {}
    for n, v in r[0].items(): d[n] = d.get(n, 0) - v
    for n, v in r[1].items(): d[n] = d.get(n, 0) + v
    return d

def save(rx, note):
    _, _, _, ms, rx = build(rx)
    fmt = lambda d: "+".join(f"{'' if v == 1 else v}{n}" for n, v in sorted(d.items())) or "0"
    json.dump(dict(target=target, note=note, n_reactions=len(rx), n_metabolites=len(ms), reactions=[[fmt(s), fmt(p)] for s, p in rx],
                   structured=[[[[n, str(v)] for n, v in sorted(s.items())], [[n, str(v)] for n, v in sorted(p.items())]] for s, p in rx]),
              open(out, "w"), indent=1)
    return rx, len(ms)

def distance_order(rx):
    touch = {}
    for k, (s, p) in enumerate(rx):
        for n in set(s) | set(p): touch.setdefault(n, set()).add(k)
    deg = {n: len(v) for n, v in touch.items()}; hubs = set(sorted(deg, key=lambda n: -deg[n])[:10]) - {target}
    dist, frontier = {target: 0}, [target]
    while frontier:
        nxt = []
        for n in frontier:
            for k in touch.get(n, ()):
                for n2 in set(rx[k][0]) | set(rx[k][1]):
                    if n2 not in dist and n2 not in hubs: dist[n2] = dist[n] + 1; nxt.append(n2)
        frontier = nxt
    return sorted((n for n in touch if not protected(n)), key=lambda n: -dist.get(n, 10**6))

t0 = time.time()
rx, _ = save(rx, "stage 3 start")[0], None
print(f"start: {len(rx)} reactions; gap: {ok(rx)} ({time.time() - t0:.0f}s)", flush=True)
changed = True
while changed:
    changed = False
    # 1. fix metabolites, far-first blocks
    cands = distance_order(rx); n = 2
    while cands:
        size = -(-len(cands) // n); chunks = [cands[i:i + size] for i in range(0, len(cands), size)]
        hit = False
        for ch in chunks:
            trial = fix(rx, ch)
            if trial and ok(trial):
                rx, m_ = save(trial, "stage 3 in progress"); cands = [c for c in cands if c not in set(ch) and any(c in s or c in p for s, p in rx)]
                n = max(n - 1, 2); hit = changed = True
                print(f"[{time.time() - t0:6.0f}s] fixed {len(ch)} metabolites: {len(rx)} reactions, {m_} metabolites", flush=True)
                break
        if not hit:
            if n >= len(cands): break
            n = min(len(cands), 2 * n)
    # 2. merge linear chains
    progress = True
    while progress:
        progress = False
        use = {}
        for k, r in enumerate(rx):
            for m, v in net(r).items():
                if v != 0: use.setdefault(m, []).append((k, v))
        for m, lst in use.items():
            if protected(m) or len(lst) != 2: continue
            (k1, v1), (k2, v2) = lst
            if (v1 > 0) == (v2 > 0): continue
            a, b = abs(v2), abs(v1)                      # a*r1 + b*r2 cancels m
            d = {}
            for kk, f in ((k1, a), (k2, b)):
                for n_, v in net(rx[kk]).items(): d[n_] = d.get(n_, 0) + f * v
            d = {n_: v for n_, v in d.items() if v != 0}
            merged = ({n_: -v for n_, v in d.items() if v < 0}, {n_: v for n_, v in d.items() if v > 0})
            trial = [r for kk, r in enumerate(rx) if kk not in (k1, k2)] + ([merged] if merged[0] != merged[1] else [])
            if ok(trial):
                rx, m_ = save(trial, "stage 3 in progress"); progress = changed = True
                print(f"[{time.time() - t0:6.0f}s] merged chain through {m}: {len(rx)} reactions, {m_} metabolites", flush=True)
                break
    # 3. delete single reactions
    for k in range(len(rx) - 1, -1, -1):
        if k >= len(rx): continue
        trial = rx[:k] + rx[k + 1:]
        if ok(trial):
            rx, m_ = save(trial, "stage 3 in progress"); changed = True
            print(f"[{time.time() - t0:6.0f}s] deleted a reaction: {len(rx)} reactions, {m_} metabolites", flush=True)
rx, m_ = save(rx, "stage 3 done: no fix, merge or deletion keeps the gap")
print(f"stage 3 done: {len(rx)} reactions, {m_} metabolites ({time.time() - t0:.0f}s)", flush=True)
