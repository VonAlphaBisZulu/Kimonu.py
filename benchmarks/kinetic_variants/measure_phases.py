"""Four-phase (S.2.2) vs two-phase (Remark S2-1) upstream sets, seed by seed, on Kimonu's topology.

Per model: (a) per-seed set differences between the two algorithms and how each differing complex
is classified (C_-3 exit / C_-4 / dropped as terminal / dropped as entry), (b) property checks on
the final modules of both algorithms with balanced-only linkage classes excluded: autonomy (no
producer outside the module), feeding (every complex reaches a complex outside the module in the
full graph), Lemma S2-4 (no globally terminal complex).
"""
import json, sys
import numpy as np
sys.path.insert(0, "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks")
sys.path.insert(0, "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/src")
import kimonu
from compat import upstream
from compat.upstream import _Topology, _sccs, upstream_set, _weak_components
from run_kimonu_crn import triplets_to_sparse, to_exact
B = "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks"

def load(M):
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    return crn, nc, A, Y, S

def tentative(C0, succ, pred, non_reactant):
    C1 = set(C0)
    while True:
        entry = {c for c in C1 if pred[c] - C1}
        if not entry: break
        C1 -= entry
        if not C1: return set(), set(), set(), set()
    C2 = C1 - non_reactant
    C3, Cm3 = set(C2), set()
    while True:
        exit_ = {c for c in C3 if succ[c] - C3}
        if not exit_: break
        C3 -= exit_; Cm3 |= exit_
        if not C3: break
    Cm4 = set()
    if C3:
        for comp in _sccs({v: succ[v] & C3 for v in C3}, set(C3)):
            if any((succ[v] & C3) - comp for v in comp): Cm4 |= comp
    return Cm3 | Cm4, C1, Cm3, Cm4

def reach_outside(module, succ):
    """complexes of `module` from which some complex outside `module` is reachable"""
    out = {c for c in module if succ[c] - module}
    changed = True
    while changed:
        changed = False
        for c in module - out:
            if succ[c] & out: out.add(c); changed = True
    return out

def props(mods, succ, pred, term, balanced):
    """counts over non-singleton modules that are not balanced-only linkage classes"""
    n = aut = feed = lem = 0; lem_bal = 0
    for m in mods:
        if len(m) < 2 or m <= balanced: continue
        n += 1
        if any(pred[c] - m for c in m): aut += 1
        if len(reach_outside(m, succ)) < len(m): feed += 1
        t = m & term
        if t: lem += 1; lem_bal += int(t <= balanced)
    return n, aut, feed, lem, lem_bal

def run(M, detail=False):
    crn, nc, A, Y, S = load(M)
    topo = _Topology(A); succ, pred = topo.succ, topo.pred
    non_reactant = {c for c in range(nc) if not succ[c]}
    term = set()
    for comp in _sccs(succ, set(range(nc))):
        if not any(succ[v] - comp for v in comp): term |= comp
    res = kimonu.concordance_from_matrices(A, S); bal = set(res.balanced)
    seeds = [set(m) - bal for m in res.modules]; seeds = [m for m in seeds if m]
    only_t = only_k = 0; kinds = {}
    tsets = []
    for sd in seeds:
        C0 = sd | bal
        t, C1, Cm3, Cm4 = tentative(C0, succ, pred, non_reactant)
        k = upstream_set(C0, topo)
        if t: tsets.append(t)
        for c in t - k:
            only_t += 1
            where = "C-3 exit" if c in Cm3 else "C-4"
            # why Kimonu dropped it: entry cascade cannot differ (same phase I); so terminal SCC of C1 in full graph
            comps = _sccs(succ, C1); comp = next(x for x in comps if c in x)
            why = "in terminal SCC of C1 (full graph)" if not any(succ[v] - comp for v in comp) else "??"
            key = (where, why, "globally terminal" if c in term else "nonterminal"); kinds[key] = kinds.get(key, 0) + 1
            if detail: print(f"   only 4-phase: {crn_name(crn, c)}  [{where}; {why}; {'TERM' if c in term else 'nonterm'}]")
        for c in k - t:
            only_k += 1
            key = ("only 2-phase", "non-reactant" if c in non_reactant else "other", "globally terminal" if c in term else "nonterminal"); kinds[key] = kinds.get(key, 0) + 1
            if detail: print(f"   only 2-phase: {crn_name(crn, c)}  [{key[1]}; {key[2]}]")
    # final modules: tentative with overlap merge; Kimonu's own
    merged = []
    for g in tsets:
        hit = [x for x in merged if x & g]
        for x in hit: merged.remove(x); g = g | x
        merged.append(g)
    kin = upstream.kinetic_modules(A, res.modules, res.balanced, Y=Y, nullity=res.nullity)
    pt = props(merged, succ, pred, term, bal); pk = props(kin.modules, succ, pred, term, bal)
    gt = max((len(x) for x in merged), default=0); gk = kin.giant_module
    print(f"{M:26s} seeds {len(seeds):5d} | only4 {only_t:4d} only2 {only_k:4d} | 4-phase mods {pt[0]:4d} aut {pt[1]} feed {pt[2]} S2-4 {pt[3]} (bal {pt[4]}) giant {gt:4d} | 2-phase mods {pk[0]:4d} aut {pk[1]} feed {pk[2]} S2-4 {pk[3]} (bal {pk[4]}) giant {gk:4d}")
    for k_, v in sorted(kinds.items()): print(f"     {v:5d}  {k_}")
    sys.stdout.flush()

def crn_name(crn, c):
    # composition of complex c from Y triplets
    comp = [(crn["metabolites"][i], v) for i, j, v in crn["Y"] if j == c] if isinstance(crn["Y"][0], list) else []
    return " + ".join(f"{v:g} {m}" if v != 1 else m for m, v in comp) or f"complex#{c}"

if __name__ == "__main__":
    for M in sys.argv[1:]:
        run(M, detail=(len(sys.argv) == 2))
