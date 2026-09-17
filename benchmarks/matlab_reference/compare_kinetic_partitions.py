"""Kinetic partitions of the three implementations on iCN900, on an exact complex bijection.

COCOA writes a stoichiometric coefficient as a `2_` prefix on the species, which earlier
comparisons here dropped, silently discarding the ~54 multi-coefficient complexes -- exactly the
boundary complexes the partitions disagree about. Rendering both sides the same way removes that.
"""
import json, sys
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks")
sys.path.insert(0, "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/src")
import kimonu
from kimonu.kinetic import _Topology, _sccs
from run_kimonu_crn import triplets_to_sparse, to_exact

B = "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks"
M = "iCN900"
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64)
Y = triplets_to_sparse(crn["Y"], nm, nc)
S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
topo = _Topology(A); succ, pred = topo.succ, topo.pred
term = set()
for comp in _sccs(succ, set(range(nc))):
    if not any(succ[v] - comp for v in comp):
        term |= comp

kim = json.load(open(f"{B}/sweep_shared/kimonu_{M}.json"))["complexes"]
def cocoa_name(comp):
    """COCOA's complex id: species sorted, coefficient as a prefix -- `2_h[c]`, `4.5_h[c]`, none for 1."""
    parts = []
    for sp, v in comp.items():
        v = float(v)
        parts.append(sp if v == 1 else f"{int(v) if v.is_integer() else v}_{sp}")
    return "+".join(sorted(parts))

ours = [cocoa_name(c["composition"]) for c in kim]
co = json.load(open(f"{B}/sweep/cocoa_eff/cocoa_eff_{M}_HiGHS_p64_seed1234.json"))
ids = ["+".join(sorted(s.split("+"))) for s in co["complex_ids"]]
assert len(ids) == len(ours)
dup = [n for n, c in Counter(ours).items() if c > 1]
print(f"bijection: unmatched cocoa {len(set(ids)-set(ours))}  unmatched kimonu {len(set(ours)-set(ids))}"
      f"  colliding names {len(dup)}")
pos = {n: i for i, n in enumerate(ours)}
C = defaultdict(set)
for i, m in enumerate(co["kinetic_module"]):
    if ids[i] in pos: C[m].add(pos[ids[i]])

non_reactant = {c for c in range(nc) if not succ[c]}
def tentative(C0):
    C1 = set(C0)
    while True:
        e = {c for c in C1 if pred[c] - C1}
        if not e: break
        C1 -= e
        if not C1: return set()
    C2 = C1 - non_reactant
    C3, Cm3 = set(C2), set()
    while True:
        x = {c for c in C3 if succ[c] - C3}
        if not x: break
        C3 -= x; Cm3 |= x
        if not C3: break
    Cm4 = set()
    for comp in _sccs({v: succ[v] & C3 for v in C3}, set(C3)):
        if any((succ[v] & C3) - comp for v in comp): Cm4 |= comp
    return Cm3 | Cm4
def merge(sets):
    out = []
    for g in sets:
        hit = [x for x in out if x & g]
        for x in hit: out.remove(x); g = g | x
        out.append(g)
    return out

res = kimonu.concordance_from_matrices(A, S); bal = set(res.balanced)
seeds = [m for m in (set(x) - bal for x in res.modules) if m]
T = merge([s for s in (tentative(m | bal) for m in seeds) if s])
kin = kimonu.kinetic_modules(A, res.modules, res.balanced, Y=Y, nullity=res.nullity)
for tag, blocks in (("tentative (S.2.2)", T), ("Kimonu (S2-1)", [set(m) for m in kin.modules]),
                    ("COCOA", list(C.values()))):
    ns = [b for b in blocks if len(b) > 1]
    cov = set().union(*ns) if ns else set()
    print(f"{tag:20s} non-singleton {len(ns):4d}  largest {max((len(b) for b in ns), default=0):4d}"
          f"  covered {len(cov):4d}  globally terminal inside {len(cov & term):4d}")
