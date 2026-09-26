"""The tentative four-phase algorithm of S.2.2, implemented against Kimonu's topology.

If Kimonu reproduces the reference's output when told to run the reference's algorithm, then both
implementations are correct and the only thing separating them is which of the paper's two
algorithms they chose. That is a different and much stronger statement than either being buggy.

  Phase I   C1 = C0 minus entry complexes, to a fixed point
  Phase II  C2 = C1 minus non-reactant complexes (no outgoing reaction anywhere)
  Phase III C3 = C2 minus exit complexes, to a fixed point; the removed ones form C_-3
  Phase IV  C4 = C3 minus its nonterminal complexes, judged in G_{C4}; removed ones form C_-4
  Final     C_up = C_-4 u C_-3
"""
import json, sys
import numpy as np
sys.path.insert(0, "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks")
sys.path.insert(0, "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/src")
import kimonu
from compat import upstream
from compat.upstream import _Topology, _sccs
from run_kimonu_crn import triplets_to_sparse, to_exact

B = "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks"
M = sys.argv[1] if len(sys.argv) > 1 else "iCN900"
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64)
Y = triplets_to_sparse(crn["Y"], nm, nc)
S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
topo = _Topology(A); succ, pred = topo.succ, topo.pred
non_reactant = {c for c in range(nc) if not succ[c]}

def tentative(C0):
    C1 = set(C0)
    while True:                                        # Phase I
        entry = {c for c in C1 if pred[c] - C1}
        if not entry: break
        C1 -= entry
        if not C1: return set()
    C2 = C1 - non_reactant                             # Phase II
    C3, Cm3 = set(C2), set()                           # Phase III
    while True:
        exit_ = {c for c in C3 if succ[c] - C3}
        if not exit_: break
        C3 -= exit_; Cm3 |= exit_
        if not C3: break
    Cm4 = set()                                        # Phase IV, terminality inside G_{C4}
    if C3:
        for comp in _sccs({v: succ[v] & C3 for v in C3}, set(C3)):
            if any((succ[v] & C3) - comp for v in comp):
                Cm4 |= comp                            # nonterminal in the induced graph
    return Cm3 | Cm4

res = kimonu.concordance_from_matrices(A, S)
bal = set(res.balanced)
seeds = [set(m) - bal for m in res.modules]
seeds = [m for m in seeds if m]
print(f"{M}: {len(seeds)} seeds, {len(bal)} balanced")

sets = [s for s in (tentative(m | bal) for m in seeds) if s]
tot = sum(len(s) for s in sets); uni = len(set().union(*sets)) if sets else 0
print(f"tentative (S.2.2)   sets {len(sets):4d}  total {tot:5d}  union {uni:5d}  overlap {tot-uni:5d}")
# merge overlapping, as both implementations do
merged = []
for g in sets:
    hit = [x for x in merged if x & g]
    for x in hit: merged.remove(x); g = g | x
    merged.append(g)
ns = [len(x) for x in merged if len(x) > 1]
print(f"  after overlap merge: blocks {len(merged):4d}  non-singleton {len(ns):4d}  "
      f"largest {max((len(x) for x in merged), default=0):4d}  covered {sum(ns):5d}")
term = {c for c in range(nc) for comp in [None] }  # placeholder
allc = _sccs(succ, set(range(nc)))
term = set()
for comp in allc:
    if not any(succ[v] - comp for v in comp): term |= comp
bad = sum(1 for x in merged for c in x if c in term)
print(f"  complexes in the result that are globally terminal (Lemma S2-4 violations): {bad}")
kin = upstream.kinetic_modules(A, res.modules, res.balanced, Y=Y, nullity=res.nullity)
sz = sorted((len(m) for m in kin.modules), reverse=True); k_ns = [x for x in sz if x > 1]
print(f"Upstream Alg (S2-1) blocks {len(sz):4d}  non-singleton {len(k_ns):4d}  largest {sz[0]:4d}  covered {sum(k_ns):5d}")
print(f"reference R stage   (corrected NA)                    non-singleton  246  largest   19  covered   583")

# direct check: what separates the two algorithms' outputs, and is any of it terminal?
tent = set().union(*merged) if merged else set()
kim_cov = set()
for m in kin.modules:
    if len(m) > 1: kim_cov |= set(m)
only_t = tent - kim_cov; only_k = kim_cov - tent
print(f"\ncoverage: tentative {len(tent)}  upstream-alg {len(kim_cov)}")
print(f"  only in tentative: {len(only_t)}  of which globally terminal {len(only_t & term)}")
print(f"  only in upstream : {len(only_k)}  of which globally terminal {len(only_k & term)}")
