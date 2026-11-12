"""Kimonu in reference-equivalent mode (kimonu.reference_route), timed stage by stage.

The stages of kimonu.reference_route run one by one so that each is timed (fig. S1). Balanced complexes and
concordance from ker N (M = A K). Kinetic modules from upstream sets seeded by the concordance modules with two or
more complexes, together with the balanced complexes, and merged where they overlap, as the reference's R stage does
(kimonu.reference.upstream_modules). ACR and ACRR by the pairwise Shinar-Feinberg rule within kinetic modules
(kimonu.reference.shinar_feinberg). Timed end to end, blocked reactions already removed.

    python run_kimonu_refmode.py <model> <out.json>   # modules: the non-singleton kinetic modules, complex indices
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import kimonu
from kimonu.network import zero_complexes
from kimonu.reference import shinar_feinberg, upstream_modules
from run_kimonu_crn import triplets_to_sparse

M, out = sys.argv[1:3]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
net = kimonu.crn_from_matrices(A, Y)              # N = Y·A, the S of the crn.json
t0 = time.perf_counter()
res = kimonu.concordance(net)
t1 = time.perf_counter()
modules, _ = upstream_modules(net, [g for g in res.modules if len(g) > 1], res.balanced)
t2 = time.perf_counter()
acr, acrr = shinar_feinberg(net, modules)
t3 = time.perf_counter()
mets = crn["metabolites"]; cut = mets.index("E_1") if "E_1" in mets else nm
zero = zero_complexes(net)
n_kinetic, giant = sum(1 for m in modules if len(m) > 1), max(len(m) for m in modules)
json.dump(dict(model=M, seconds=t3 - t0, seconds_concordance=t1 - t0, seconds_kinetic=t2 - t1, seconds_acr=t3 - t2,
               n_complexes=nc - len(zero), balanced=res.n_balanced, conc=sum(1 for g in res.modules if len(g) > 1),
               n_kinetic=n_kinetic, giant=giant,
               acr=len(acr), acrr=len(acrr), acr_free=sum(1 for i in acr if i < cut),
               acrr_free=sum(1 for i, j in acrr if i < cut and j < cut),
               modules=[sorted(map(int, m)) for m in modules if len(m) > 1]), open(out, "w"), indent=1)
print(f"{M}: {t3 - t0:.2f} s, kinetic modules {n_kinetic}, giant {giant}")
