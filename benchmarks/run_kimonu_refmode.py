"""Kimonu in reference-equivalent mode: the reference's definitions, computed from null spaces.

Balanced complexes and concordance from ker N (M = A K); kinetic modules from upstream sets seeded by
the multi-complex concordance classes only and merged where they overlap, as the reference's R stage
does (compat/upstream.py with Y=None: no Proposition S3-4 or deficiency merge, which the reference does
not apply); ACR/ACRR by the reference's combinatorial rule (Corollary S2-10). Timed end to end, blocked
reactions already removed.

    python run_kimonu_refmode.py <model> <out.json>
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import kimonu
from compat import upstream
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.s210_criterion import s210

M, out = sys.argv[1:3]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc); S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
t0 = time.perf_counter()
res = kimonu.concordance_from_matrices(A, S)
t1 = time.perf_counter()
kin = upstream.kinetic_modules(A, [g for g in res.modules if len(g) > 1], res.balanced, Y=None)
t2 = time.perf_counter()
acr, acrr = s210(Y, kin.modules)
t3 = time.perf_counter()
mets = crn["metabolites"]; cut = mets.index("E_1") if "E_1" in mets else nm
zero = {c for c in range(nc) if Y.tocsc()[:, c].nnz == 0}
json.dump(dict(model=M, seconds=t3 - t0, seconds_concordance=t1 - t0, seconds_kinetic=t2 - t1, seconds_acr=t3 - t2,
               n_complexes=nc - len(zero), balanced=res.n_balanced, conc=sum(1 for g in res.modules if len(g) > 1),
               n_kinetic=sum(1 for m in kin.modules if len(m) > 1), giant=kin.giant_module,
               acr=len(acr), acrr=len(acrr), acr_free=sum(1 for i in acr if i < cut),
               acrr_free=sum(1 for i, j in acrr if i < cut and j < cut)), open(out, "w"), indent=1)
print(f"{M}: {t3 - t0:.2f} s, kinetic modules {sum(1 for m in kin.modules if len(m) > 1)}, giant {kin.giant_module}")
