"""Pairs an exhaustive reference run rejected that Kimonu accepts, as 1-based rows for retest_rejected.m.

    python matlab_reference/retest_pairs.py <model>   ->  reference_final/retest/pairs_<model>.txt
"""
import json, sys, os
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, B)
from compare import partition, key, concordant_pairs
from run_kimonu_crn import triplets_to_sparse
m = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{m}_fixed.crn.json")); nc, nm = crn["n_complexes"], crn["n_metabolites"]
Y = triplets_to_sparse(crn["Y"], nm, nc).tocsc()
idx = {key({crn["metabolites"][i]: float(Y[i, c]) for i in Y[:, c].nonzero()[0]}): c for c in range(nc)}
r = json.load(open(f"{B}/reference_final/reference_{m}.json")); k = json.load(open(f"{B}/sweep_shared/kimonu_{m}.json"))
rp, kp = concordant_pairs(partition(r)[1]), concordant_pairs(partition(k)[1])
extra = rp - kp
rows = sorted((min(idx[a], idx[b]) + 1, max(idx[a], idx[b]) + 1) for a, b in kp - rp)
out = f"{B}/reference_final/retest/pairs_{m}.txt"
open(out, "w").write("".join(f"{a} {b}\n" for a, b in rows))
print(f"{m}: reference pairs {len(rp)}, Kimonu {len(kp)}, rejected-by-reference {len(rows)}, accepted-not-in-Kimonu {len(extra)} -> {out}")
