"""Pairs an exhaustive reference run rejected that Kimonu accepts, as 1-based rows for retest_rejected.m.

    python matlab_reference/retest_pairs.py <model>   ->  retest/pairs_<model>.txt in results/reference/exhaustive/

The reference run and Kimonu's September run (upstream_set_route_2026-09/kimonu_runs/) are read from KIMONU_RECORDS,
the pairs written to KIMONU_OUT (records.py).
"""
import json, sys, os
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, B)
import records                                         # paths of the records (KIMONU_RECORDS, KIMONU_OUT)
from compare import partition, key, concordant_pairs
from run_kimonu_crn import triplets_to_sparse
m = sys.argv[1]
crn = json.load(open(records.network(m))); nc, nm = crn["n_complexes"], crn["n_metabolites"]
Y = triplets_to_sparse(crn["Y"], nm, nc).tocsc()
idx = {key({crn["metabolites"][i]: float(Y[i, c]) for i in Y[:, c].nonzero()[0]}): c for c in range(nc)}
r = json.load(open(records.read(records.REF_EXHAUSTIVE, f"reference_{m}.json"))); k = json.load(open(records.read(records.UPSTREAM_ROUTE, "kimonu_runs", f"kimonu_{m}.json")))
rp, kp = concordant_pairs(partition(r)[1]), concordant_pairs(partition(k)[1])
extra = rp - kp
rows = sorted((min(idx[a], idx[b]) + 1, max(idx[a], idx[b]) + 1) for a, b in kp - rp)
out = records.write(records.REF_EXHAUSTIVE, "retest", f"pairs_{m}.txt")
open(out, "w").write("".join(f"{a} {b}\n" for a, b in rows))
print(f"{m}: reference pairs {len(rp)}, Kimonu {len(kp)}, rejected-by-reference {len(rows)}, accepted-not-in-Kimonu {len(extra)} -> {out}")
