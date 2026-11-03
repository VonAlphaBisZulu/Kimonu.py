"""Seed invariance of the iterated route: kinetic modules, ACR and ACRR at three random interior fluxes r*.

The relations hold for generic rate constants, so a different seed (a different r*) may change them only with
probability zero. Runs kimonu.kinetic_modules + acr_acrr with detached uptakes at seeds 0, 1, 2 and reports
whether the module partition, the ACR metabolites and the ACRR groups are identical.

    python kinetic_variants/seed_invariance.py <model>   -> $SEED_OUT/<model>.json
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from run_kimonu_crn import triplets_to_sparse, to_exact

OUT = os.environ.get("SEED_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_seeds")
M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
runs = []
for seed in (0, 1, 2):
    t = time.time()
    kin = kimonu.kinetic_modules(A, Y, N, detach_uptakes=True, seed=seed)
    ar = kimonu.acr_acrr(Y, kin.modules, known_acr=kin.known_acr)
    runs.append(dict(seed=seed, seconds=time.time() - t,
                     modules=sorted(sorted(map(int, m)) for m in kin.modules if len(m) > 1),
                     acr=sorted(map(int, ar.acr)), acrr=sorted(sorted(map(int, p)) for p in ar.acrr)))
    print(M, seed, len(runs[-1]["modules"]), len(runs[-1]["acr"]), len(runs[-1]["acrr"]), flush=True)
same = {k: all(r[k] == runs[0][k] for r in runs) for k in ("modules", "acr", "acrr")}
os.makedirs(OUT, exist_ok=True)
json.dump(dict(model=M, identical=same, runs=[{k: v for k, v in r.items() if k in ("seed", "seconds")} | {
    "n_modules": len(r["modules"]), "n_acr": len(r["acr"]), "n_acrr": len(r["acrr"])} for r in runs]), open(f"{OUT}/{M}.json", "w"), indent=1)
print(M, "identical", same)
