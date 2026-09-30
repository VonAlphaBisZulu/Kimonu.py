"""Integer-flux tangent readout against the previous one: identical results and timing, per model.

    python kinetic_variants/fast_vs_old.py <model>   -> kinetic_variants/fast_<model>.json
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.tangent_fast import one_shot_fast
M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc); S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
old = json.load(open(f"{B}/kinetic_variants/tangent_{M}.json"))
out = dict(model=M)
for mode, det in (("attached", False), ("detached", True)):
    t0 = time.time(); mods, acr, groups, dim = one_shot_fast(A, Y, S, det, want_groups=True); t = time.time() - t0
    o = old[mode]
    same = (len(acr) == o["acr_tangent"] and sorted((len(g) for g in groups), reverse=True) == o["acrr_groups_tangent"] and dim == o["tangent_dim"])
    out[mode] = dict(seconds=t, seconds_old=o["seconds_tangent"], same_counts=same)
    print(f"{M:28s} {mode:8s} {t:8.1f}s vs {o['seconds_tangent']:8.1f}s ({o['seconds_tangent'] / t:5.1f}x) | same ACR, groups, dim T: {same}", flush=True)
json.dump(out, open(f"{B}/kinetic_variants/fast_{M}.json", "w"), indent=1)
