"""The iterated route with the Descartes step against the tangent space, per model and uptake mode.

    python kinetic_variants/descartes_vs_tangent.py <model> [...]   -> kinetic_variants/descartes_<model>.json
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import kimonu
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.tangent_one_shot import one_shot_exact
from kinetic_variants.descartes import kinetic_modules_descartes

def pairs(groups): return {frozenset((a, b)) for g in groups for a in g for b in g if a < b}
for M in sys.argv[1:]:
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc); S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    zero = {c for c in range(nc) if Y.tocsc()[:, c].nnz == 0}; mets = crn["metabolites"]
    out = dict(model=M)
    for mode, det in (("attached", False), ("detached", True)):
        t0 = time.time(); mods_d, acr_d, groups_d, fire = kinetic_modules_descartes(A, Y, S, detach_uptakes=det); t_d = time.time() - t0
        t0 = time.time(); mods_t, acr_t, groups_t, dim = one_shot_exact(A, Y, S, det, want_groups=True); t_t = time.time() - t0
        cd = pairs([sorted(set(m) - zero) for m in mods_d]); ct = pairs([sorted(set(m) - zero) for m in mods_t])
        ad, at = set(acr_d), set(acr_t); pd, pt = pairs([sorted(g) for g in groups_d]), pairs([sorted(g) for g in groups_t])
        out[mode] = dict(firings=[mets[i] for i in fire], seconds_descartes=t_d, seconds_tangent=t_t,
                         local_only_couplings=len(ct - cd), local_only_acr=sorted(mets[i] for i in at - ad), local_only_acrr_pairs=len(pt - pd),
                         descartes_not_in_tangent=dict(couplings=len(cd - ct), acr=len(ad - at), acrr=len(pd - pt)))
        o = out[mode]
        print(f"{M:12s} {mode:8s} firings {len(fire):3d} | local-only now: couplings {o['local_only_couplings']:5d} ACR {o['local_only_acr']} ACRR pairs {o['local_only_acrr_pairs']:5d} | "
              f"not in tangent (must be 0): {o['descartes_not_in_tangent']} | {t_d:.0f}s vs {t_t:.0f}s", flush=True)
    json.dump(out, open(f"{B}/kinetic_variants/descartes_{M}.json", "w"), indent=1)
