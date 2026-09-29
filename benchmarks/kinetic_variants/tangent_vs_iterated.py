"""Jacobian route (tangent space T) against the iterated psi route, both uptake modes, per model.

Both are exact and use the same interior flux (seed 0). The iterated result is always sound, so every
coupling, ACR metabolite and ACRR pair it reports must also appear in T; anything only T reports is
"local only": fixed along each steady-state branch, not yet certified across branches.

    python kinetic_variants/tangent_vs_iterated.py <model> [...]   -> kinetic_variants/tangent_<model>.json
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import kimonu
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.tangent_one_shot import one_shot_exact

def pairs(groups):
    return {frozenset((a, b)) for g in groups for a in g for b in g if a < b}

for M in sys.argv[1:]:
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    zero = {c for c in range(nc) if Y.tocsc()[:, c].nnz == 0}
    nullity = kimonu.concordance_from_matrices(A, S).nullity
    out = dict(model=M, n_metabolites=nm, n_conservation=nm - (nr - nullity))
    for mode, det in (("attached", False), ("detached", True)):
        t0 = time.time(); it = kimonu.kinetic_modules(A, Y, S, detach_uptakes=det); t_it = time.time() - t0
        ar = kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr)
        t0 = time.time(); mods, acr, acrr_groups, dim = one_shot_exact(A, Y, S, det, want_groups=True); t_tan = time.time() - t0
        it_mods = [set(m) - zero for m in it.modules]
        coup_it, coup_tan = pairs([sorted(m) for m in it_mods]), pairs([sorted(set(m) - zero) for m in mods])
        acrr_it = {frozenset(p) for p in ar.acrr}; acrr_tan = pairs([sorted(g) for g in acrr_groups])
        out[mode] = dict(
            tangent_dim=dim, rounds=it.rounds, seconds_iterated=t_it, seconds_tangent=t_tan,
            couplings_iterated=len(coup_it), couplings_tangent=len(coup_tan),
            iterated_not_in_tangent=len(coup_it - coup_tan), local_only_couplings=len(coup_tan - coup_it),
            acr_iterated=len(ar.acr), acr_tangent=len(acr), local_only_acr=sorted(crn["metabolites"][i] for i in set(acr) - set(ar.acr)),
            acr_iterated_not_in_tangent=len(set(ar.acr) - set(acr)),
            acrr_pairs_iterated=len(acrr_it), acrr_pairs_tangent=len(acrr_tan),
            acrr_iterated_not_in_tangent=len(acrr_it - acrr_tan), local_only_acrr_pairs=len(acrr_tan - acrr_it),
            acrr_groups_tangent=sorted((len(g) for g in acrr_groups), reverse=True))
        o = out[mode]
        print(f"{M:28s} {mode:8s} T {dim:4d} (cons {out['n_conservation']:4d}) | couplings it {o['couplings_iterated']:7d} tan {o['couplings_tangent']:7d} local-only {o['local_only_couplings']:5d} missing {o['iterated_not_in_tangent']} | "
              f"ACR it {o['acr_iterated']:3d} tan {o['acr_tangent']:3d} | ACRR pairs it {o['acrr_pairs_iterated']:6d} tan {o['acrr_pairs_tangent']:6d} local-only {o['local_only_acrr_pairs']:5d} | {t_it:.0f}s vs {t_tan:.0f}s", flush=True)
    json.dump(out, open(f"{B}/kinetic_variants/tangent_{M}.json", "w"), indent=1)
