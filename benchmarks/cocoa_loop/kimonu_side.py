"""Kimonu's half of the loop comparison: exact concordance, and the iterated route in both uptake modes.

The concordance modules written here are what cocoa_side.jl hands to COCOA's own kinetic stage, so the
two kinetic routes start from the same exact partition and only the loop differs. Complexes are keyed by
composition (sorted metabolite:coefficient pairs), which is how the two tools are matched; the zero complex
has no composition and is left out, as COCOA leaves it out of every result.

    python kimonu_side.py <model> <out.json> [crn.json]   (default: benchmarks/unblocked/<model>_fixed.crn.json)
"""
import json, os, sys, time
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import numpy as np
import kimonu
from kimonu.io import substrate_complexes
from run_kimonu_crn import triplets_to_sparse, to_exact


def key(comp):
    return "|".join(f"{m}:{c:.9g}" for m, c in sorted(comp.items()))


def main(M, out, crn_path=None):
    crn = json.load(open(crn_path or f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    mets = crn["metabolites"]
    Yc = Y.tocsc()
    comp = [{mets[int(i)]: float(v) for i, v in zip(Yc[:, j].indices, Yc[:, j].data)} for j in range(nc)]
    keep = [j for j in range(nc) if comp[j]]                     # the zero complex has no composition
    t0 = time.perf_counter()
    conc = kimonu.concordance_from_matrices(A, S)
    res = {"model": M, "seconds_concordance": time.perf_counter() - t0,
           "balanced": [key(comp[j]) for j in conc.balanced if comp[j]],
           "concordance_modules": [[key(comp[j]) for j in g if comp[j]] for g in conc.modules],
           "substrate": sorted(key(comp[c]) for c in set(substrate_complexes(A).values()) if comp[c])}
    res["concordance_modules"] = [g for g in res["concordance_modules"] if g]
    for mode, det in (("attached", False), ("detached", True)):
        t0 = time.perf_counter()
        kin = kimonu.iterated_route(A, Y, S, detach_uptakes=det, metabolite_ids=mets)
        res[mode] = {"seconds": time.perf_counter() - t0, "rounds": kin.rounds,
                     "modules": [[key(comp[j]) for j in g if comp[j]] for g in kin.modules if len([j for j in g if comp[j]]) > 1],
                     "acr": sorted(kin.acr), "acrr": sorted(sorted(p) for p in kin.acrr)}
    json.dump(res, open(out, "w"))
    print(M, "done", {m: len(res[m]["modules"]) for m in ("attached", "detached")}, flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
