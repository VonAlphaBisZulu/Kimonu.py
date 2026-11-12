"""How much of each network the relations cover: complexes inside kinetic modules and metabolites with ACR or ACRR.

uptake_modes.py counts pairs (coupled complex pairs, ACRR pairs). This script counts members, for the two shares
of Figure 3: the complexes that lie in a kinetic module of two or more complexes, and the metabolites that have
ACR or belong to an ACRR group, over all metabolites and over the free metabolites (those of the original model,
without enzymes and enzyme-substrate complexes). A free metabolite counts as an ACRR member when its group holds at
least one other free metabolite, in line with the ACRR pairs among free metabolites of uptake_modes.py.

The iterated route (kimonu.iterated_route, from src) runs on every model. The tangent route runs only
where the uptake_modes record reports a gap, since elsewhere both routes return the same modules, ACR and ACRR. The
counts of this run must reproduce the record of uptake_modes.py (modules, giant, pairs, acr, acrr_pairs and the
two free counts); that is the check of this script.

    python kinetic_variants/member_counts.py <model> [attached]   -> kinetic_variants/member_counts/<model>.json
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.uptake_modes import groups_of

KEYS = ("modules", "giant", "pairs", "acr", "acrr_pairs", "acr_free", "acrr_pairs_free")


def members(modules, acr, acrr, free):
    sizes = sorted((len(m) for m in modules if len(m) > 1), reverse=True)
    acr = set(acr); pairs = {tuple(sorted(p)) for p in acrr}
    groups = groups_of(pairs)[0]
    free_groups = [g & free for g in groups if len(g & free) > 1]
    in_groups = set().union(*groups) if groups else set()
    in_free_groups = set().union(*free_groups) if free_groups else set()
    record = dict(modules=len(sizes), giant=max(sizes, default=1), pairs=sum(s * (s - 1) // 2 for s in sizes),
                  acr=len(acr), acrr_pairs=len(pairs), acr_free=len(acr & free),
                  acrr_pairs_free=sum(a in free and b in free for a, b in pairs))
    return record, dict(
        module_sizes=sizes, complexes_in_modules=sum(sizes),
        acr=len(acr), acr_free=len(acr & free),
        acrr_group_sizes=sorted((len(g) for g in groups), reverse=True), acrr_members=len(in_groups),
        acrr_free_group_sizes=sorted((len(g) for g in free_groups), reverse=True), acrr_members_free=len(in_free_groups),
        fixed_members=len(acr | in_groups), fixed_members_free=len((acr & free) | in_free_groups))


if __name__ == "__main__":
    M = sys.argv[1]; name = "attached" if "attached" in sys.argv[2:] else "detached"; detach = name == "detached"
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
    free = set(range(mets.index("E_1") if "E_1" in mets else nm))
    f = f"{B}/kinetic_variants/uptake_modes/{M}.json"; one = f.replace(".json", f".{name}.json")
    stored = json.load(open(one if os.path.exists(one) else f))[name]
    out = dict(model=M, uptakes=name, n_complexes=nc, n_metabolites=nm, n_free=len(free), gap=stored["gap"])
    t0 = time.time()
    it = kimonu.iterated_route(A, Y, N, detach_uptakes=detach)
    record, out["it"] = members(it.modules, it.acr, it.acrr, free)
    assert record == {k: stored["it"][k] for k in KEYS}, (M, "iterated", record, stored["it"])
    out["rounds"], out["seconds_iterated"] = it.rounds, time.time() - t0
    out["tan"] = None                             # no gap: the tangent route returns the same lists
    if stored["gap"]:
        t0 = time.time()
        tan = kimonu.tangent_route(A, Y, N, detach_uptakes=detach)
        record, out["tan"] = members(tan.modules, tan.acr, tan.acrr, free)
        assert record == {k: stored["tan"][k] for k in KEYS}, (M, "tangent", record, stored["tan"])
        out["seconds_tangent"] = time.time() - t0
    else:
        assert all(stored["it"][k] == stored["tan"][k] for k in KEYS), M
    os.makedirs(f"{B}/kinetic_variants/member_counts", exist_ok=True)
    json.dump(out, open(f"{B}/kinetic_variants/member_counts/{M}{'' if detach else '.attached'}.json", "w"), indent=1)
    print(M, name, json.dumps({k: v for k, v in out.items() if k not in ("it", "tan")}),
          {k: v for k, v in out["it"].items() if "sizes" not in k}, flush=True)
