"""Attached against detached uptakes: both routes, the gap relations and their dimension, per model.

For each uptake mode the iterated route (kimonu.iterated_route) and the tangent route
(kimonu.tangent_route), both from src, are compared with the counts of hierarchy_data.py (coupled complex
pairs over all complexes, ACR metabolites, ACRR pairs, and the same restricted to free metabolites). The
gap relations are listed as in two_point_gap.py prepare: every ACR metabolite only the tangent route
reports, and for each tangent ACRR group one link per iterated group it joins (rep_1 - rep_k). g_rel is the
dimension these relation vectors add to C_iter, the span of the iterated route's fixed log-combinations
(exact ranks, as in injectivity_gap.py). The attached block reproduces the October records hierarchy/<model>.json (commit 7fafd4a)
and, on the gap models, the 31/6/14 relations and 21/3/7 directions; that is the check of this script.

    python kinetic_variants/uptake_modes.py <model>            -> results/kimonu/uptake_modes/<model>.json
    python kinetic_variants/uptake_modes.py <model> detached   -> results/kimonu/uptake_modes/<model>.detached.json

Written to KIMONU_OUT (records.py).
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import records                                         # paths of the records (KIMONU_RECORDS, KIMONU_OUT)
import kimonu
from run_kimonu_crn import triplets_to_sparse
from kinetic_variants.injectivity_gap import exact_V, kernel_dim
from fractions import Fraction


def pairs(groups):
    return {(a, b) for g in groups for a in g for b in g if a < b}


def groups_of(links):
    """ACRR groups (sets of metabolites) from their pairs, by union-find"""
    parent = {}
    def find(a):
        while parent.setdefault(a, a) != a:
            a = parent[a]
        return a
    for a, b in links:
        parent[find(a)] = find(b)
    out = {}
    for a in parent:
        out.setdefault(find(a), set()).add(a)
    return list(out.values()), find, parent


def mode(crn, net, detach, free):
    nm = len(net.metabolite_ids)
    t0 = time.time()
    it = kimonu.iterated_route(net, detach_uptakes=detach)
    t_it = time.time() - t0
    t0 = time.time()
    tan = kimonu.tangent_route(net, detach_uptakes=detach)
    t_tan = time.time() - t0
    cp_it, cp_t = pairs([sorted(m) for m in it.modules]), pairs([sorted(m) for m in tan.modules])
    acr_it, acr_t = set(it.acr), set(tan.acr)
    acrr_it = {tuple(sorted(p)) for p in it.acrr}; acrr_t = {tuple(sorted(p)) for p in tan.acrr}
    nontriv = lambda gs: [len(g) for g in gs if len(g) > 1]
    count = lambda mods, acr, acrr: dict(
        modules=len(nontriv(mods)), giant=max(nontriv(mods), default=1), pairs=len(pairs([sorted(m) for m in mods])),
        acr=len(acr), acrr_pairs=len(acrr), acr_free=len(acr & free),
        acrr_pairs_free=sum(a in free and b in free for a, b in acrr))
    out = dict(seconds_iterated=t_it, seconds_tangent=t_tan, tangent_dim=tan.dimension,
               it=count(it.modules, acr_it, acrr_it), tan=count(tan.modules, acr_t, acrr_t),
               violations=dict(pairs=len(cp_it - cp_t), acr=len(acr_it - acr_t), acrr_pairs=len(acrr_it - acrr_t)))
    out["gap"] = any(out["it"][k] != out["tan"][k] for k in ("pairs", "acr", "acrr_pairs"))
    # gap relations, as in two_point_gap.py prepare
    mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
    _, find, parent = groups_of(acrr_it)
    groups_t, _, _ = groups_of(acrr_t)
    rel = [dict(kind="ACR", c={mets[i]: 1}) for i in sorted(acr_t - acr_it)]
    for g in groups_t:
        parts = {}
        for i in g:
            parts.setdefault(find(i) if i in parent else ("solo", i), []).append(i)
        reps = sorted(min(p) for p in parts.values())
        for k in reps[1:]:
            rel.append(dict(kind="ACRR", group=[mets[i] for i in sorted(g)], c={mets[reps[0]]: 1, mets[k]: -1}))
    out["relations"] = rel
    out["n_relations"] = len(rel)
    if rel:
        t0 = time.time()
        _, Vrows = exact_V(crn["Y"], nm, it.modules, it.acr)
        C = [{mets.index(k): Fraction(v) for k, v in r["c"].items()} for r in rel]
        dKV = kernel_dim([Vrows], nm)
        out.update(dim_kerV=dKV, g_rel=dKV - kernel_dim([Vrows, C], nm), seconds_g_rel=time.time() - t0)
    return out


if __name__ == "__main__":
    M = sys.argv[1]
    crn = json.load(open(records.network(M)))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    net = kimonu.crn_from_matrices(A, Y)                 # N = Y·A; metabolites by index
    mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
    free = set(range(mets.index("E_1") if "E_1" in mets else nm))
    out = dict(model=M, n_complexes=nc, n_metabolites=nm, n_free=len(free))
    only = sys.argv[2:]                           # one mode only, for the slow models
    for name, det in (("attached", False), ("detached", True)):
        if only and name not in only:
            continue
        out[name] = mode(crn, net, det, free)
        print(M, name, json.dumps({k: v for k, v in out[name].items() if k != "relations"}), flush=True)
    json.dump(out, open(records.write(records.UPTAKE_MODES, f"{M}{''.join('.' + o for o in only)}.json"), "w"), indent=1)
