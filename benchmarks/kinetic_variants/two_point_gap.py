"""The relations only the tangent route reports, decided with the two-point criterion where possible.

    python kinetic_variants/two_point_gap.py <model> prepare               -> <out>/<model>_relations.json, _V.npz
    python kinetic_variants/two_point_gap.py <model> relation <i> [search s] [milp s]   -> <out>/<model>_rel<i>.json

prepare: iterated route and tangent route (attached uptakes), and one relation per missing link: each ACR
metabolite only the tangent route reports, and for each tangent ACRR group the links between the iterated
groups it joins (rep_1 - rep_k). Proving all links of a group proves the group.
relation: first `search` (a counterexample disproves the relation), then `certify` with the iterated route's
valid equalities (infeasible proves it). Anything else is undecided.
"""
import json, os, sys, time
import numpy as np, scipy.sparse as sp
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.two_point import substrate_composition, search, certify, iterated_constraints, realise

OUT = os.environ.get("TWO_POINT_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_twopoint")


def load(M):
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    Nf = triplets_to_sparse(crn["S"], nm, nr)
    return crn, A, Y, Nf


def prepare(M):
    import kimonu
    from kinetic_variants.tangent_one_shot import one_shot_exact
    crn, A, Y, Nf = load(M); mets = crn["metabolites"]; N = to_exact(Nf)
    V, classes, it = iterated_constraints(A, Y, N)
    ar = kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr)
    _, acr_t, groups_t, _ = one_shot_exact(A, Y, N, False, want_groups=True)
    # iterated ACRR groups (union of its pairs)
    parent = {}
    def find(a):
        while parent.setdefault(a, a) != a: a = parent[a]
        return a
    for a, b in ar.acrr: parent[find(a)] = find(b)
    rel = [dict(kind="ACR", c={mets[i]: 1}) for i in sorted(set(acr_t) - set(ar.acr))]
    for g in groups_t:
        parts = {}
        for i in g: parts.setdefault(find(i) if i in parent else ("solo", i), []).append(i)
        reps = [min(p) for p in parts.values()]
        if len(reps) > 1:
            for k in reps[1:]:
                rel.append(dict(kind="ACRR", group=[mets[i] for i in sorted(g)], c={mets[reps[0]]: 1, mets[k]: -1}))
    os.makedirs(OUT, exist_ok=True)
    sp.save_npz(f"{OUT}/{M}_V.npz", V.tocsr())
    json.dump(dict(model=M, classes=classes, relations=rel, n_V=V.shape[0],
                   n_classes=len(classes), n_reactions=Nf.shape[1], n_metabolites=Nf.shape[0]),
              open(f"{OUT}/{M}_relations.json", "w"), indent=1)
    print(f"{M}: {len(rel)} relations ({sum(r['kind'] == 'ACR' for r in rel)} ACR), V {V.shape}, {len(classes)} sign classes")


def relation(M, i, search_s=1200, milp_s=10800):
    crn, A, Y, Nf = load(M); mets = crn["metabolites"]
    info = json.load(open(f"{OUT}/{M}_relations.json")); r = info["relations"][i]
    V = sp.load_npz(f"{OUT}/{M}_V.npz")
    G = substrate_composition(A, Y)
    c = np.zeros(len(mets))
    for k, v in r["c"].items(): c[mets.index(k)] = v
    out = dict(model=M, index=i, relation=r)
    t0 = time.time()
    hit = search(Nf, G, c, V=V, starts=50, seed=i, time_limit=search_s)
    out["search_seconds"] = time.time() - t0
    if hit is not None:
        dxi, r1, r2, res = hit
        out.update(verdict="disproven", residual=res, c_dot_dxi=float(c @ dxi),
                   dxi={mets[k]: float(dxi[k]) for k in np.flatnonzero(np.abs(dxi) > 1e-9)})
    else:
        st, sec, dxi = certify(Nf, G, c, V=V, classes=info["classes"], time_limit=milp_s, threads=8)
        out.update(milp_status=st, milp_seconds=sec)
        if "infeasible" in st:
            out["verdict"] = "proven"
        else:
            out["verdict"] = "undecided"
            if dxi is not None:      # a sign pattern exists: try to realise it along its ray
                out["sign_pattern_realised"] = any(realise(Nf, G, t * dxi) is not None for t in (0.01, 0.1, 1, 10))
    print(json.dumps({k: v for k, v in out.items() if k != "dxi"}), flush=True)
    json.dump(out, open(f"{OUT}/{M}_rel{i}.json", "w"), indent=1)


if __name__ == "__main__":
    M, cmd = sys.argv[1], sys.argv[2]
    if cmd == "prepare":
        prepare(M)
    else:
        relation(M, int(sys.argv[3]), *(float(x) for x in sys.argv[4:6]))
