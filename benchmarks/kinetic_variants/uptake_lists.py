"""Which ACR metabolites detaching removes, and where they go: iterated route, both uptake modes, per model.

Records the iterated route's ACR metabolites and ACRR groups (names) for attached and detached uptakes, and
for every metabolite that is ACR attached but not detached whether it lies in a detached ACRR group (its
level then follows the supply while its ratio to other metabolites stays fixed).

    python kinetic_variants/uptake_lists.py <model> <out_dir>   -> <out_dir>/<model>.json
"""
import json, os, sys
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from run_kimonu_crn import triplets_to_sparse, to_exact

M, out_dir = sys.argv[1], sys.argv[2]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
free = set(range(mets.index("E_1") if "E_1" in mets else nm))
out = dict(model=M)
for name, det in (("attached", False), ("detached", True)):
    it = kimonu.iterated_route(A, Y, N, detach_uptakes=det)
    parent = {}
    def find(a):
        while parent.setdefault(a, a) != a:
            a = parent[a]
        return a
    for a, b in it.acrr:
        parent[find(a)] = find(b)
    groups = {}
    for a in parent:
        groups.setdefault(find(a), []).append(a)
    out[name] = dict(acr=sorted(int(i) for i in it.acr), acrr_groups=[sorted(int(i) for i in g) for g in groups.values()])
att, det = set(out["attached"]["acr"]), set(out["detached"]["acr"])
in_group = {i for g in out["detached"]["acrr_groups"] for i in g}
lost = sorted(att - det)
out["summary"] = dict(acr_attached=len(att), acr_detached=len(det), acr_attached_free=len(att & free),
                      acr_detached_free=len(det & free), lost=len(lost), lost_free=len(set(lost) & free),
                      lost_in_detached_acrr=len(set(lost) & in_group), lost_free_in_detached_acrr=len(set(lost) & in_group & free),
                      gained=len(det - att), lost_names=[mets[i] for i in lost if i in free][:40])
os.makedirs(out_dir, exist_ok=True)
json.dump(out, open(f"{out_dir}/{M}.json", "w"))
print(M, json.dumps(out["summary"]), flush=True)
