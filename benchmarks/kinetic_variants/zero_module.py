"""Complexes with a fixed concentration (the complex counterpart of ACR), per uptake mode.

psi_c is fixed at every positive steady state of given rate constants exactly when complex c has a zero row of
Y^T B (B = free_directions of the iterated route's modules). Attached, these are the members of the zero complex's
module (psi_0 = 1); detached, the zero complex is left out of the readout and they form a module of their own.
Reports both counts per mode (they must agree attached), the fixed complexes by composition and the
single-metabolite ones among them, which are ACR metabolites.

    python kinetic_variants/zero_module.py <model>   -> $ZM_OUT/<model>.json
"""
import json, os, sys
import numpy as np
from scipy import sparse
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.linalg import matrix_product
from kimonu.readout import free_directions
from run_kimonu_crn import triplets_to_sparse

OUT = os.environ.get("ZM_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_accr")
M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
net = kimonu.crn_from_matrices(A, Y)                 # N = Y·A, the S of the crn.json
Yc = sparse.csc_matrix(Y)
mets = [m if isinstance(m, str) else m.get("id") for m in crn["metabolites"]]
comp = lambda c: {mets[i]: float(v) for i, v in zip(Yc.indices[Yc.indptr[c]:Yc.indptr[c + 1]], Yc.data[Yc.indptr[c]:Yc.indptr[c + 1]])}
zero = [c for c in range(nc) if Yc.indptr[c + 1] == Yc.indptr[c]]
out = dict(model=M, n_complexes=nc, zero=zero)
for mode, detach in (("attached", False), ("detached", True)):
    kin = kimonu.iterated_route(net, detach_uptakes=detach)
    mod = next((sorted(map(int, m)) for m in kin.modules if zero and zero[0] in m), [])
    members = [c for c in mod if c not in zero and c < nc]
    rows = free_directions(net, kin.modules)
    YtB = matrix_product(net.Yt, rows)
    fixed = [c for c in range(nc) if c not in zero and YtB[c] == ()]
    if not detach:
        assert fixed == members, (M, len(fixed), len(members))
    out[mode] = dict(n_zero_module=len(members), n_fixed=len(fixed), fixed=[comp(c) for c in fixed],
                     single=[next(iter(comp(c))) for c in fixed if len(comp(c)) == 1 and list(comp(c).values()) == [1.0]],
                     n_acr=len(kin.acr))
    print(M, mode, "zero module", len(members), "fixed", len(fixed), flush=True)
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"), indent=1)
