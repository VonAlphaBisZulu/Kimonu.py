"""S8 addendum: both routes at the LP point rounded at scale 10^4 and 10^5, against the record of scale_check.py at 10^6.

    python kinetic_variants/scale_check_1e4.py iHN637 <dir>/iHN637_1e4.json     # <dir> holds scale_check.py's iHN637.json
"""
import json, os, sys, time
import numpy as np
B = "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts"
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu import steady_state
from run_kimonu_crn import triplets_to_sparse
M, dest = sys.argv[1], sys.argv[2]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
net = kimonu.crn_from_matrices(triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64), triplets_to_sparse(crn["Y"], nm, nc))
rec = json.load(open(f"{os.path.dirname(dest)}/{M}.json"))
out = dict(model=M)
for scale in (10**4, 10**5):
    steady_state._LP_SCALE = scale
    flux = kimonu.interior_flux(net, 0)
    for route in ("tangent", "iterated"):
        t0 = time.perf_counter()
        r = (kimonu.tangent_route if route == "tangent" else kimonu.iterated_route)(net, detach_uptakes=True, flux=flux)
        s = time.perf_counter() - t0
        a = rec[f"{route}_scale_1e6"]
        same = (sorted(sorted(int(c) for c in x) for x in r.modules if len(x) > 1) == a["modules"]
                and sorted(int(x) for x in r.acr) == a["acr_list"] and sorted(sorted(int(x) for x in p) for p in r.acrr) == a["acrr_list"])
        out[f"{route}_scale_{scale}"] = dict(seconds=s, same_as_1e6=same, max_bits=max(int(x).bit_length() for x in flux))
        print(M, route, scale, f"{s:.1f} s", "same", same, flush=True)
json.dump(out, open(dest, "w"), indent=1)
