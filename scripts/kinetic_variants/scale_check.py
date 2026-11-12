"""S8 check: the tangent route (and the iterated route) at the LP point rounded at scale 10^6 (the package default),
at scale 100, and unrounded (the LP point's float64 values as exact binary fractions, as before 2026-09-30).

    python kinetic_variants/scale_check.py <model> <out.json>     # from scripts/, PYTHONPATH with StrainDesign 1.19.2

Run on 2026-11-11 for iNJ661 and iHN637 (sadreyevlab, 2 cores, not exclusive) into
/scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/scale_check/<model>.json, read by figures/paper_numbers.py (S8).
iNJ661 was timed again on an exclusive node, one core, three runs (kinetic_variants/scale_check_exclusive.sbatch,
job 572320) into scale_check/exclusive/iNJ661_rep{1,2,3}.json; the paper's S8 times are their medians.
"""
import json, os, sys, time
from fractions import Fraction
from math import lcm
import numpy as np
from scipy.optimize import linprog
B = "/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/scripts"
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu import steady_state
from kimonu.linalg import nullspace, to_float
from run_kimonu_crn import triplets_to_sparse

M, dest = sys.argv[1], sys.argv[2]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
net = kimonu.crn_from_matrices(triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64), triplets_to_sparse(crn["Y"], nm, nc))


def unrounded(seed=0):
    """interior_flux's LP point without rounding: free coordinates b_j = r_i / K_ij with r_i the exact float64"""
    Nf = to_float(net.N); m, q = Nf.shape
    rng = np.random.default_rng(seed)
    lower = 1 + rng.random(q)
    lp = linprog(rng.random(q), A_eq=Nf, b_eq=np.zeros(m), bounds=[(lo, None) for lo in lower], method="highs")
    K = nullspace(net.N)
    free = {}
    for i, row in enumerate(K):
        if len(row) == 1:
            (j, v), = row
            free.setdefault(j, Fraction(float(lp.x[i])) / v)
    flux = [sum((v * free.get(j, 0) for j, v in row), Fraction(0)) for row in K]
    assert all(x > 0 for x in flux)
    d = lcm(*(x.denominator for x in flux))
    return [x * d for x in flux]


def summary(r, seconds, flux):
    return dict(seconds=seconds, n_modules=r.n_modules, giant=r.largest_module, acr=len(r.acr), acrr_pairs=len(r.acrr),
                acrr_groups=len(getattr(r, "acrr_groups", []) or []), max_bits=max(int(x).bit_length() for x in flux),
                modules=sorted(sorted(int(c) for c in s) for s in r.modules if len(s) > 1),
                acr_list=sorted(int(x) for x in r.acr), acrr_list=sorted(sorted(int(x) for x in p) for p in r.acrr))


out = dict(model=M, kimonu_lp_scale=steady_state._LP_SCALE)
fluxes = {}
for name, scale in (("scale_1e6", 10**6), ("scale_100", 100)):
    steady_state._LP_SCALE = scale
    fluxes[name] = kimonu.interior_flux(net, 0)
steady_state._LP_SCALE = 10**6
fluxes["unrounded"] = unrounded(0)
for route in ("tangent", "iterated"):
    for name, flux in fluxes.items():
        t0 = time.perf_counter()
        r = (kimonu.tangent_route if route == "tangent" else kimonu.iterated_route)(net, detach_uptakes=True, flux=flux)
        s = time.perf_counter() - t0
        out[f"{route}_{name}"] = summary(r, s, flux)
        print(M, route, name, f"{s:.1f} s", out[f"{route}_{name}"]["n_modules"], out[f"{route}_{name}"]["acr"],
              out[f"{route}_{name}"]["acrr_pairs"], "bits", out[f"{route}_{name}"]["max_bits"], flush=True)
for route in ("tangent", "iterated"):
    a = out[f"{route}_scale_1e6"]
    for name in ("scale_100", "unrounded"):
        b = out[f"{route}_{name}"]
        out[f"{route}_{name}"]["same_as_1e6"] = all(a[k] == b[k] for k in ("modules", "acr_list", "acrr_list"))
json.dump(out, open(dest, "w"), indent=1)
