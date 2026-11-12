"""Does Figure 4 (the reference-equivalent mode) need the upstream sets? One round of ker(W), with and without one log closure, against
Kimonu's reference-equivalent mode (kimonu.reference_route: concordance, upstream sets merged where they
overlap, ACR/ACRR by the reference's pairwise rule).

Variants, each read out with the reference's pairwise rule (Shinar-Feinberg) and with the kernel rule (zero and equal rows
of ker((Y Delta)^T)):
    ref       the reference-equivalent mode (upstream sets)
    w1        the coupled sets of ker(W), one round, no closure
    w1c       w1 followed by one closure in log space (one readout of ker((Y Delta)^T)), no feedback
    full      the iterated route to its fixed point (kimonu.iterated_route)
Detached uptakes; the zero complex is left out of every module.

    python single_round.py <model>   -> $SR_OUT/<model>.json
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu import free_directions, readout
from kimonu.io import zero_complexes
from kimonu.iterated import coupled_columns
from kimonu.reference import shinar_feinberg
from run_kimonu_crn import triplets_to_sparse

OUT = os.environ.get("SR_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_single_round")
M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
net = kimonu.crn_from_matrices(A, Y)                 # N = Y·A, the S of the crn.json
mets = crn["metabolites"]; cut = mets.index("E_1") if "E_1" in mets else nm
skip = zero_complexes(net)
clean = lambda ms: sorted(sorted(int(c) for c in s if c not in skip) for s in ms if len(set(s) - skip) > 1)

t0 = time.time()
ref = clean(kimonu.reference_route(net).modules)
t_ref = time.time() - t0
t0 = time.time()
flux = kimonu.interior_flux(net, 0)
w1 = [{c for c in g if c < nc and c not in skip} for g in coupled_columns(net, flux, True)]
w1 = [s for s in w1 if s]
w1c = readout(net, free_directions(net, w1)).modules
t_w = time.time() - t0
full = kimonu.iterated_route(net, detach_uptakes=True, flux=flux).modules
variants = dict(ref=ref, w1=clean(w1), w1c=clean(w1c), full=clean(full))


def relations(mods):
    a, r = shinar_feinberg(net, mods)
    R = readout(net, free_directions(net, [set(m) for m in mods]))
    pairs = sorted(tuple(sorted((x, y))) for g in R.acrr_groups for i, x in enumerate(g) for y in g[i + 1:])
    free = lambda acr, acrr: dict(acr=sorted(int(i) for i in acr if i < cut),
                                  acrr=sorted([int(i), int(j)] for i, j in acrr if i < cut and j < cut))
    return dict(pairwise=free(a, r), kernel=free(R.acr, pairs))


out = dict(model=M, seconds_ref=t_ref, seconds_w1=t_w)
rel = {k: relations(v) for k, v in variants.items()}
for k, v in variants.items():
    out[k] = dict(n_kinetic=len(v), giant=max((len(m) for m in v), default=0),
                  modules_equal_ref=v == ref,
                  pairwise_equal_ref=rel[k]["pairwise"] == rel["ref"]["pairwise"],
                  kernel_equal_ref_pairwise=rel[k]["kernel"] == rel["ref"]["pairwise"],
                  acr_free_pairwise=len(rel[k]["pairwise"]["acr"]), acrr_free_pairwise=len(rel[k]["pairwise"]["acrr"]),
                  acr_free_kernel=len(rel[k]["kernel"]["acr"]), acrr_free_kernel=len(rel[k]["kernel"]["acrr"]))
    print(M, k, out[k], flush=True)
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"))
