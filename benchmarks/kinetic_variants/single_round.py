"""Does Figure 4 (the reference-equivalent mode) need the upstream sets? One round of ker(W), with and without one log closure, against
Kimonu's reference-equivalent mode (run_kimonu_refmode.py: concordance, upstream sets merged where they
overlap, ACR/ACRR by the reference's pairwise rule).

Variants, each read out with the reference's pairwise rule (s210) and with the kernel rule (zero and equal rows
of ker((Y Delta)^T)):
    ref       the reference-equivalent mode (upstream sets)
    w1        the coupled sets of ker(W), one round, no closure
    w1c       w1 followed by one closure in log space (one readout of ker((Y Delta)^T)), no feedback
    full      the iterated route to its fixed point (kimonu.kinetic_modules)
Detached uptakes; the zero complex is left out of every module.

    python single_round.py <model>   -> $SR_OUT/<model>.json
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from compat import upstream
from kimonu.acr import free_directions, readout
from kimonu.io import zero_complexes
from kimonu.kinetic import _closure, _psi_coupled
from kimonu.linalg import positive_flux
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.s210_criterion import s210

OUT = os.environ.get("SR_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_single_round")
M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
mets = crn["metabolites"]; cut = mets.index("E_1") if "E_1" in mets else nm
skip = zero_complexes(Y)
clean = lambda ms: sorted(sorted(int(c) for c in s if c not in skip) for s in ms if len(set(s) - skip) > 1)

t0 = time.time()
res = kimonu.concordance_from_matrices(A, N)
ref = clean(upstream.kinetic_modules(A, [g for g in res.modules if len(g) > 1], res.balanced, Y=None).modules)
t_ref = time.time() - t0
t0 = time.time()
flux = positive_flux(N, 0)
w1 = [s for s in _psi_coupled(A, Y, N, flux, True)[0] if s]
w1c, _ = _closure(Y, w1, skip)
t_w = time.time() - t0
full = kimonu.kinetic_modules(A, Y, N, detach_uptakes=True).modules
variants = dict(ref=ref, w1=clean(w1), w1c=clean(w1c), full=clean(full))


def relations(mods):
    a, r = s210(Y, mods)
    R = readout(Y, free_directions(Y, [set(m) for m in mods]), skip)
    pairs = sorted(tuple(sorted((x, y))) for g in R.acrr for i, x in enumerate(g) for y in g[i + 1:])
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
