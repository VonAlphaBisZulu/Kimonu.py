"""One timed run of the tangent route, for the runtime benchmarks.

    python tangent_cell.py <network> <out.json> [--attached] [--seed N]

<network> is a benchmark model (unblocked/<network>_fixed.crn.json), a small network of local_acr_symbolic.py, a
ladder step envz<k> or e_coli_core. Only kimonu.tangent_route is timed. Uptakes are detached unless --attached.
--seed sets the seed of the interior flux r*.
"""
import json, os, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from run_kimonu_crn import triplets_to_sparse

name, dest = sys.argv[1], sys.argv[2]
detached = "--attached" not in sys.argv
seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 0
path = f"{B}/unblocked/{name}_fixed.crn.json"
if os.path.exists(path):
    crn = json.load(open(path))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
else:
    from kinetic_variants import local_acr_symbolic as L
    if name in L.NETWORKS:
        A, Y, N, _ = L.matrices(*L.NETWORKS[name][:2])
    elif name.startswith("envz"):
        A, Y, N, _ = L.matrices(*L.envz_copies(int(name[4:])))
    else:
        A, Y, N, _ = L.ecoli_core(name)
net = kimonu.crn_from_matrices(A, Y)                     # N = Y·A
t0 = time.perf_counter()
r = kimonu.tangent_route(net, detach_uptakes=detached, seed=seed)
seconds = time.perf_counter() - t0
modules = sorted(sorted(int(c) for c in s) for s in r.modules if len(s) > 1)
acrr = sorted(sorted(int(x) for x in p) for p in r.acrr)
os.makedirs(os.path.dirname(dest), exist_ok=True)
json.dump(dict(network=name, detached=detached, seed=seed, species=int(Y.shape[0]), complexes=int(A.shape[0]), seconds=seconds,
               n_modules=r.n_modules, giant_module=r.largest_module, modules=modules,
               acr=sorted(int(x) for x in r.acr), acrr=acrr), open(dest, "w"))
print(f"{name} {'detached' if detached else 'attached'} {seconds:.2f} s, {r.n_modules} modules, "
      f"{len(r.acr)} ACR, {len(acrr)} ACRR pairs", flush=True)
