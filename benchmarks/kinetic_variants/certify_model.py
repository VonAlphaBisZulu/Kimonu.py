"""Sign certificate for one metabolite of a full model: infeasible means fixed at every positive steady state.

    python kinetic_variants/certify_model.py <model> <metabolite> [time limit s]
"""
import json, os, sys, numpy as np, scipy.sparse as sp
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, B)
from run_kimonu_crn import triplets_to_sparse
from kinetic_variants.branch_certificate import certify
M, t = sys.argv[1], sys.argv[2]; tl = int(sys.argv[3]) if len(sys.argv) > 3 else 3600
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc); N = triplets_to_sparse(crn["S"], nm, nr)
G = ((A < 0).astype(float).T @ Y.T).tocsr()
c = np.zeros(nm); c[crn["metabolites"].index(t)] = 1
print(M, t, certify(N, G, c, time_limit=tl, threads=4), flush=True)
