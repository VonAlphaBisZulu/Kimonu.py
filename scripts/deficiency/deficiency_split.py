"""Deficiency of the 33 benchmark models and its relation to Kimonu's concordance partition.

    python deficiency/deficiency_split.py [model ...]   # from scripts/ -> networks/deficiency/deficiency_33.json of KIMONU_OUT
"""
import json, sys, glob, os, time
from fractions import Fraction
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
import straindesign as sd
from kimonu.linalg import exact, matrix_product, entries
from kimonu.readout import row_readout
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import records                                         # networks and figdata.json from KIMONU_RECORDS, output to KIMONU_OUT

fig = json.load(open(records.read(records.SOURCE_DATA, 'figdata.json')))['kimonu_ref']
check_rank = set(sys.argv[1:])
rows = []
for f in sorted(glob.glob(records.read(records.UNBLOCKED, '*_fixed.crn.json'))):
    m = os.path.basename(f)[:-len('_fixed.crn.json')]
    d = json.load(open(f)); n, q = d['n_complexes'], d['n_reactions']
    A = sparse.csr_matrix(([v for _, _, v in d['A']], ([i - 1 for i, _, _ in d['A']], [j - 1 for _, j, _ in d['A']])), shape=(n, q))
    # zero complex: rows of A with no Y entries; a complex row may be empty in Y (zero complex)
    ycols = {j - 1 for _, j, _ in d['Y']}
    nonzero_cx = sum(1 for c in range(n) if c in ycols)
    # linkage classes: components of the undirected complex graph, edges = reactions
    G = (abs(A) @ abs(A).T).tocsr()
    ell, _ = connected_components(G, directed=False)
    used = np.asarray(abs(A).sum(axis=1)).ravel() > 0      # complexes that occur in some reaction
    ell -= int((~used).sum())                               # isolated unused rows are not linkage classes
    n_used = int(used.sum())
    N = exact([(i - 1, j - 1, Fraction(str(v))) for i, j, v in d['S']], (d['n_metabolites'], q))
    t = time.time(); K = sd.sparse_nullspace(N); dK = K.shape[1]
    rankN = q - dK; delta = n_used - ell - rankN
    Mrows = matrix_product(A.astype(np.int64), K)
    zero, classes, _ = row_readout(Mrows)
    nonsingle = sum(1 for c in classes if len(c) > 1)
    rankM = None
    if m in check_rank:
        trip = [(i, j, v) for i, r in enumerate(Mrows) for j, v in r]
        rankM = dK - sd.sparse_nullspace(exact(trip, (n, dK))).shape[1]
    fd = fig.get(m, {})
    rows.append(dict(model=m, n=n_used, zero_cx=n_used - sum(1 for c in range(n) if c in ycols and used[c]), ell=ell, rankN=rankN, delta=delta,
                     dimK=dK, balanced=len(zero), conc_nonsingle=nonsingle, directions=len(classes), rankM=rankM,
                     fig_bal=fd.get('balanced'), fig_conc=fd.get('conc'), sec=round(time.time() - t, 1)))
    print(json.dumps(rows[-1]), flush=True)
json.dump(rows, open(records.write(records.DEFICIENCY, 'deficiency_33.json'), 'w'), indent=1)
