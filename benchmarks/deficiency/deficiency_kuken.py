"""Deficiency and exact concordance on Küken et al. 2022's own preprocessed networks (scenario i).

Their result files hold the network after their preprocessing (split, blocked removed at 1e-6): A, Y, S.
We recompute balanced complexes and concordance exactly on that same network, so any difference is
the method, not the input.
"""
import json, glob, os, sys, time
from fractions import Fraction
import numpy as np
import scipy.io as sio
from scipy import sparse
from scipy.sparse.csgraph import connected_components
import straindesign as sd
from kimonu.linalg import exact, matrix_product
from kimonu.readout import row_readout

D = '/scratch/sadreyev_lab/pschneider/concordant_complexes/Results_concordant_original_irreversibility_considered'
rows = []
for f in sorted(glob.glob(f'{D}/*.mat')):
    name = os.path.basename(f)[:-4]
    d = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
    m = d['Results_balanced'].MODEL_r
    A = sparse.csr_matrix(np.asarray(m.A, dtype=float)); Y = sparse.csr_matrix(np.asarray(m.Y, dtype=float))
    Sm = sparse.csr_matrix(np.asarray(m.S, dtype=float))
    # remove reactions that cannot carry flux under the stored bounds (all lb = 0 here), then unused complexes
    lb = np.asarray(m.lb, float); ub = np.asarray(m.ub, float); q0 = A.shape[1]
    assert (lb >= 0).all(); n_lbpos = int((lb > 0).sum())       # Küken: lb set to 0 after the split; report if not
    from scipy.optimize import linprog
    keep = []
    for j in range(q0):
        c = np.zeros(q0); c[j] = -1
        r = linprog(c, A_eq=Sm, b_eq=np.zeros(Sm.shape[0]), bounds=list(zip(lb, ub)), method='highs')
        if r.status == 0 and -r.fun > 1e-9: keep.append(j)
    n_blocked = q0 - len(keep)
    A, Sm = A[:, keep], Sm[:, keep]
    usedc = np.asarray(abs(A).sum(axis=1)).ravel() > 0
    A, Y = A[usedc], Y[:, usedc]
    n, q = A.shape
    colsum = np.asarray(A.sum(axis=0)).ravel()
    one_sided = int((np.asarray(abs(A).sum(axis=0)).ravel() == 1).sum())       # exchanges without a zero complex
    ya_ok = abs(Y @ A - Sm).max() < 1e-9
    # exact N = Y·A from the decimal coefficients of Y
    YA = (Y @ A).tocoo()
    N = exact([(int(i), int(j), Fraction(str(round(v, 12)))) for i, j, v in zip(YA.row, YA.col, YA.data) if abs(v) > 1e-12], YA.shape)
    t = time.time(); K = sd.sparse_nullspace(N); dK = K.shape[1]; rankN = q - dK
    G = (abs(A) @ abs(A).T).tocsr(); ell, _ = connected_components(G, directed=False)
    Ai = sparse.csr_matrix(np.rint(A.toarray()).astype(np.int64))
    Mrows = matrix_product(Ai, K)
    zero, classes, _ = row_readout(Mrows)
    nonsingle = [c for c in classes if len(c) > 1]
    trip = [(i, j, v) for i, r in enumerate(Mrows) for j, v in r]
    rankM = dK - sd.sparse_nullspace(exact(trip, (n, dK))).shape[1] if trip else 0
    kB = len(np.atleast_1d(d['B'])); kcls = len(np.atleast_1d(d['class'])); kcwb = len(np.atleast_1d(d['class_with_balanced']))
    ksizes = [len(np.atleast_1d(c)) for c in np.atleast_1d(d['class'])]
    r = dict(model=name, n_stored=int(len(usedc)), blocked_removed=n_blocked, lb_positive=n_lbpos, n=n, q=q, one_sided_exchanges=one_sided, YA_eq_S=bool(ya_ok), ell=int(ell), rankN=rankN,
             delta=n - int(ell) - rankN, rankM=rankM, dimK=dK, balanced=len(zero), modules_nonsingle=len(nonsingle),
             directions=len(classes), singletons=len(classes) - len(nonsingle),
             kuken_balanced=kB, kuken_modules=kcls, kuken_m=kcwb, kuken_in_modules=sum(ksizes),
             kuken_index=1 - kcwb / n, index_ours_same_def=1 - (len(nonsingle) + (1 if zero else 0)) / n,
             index_all_classes=1 - (len(classes) + (1 if zero else 0)) / n, floor_index=1 - (n - int(ell) - rankN) / n,
             sec=round(time.time() - t, 1))
    rows.append(r); print(json.dumps(r), flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'deficiency_kuken12_unblocked.json'), 'w'), indent=1)
