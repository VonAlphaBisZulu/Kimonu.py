"""Write F2C2-shaped coupling matrices from the exact nullspace, for the reference ACR stage.

Get_excel_tables_kinetic_modules.R reads q$coupling.matrix and keeps only code 1 (fully coupled),
zeroing 2/3/4. Full coupling is proportionality of rows of K = ker(N), so code 1 is all it needs
and all we write; everything else is 0.
"""
import json, sys, os
from fractions import Fraction
import numpy as np, scipy.sparse as sp, scipy.io as sio
import straindesign as sd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import records                                         # paths of the records (KIMONU_RECORDS, KIMONU_OUT)
from run_kimonu_crn import triplets_to_sparse, to_exact
OUT=records.write(records.COUPLING_MATRICES, folder=True)    # networks/reference_coupling_matrices of KIMONU_OUT
for M in sys.argv[1:]:
    crn=json.load(open(records.network(M))); nr,nm=crn["n_reactions"],crn["n_metabolites"]
    K=sd.sparse_nullspace(to_exact(triplets_to_sparse(crn["S"],nm,nr)))
    rows={}
    if isinstance(K, sd.ExactCOO):
        for r,c,v in zip(K.rows,K.cols,K.data): rows.setdefault(int(r),[]).append((int(c),int(v)))
    else:
        Kc=K.tocsr()
        for i in range(Kc.shape[0]):
            rr=Kc.getrow(i); rows[i]=[(int(c),int(v)) for c,v in zip(rr.indices,rr.data)]
    groups={}
    for i in range(nr):
        nz=sorted(rows.get(i,[]))
        if not nz: continue
        c0,v0=nz[0]; key=(tuple(c for c,_ in nz), tuple(Fraction(v,v0) for _,v in nz))
        groups.setdefault(key,[]).append(i)
    I=[];J=[]
    for g in groups.values():
        if len(g)<2: continue
        for a in g:
            for b in g:
                if a!=b: I.append(a); J.append(b)
    cm=sp.csc_matrix((np.ones(len(I)),(I,J)),shape=(nr,nr))
    sio.savemat(f"{OUT}/{M}_stoichiometric_coupling_fixed.mat", {"coupling_matrix":cm}, do_compression=True)
    print(f"{M:26s} reactions {nr:6d}  fully-coupled classes {sum(1 for g in groups.values() if len(g)>1):4d}  edges {len(I)//2:6d}  -> {OUT}/{M}_stoichiometric_coupling_fixed.mat")
