"""ACR/ACRR by the reference's combinatorial rule (Corollary S2-10 as implemented in
Get_excel_tables_kinetic_modules.R) on Kimonu's own kinetic modules, for comparability.

ACR: two complexes of one module whose compositions differ in exactly one species.
ACRR: two complexes differing in exactly two species, each complex having exactly one further
species in common with the other (the A+X vs A+Y pattern; no sign or coefficient test).

    python kinetic_variants/s210_criterion.py <model> [...]   -> kinetic_variants/s210_<model>.json
"""
import json, sys, os
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import kimonu
from compat import upstream
from run_kimonu_crn import triplets_to_sparse, to_exact

def s210(Y, modules):
    Y = Y.tocsc(); acr, acrr = set(), set()
    for mod in modules:
        mem = sorted(mod)
        if len(mem) < 2: continue
        cols = {c: dict(zip(Y[:, c].nonzero()[0].tolist(), Y[:, c].data.tolist())) for c in mem}
        for a in range(len(mem)):
            ca = cols[mem[a]]
            for b in range(a + 1, len(mem)):
                cb = cols[mem[b]]
                keys = set(ca) | set(cb)
                diff = [k for k in keys if ca.get(k, 0) != cb.get(k, 0)]
                if len(diff) == 1:
                    acr.add(diff[0])
                elif len(diff) == 2 and len(set(ca) - set(diff)) == 1 and len(set(cb) - set(diff)) == 1:
                    acrr.add(tuple(sorted(diff)))
    return acr, acrr

if __name__ == "__main__":
    for M in sys.argv[1:]:
        crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
        nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
        A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
        S = to_exact(triplets_to_sparse(crn["S"], nm, nr))
        res = kimonu.concordance_from_matrices(A, S)
        kin = upstream.kinetic_modules(A, res.modules, res.balanced, Y=Y, nullity=res.nullity)
        acr, acrr = s210(Y, kin.modules)
        mets = crn["metabolites"]; cut = mets.index("E_1") if "E_1" in mets else len(mets)
        out = dict(model=M, acr=sorted(mets[i] for i in acr), acrr=sorted([mets[i], mets[j]] for i, j in acrr),
                   acr_free=sum(1 for i in acr if i < cut), acrr_free=sum(1 for i, j in acrr if i < cut and j < cut),
                   n_kinetic=sum(1 for m in kin.modules if len(m) > 1), giant=kin.giant_module)
        json.dump(out, open(f"{B}/kinetic_variants/s210_{M}.json", "w"), indent=1)
        print(f"{M:24s} ACR {len(acr):3d} (free {out['acr_free']})  ACRR {len(acrr):5d} (free {out['acrr_free']})  kinetic modules {out['n_kinetic']}", flush=True)
