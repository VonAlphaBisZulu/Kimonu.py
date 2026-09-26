"""Classify the complexes the ψ-nullspace (detached) couples beyond the upstream algorithm (overlap merge).

For every ψ block that is not an upstream module, each member that is not in the upstream module it
joins is labelled: uptake-fed (produced by a reaction out of the zero complex), terminal (in a terminal
strong component of the complex graph), or other. Writes kinetic_variants/psi_extras.json.
"""
import json, sys, os
import numpy as np, scipy.sparse as sp
from fractions import Fraction
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
src = open(os.path.join(B, "kinetic_variants/psi_nullspace.py")).read().split("DETACH = ")[0]; exec(src)
from compat.upstream import _Topology, _sccs
out = {}
for M in sys.argv[1:]:
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc).tocsc(); Sx = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    res = kimonu.concordance_from_matrices(A, Sx); ko = upstream.kinetic_modules(A, res.modules, res.balanced, Y=None)
    Kx = sd.sparse_nullspace(Sx)
    if not isinstance(Kx, sd.ExactCOO):
        Kc = Kx.tocoo(); Kx = type("K", (), dict(rows=Kc.row, cols=Kc.col, data=Kc.data, denom=1))()
    sub = {int(j): int(i) for i, j, v in zip(*sp.find(A)) if v < 0}
    Nent = {}
    for i, j, v in crn["S"]: Nent.setdefault(int(j) - 1, []).append((int(i) - 1, Fraction(str(v))))
    zero = [c for c in range(nc) if Y[:, c].nnz == 0]
    v = exact_positive_flux(crn, Sx, Kx, 1); W = {}; extra = {}
    for r_, s in sub.items():
        col = extra.setdefault(r_, nc + len(extra)) if s in zero else s
        for i, a in Nent.get(r_, []): W[(i, col)] = W.get((i, col), Fraction(0)) + v[r_] * a
    ncol = nc + len(extra)
    psi = kimonu.concordance_from_matrices(sp.identity(ncol, format="csr", dtype=np.int64),
                                           sd.RationalMatrix.from_fractions(((i, j, x) for (i, j), x in W.items() if x != 0), (nm, ncol)))
    mods = [set(c for c in m if c < nc and c not in zero) for m in psi.modules]
    t = _Topology(A); term = set()
    for comp in _sccs(t.succ, set(range(nc))):
        if not any(t.succ[x] - comp for x in comp): term |= comp
    fed = {int(i) for i, j, val in zip(*sp.find(A)) if val > 0 and sub.get(int(j)) in zero}
    of = {c: frozenset(m) for m in ko.modules for c in m}
    counts = {"uptake-fed": 0, "terminal": 0, "other": 0}; examples = {}
    for m in mods:
        if len(m) < 2 or frozenset(m) in {frozenset(x) for x in ko.modules}: continue
        core = max((of[c] for c in m), key=len)
        for c in m:
            if c in core: continue
            kind = "uptake-fed" if c in fed else "terminal" if c in term else "other"
            counts[kind] += 1
            nm_ = "+".join(crn["metabolites"][i] for i in Y[:, c].nonzero()[0]) or "∅"
            examples.setdefault(kind, []).append(nm_)
    out[M] = dict(counts=counts, examples={k: v[:4] for k, v in examples.items()})
    print(f"{M:26s} {counts}  other e.g. {examples.get('other', [])[:3]}", flush=True)
json.dump(out, open(f"{B}/kinetic_variants/psi_extras.json", "w"), indent=1)
