"""ψ-nullspace coupling followed by the log-space closure (Proposition S3-4 with known-ACR augmentation),
against Kimonu's current default (upstream algorithm + S3-4 + S4-1 + S4-5).

Closure: merge coupled sets whose signatures Lᵀy agree, L = ker(YΔᵀ) with YΔ over all current sets
plus unit vectors of metabolites already known to have ACR; iterate until nothing changes. The
deficiency-one lemma is not applied here; where Kimonu's default uses it, it is reported.

    python kinetic_variants/psi_closure.py <model> [...]   -> kinetic_variants/closure_<model>.json
"""
import json, sys, os, time
import numpy as np, scipy.sparse as sp
from fractions import Fraction
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
exec(open(os.path.join(B, "kinetic_variants/psi_nullspace.py")).read().split("DETACH = ")[0])
from kimonu.acr import coupling_signatures, acr_acrr

def psi_blocks(crn, A, Y, Sx, nc, nr, nm, draw=1):
    Kx = sd.sparse_nullspace(Sx)
    if not isinstance(Kx, sd.ExactCOO):
        Kc = Kx.tocoo(); Kx = type("K", (), dict(rows=Kc.row, cols=Kc.col, data=Kc.data, denom=1))()
    sub = {int(j): int(i) for i, j, v in zip(*sp.find(A)) if v < 0}
    Nent = {}
    for i, j, v in crn["S"]: Nent.setdefault(int(j) - 1, []).append((int(i) - 1, Fraction(str(v))))
    zero = [c for c in range(nc) if Y[:, c].nnz == 0]
    v = exact_positive_flux(crn, Sx, Kx, draw); W = {}; extra = {}
    for r_, s in sub.items():
        col = extra.setdefault(r_, nc + len(extra)) if s in zero else s
        for i, a in Nent.get(r_, []): W[(i, col)] = W.get((i, col), Fraction(0)) + v[r_] * a
    ncol = nc + len(extra)
    psi = kimonu.concordance_from_matrices(sp.identity(ncol, format="csr", dtype=np.int64),
                                           sd.RationalMatrix.from_fractions(((i, j, x) for (i, j), x in W.items() if x != 0), (nm, ncol)))
    return [set(c for c in m if c < nc and c not in zero) for m in psi.modules]

def closure(Y, sets, reactant=None):
    """Singletons take part when they are the substrate of some reaction: S3-4 can attach a lone
    complex to a set. Terminal complexes never enter a rate law and stay out."""
    sets = [s for s in sets if len(s) > 1 or (len(s) == 1 and reactant is not None and next(iter(s)) in reactant)]; known = set(); rounds = 0
    while True:
        rounds += 1
        sig = coupling_signatures(Y, sets, known)
        if sig is None: break
        g = {}
        for s, k in zip(sets, sig): g.setdefault(k, set()).update(s)
        new = list(g.values())
        known |= set(acr_acrr(Y, new, known_acr=known).acr)
        if len(new) == len(sets): break
        sets = new
    return sets, sorted(known), rounds

for M in sys.argv[1:]:
    try:
        crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
        A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc); Sx = to_exact(triplets_to_sparse(crn["S"], nm, nr))
        mets = crn["metabolites"]; cut = mets.index("E_1") if "E_1" in mets else len(mets)
        t0 = time.time(); blocks_ = psi_blocks(crn, A, Y, Sx, nc, nr, nm)
        reactant = {int(i) for i, j, v in zip(*sp.find(A)) if v < 0}
        mods, known, rounds = closure(Y, blocks_, reactant); mods = [m for m in mods if len(m) > 1]; t_new = time.time() - t0
        a_new = acr_acrr(Y, mods, mets, known_acr=[mets.index(x) for x in known] if known and isinstance(known[0], str) else known)
        t0 = time.time(); res = kimonu.concordance_from_matrices(A, Sx)
        kd = upstream.kinetic_modules(A, res.modules, res.balanced, Y=Y, nullity=res.nullity); t_old = time.time() - t0
        a_old = acr_acrr(Y, kd.modules, mets, known_acr=kd.known_acr)
        P = {frozenset(m) for m in mods}; Q = {frozenset(m) for m in kd.modules if len(m) > 1}
        free = lambda pairs: sum(1 for p in pairs if all(mets.index(x) < cut for x in p))
        out = dict(model=M, new=dict(n_modules=len(P), giant=max((len(m) for m in P), default=0), acr=len(a_new.acr), acrr=len(a_new.acrr),
                                     acr_free=sum(1 for x in a_new.acr if mets.index(x) < cut), acrr_free=free(a_new.acrr), rounds=rounds, seconds=t_new),
                   default=dict(n_modules=len(Q), giant=kd.giant_module, acr=len(a_old.acr), acrr=len(a_old.acrr),
                                acr_free=sum(1 for x in a_old.acr if mets.index(x) < cut), acrr_free=free(a_old.acrr), seconds=t_old),
                   relation=relation(Q, P), acr_superset=set(a_old.acr) <= set(a_new.acr),
                   acrr_superset={frozenset(p) for p in a_old.acrr} <= {frozenset(p) for p in a_new.acrr})
        json.dump(out, open(f"{B}/kinetic_variants/closure_{M}.json", "w"), indent=1)
        n, d = out["new"], out["default"]
        print(f"{M:26s} new: mods {n['n_modules']:5d} giant {n['giant']:4d} ACR {n['acr']:3d} ACRR {n['acrr']:6d} (free {n['acrr_free']:4d}) | default: mods {d['n_modules']:5d} giant {d['giant']:4d} ACR {d['acr']:3d} ACRR {d['acrr']:6d} (free {d['acrr_free']:4d}) | default {out['relation']:9s} new | ACR⊆ {out['acr_superset']} ACRR⊆ {out['acrr_superset']} | {n['seconds']:.0f}s vs {d['seconds']:.0f}s", flush=True)
    except Exception as e:
        print(f"{M:26s} FAILED {type(e).__name__}: {e}", flush=True)
