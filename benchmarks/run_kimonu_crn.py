#!/usr/bin/env python
"""Run Kimonu on a pre-decomposed network (crn.json), rather than on a cobra model.

Used for the reference implementation's own elementary-step-split models, which already carry A
and Y. Analysing them directly means the comparison against the published numbers is made on the
paper's own input, with no preprocessing of ours in between.

Usage:
    python benchmarks/run_kimonu_crn.py <crn.json> <out.json> [--repeats N]
"""
import argparse
import json
import time
from fractions import Fraction
from math import gcd

import numpy as np
from scipy import sparse

import kimonu
from kimonu.io import exact_matrix
from kimonu.kinetic import complex_graph, _sccs


def triplets_to_sparse(t, nrows, ncols, dtype=float):
    if not t:
        return sparse.csr_matrix((nrows, ncols), dtype=dtype)
    arr = np.asarray(t, dtype=float)
    return sparse.csr_matrix((arr[:, 2].astype(dtype), (arr[:, 0].astype(int) - 1,
                                                        arr[:, 1].astype(int) - 1)),
                             shape=(nrows, ncols), dtype=dtype)


def to_exact(M):
    """Exact rational form of ``M`` for the kernel computation."""
    M = M.tocoo()
    return exact_matrix(zip(M.row, M.col, M.data), M.shape)



def default_solver():
    """Prefer a free solver that is exact here.

    SCIP reproduces Gurobi's blocked set exactly on iSynCJ816, iJB785, iHN637 and iLJ478 at the
    same speed; GLPK reports one reaction blocked on iLJ478 that both others can drive nonzero,
    and a false positive is the damaging direction -- the reaction is deleted, the flux space
    shrinks, and complexes are then called balanced that are not.
    """
    import straindesign as sd
    for cand in ("scip", "gurobi", "cplex", "glpk"):
        if cand in sd.avail_solvers:
            return cand
    raise RuntimeError("no LP solver available")


def unblocked_reactions(crn, solver):
    """Indices of reactions that can carry flux, by FVA on the network the crn describes."""
    import cobra
    import straindesign as sd

    # cobra attaches its configured solver to every new Model, which is Gurobi by default here and
    # would take a licence seat even though the FVA below runs on `solver`.
    previous = cobra.Configuration().solver
    cobra.Configuration().solver = solver if solver in ("glpk", "gurobi", "cplex") else "glpk"
    model = cobra.Model("split")
    mets = [cobra.Metabolite(m) for m in crn["metabolites"]]
    model.add_metabolites(mets)
    terms = {}
    for i, j, v in crn["S"]:
        terms.setdefault(int(j) - 1, {})[int(i) - 1] = float(v)
    rxns = [cobra.Reaction(f"r{j}") for j in range(crn["n_reactions"])]
    for j, r in enumerate(rxns):
        r.lower_bound = float(crn["lb"][j])
        r.upper_bound = float(crn["ub"][j])
    model.add_reactions(rxns)
    for j, t in terms.items():
        rxns[j].add_metabolites({mets[i]: c for i, c in t.items()})
    try:
        ranges = sd.fva(model, solver=solver)
    finally:
        cobra.Configuration().solver = previous
    return [j for j in range(crn["n_reactions"])
            if not (abs(ranges.iat[j, 0]) <= 1e-9 and abs(ranges.iat[j, 1]) <= 1e-9)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("crn")
    ap.add_argument("out")
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--assume-unblocked", action="store_true",
                    help="the network already has its blocked reactions removed, so skip the FVA. "
                         "Used when every tool is handed the same pre-unblocked network, so the "
                         "shared preprocessing is costed once rather than inside one arm.")
    ap.add_argument("--solver", default=None,
                    help="LP solver for the blocked-reaction FVA; default picks the best free one.")
    args = ap.parse_args()

    with open(args.crn) as fh:
        crn = json.load(fh)

    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64)
    Y = triplets_to_sparse(crn["Y"], nm, nc)
    Sf = triplets_to_sparse(crn["S"], nm, nr)

    # "Balanced" means zero activity over the feasible set, which the kernel test reproduces only
    # when no reaction is blocked: span{v : Nv = 0, v >= 0} equals ker(N) minus the blocked
    # coordinates. Dropping them first is what makes the two definitions agree.
    n_blocked, t_blocked = 0, 0.0
    if not args.assume_unblocked:
        t0 = time.perf_counter()
        keep = unblocked_reactions(crn, args.solver or default_solver())
        t_blocked = time.perf_counter() - t0
        n_blocked = nr - len(keep)
        if n_blocked:
            A = A[:, keep]
            Sf = Sf[:, keep]
            nr = len(keep)
    S = to_exact(Sf)

    times, t_conc, t_kin, t_acr = [], [], [], []
    for _ in range(args.repeats):
        t0 = time.perf_counter()
        result = kimonu.concordance_from_matrices(A, S)
        t1 = time.perf_counter()
        kin = kimonu.kinetic_modules(A, result.modules, result.balanced, Y=Y,
                                     nullity=result.nullity)
        t2 = time.perf_counter()
        acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=crn["metabolites"],
                              known_acr=kin.known_acr)
        t3 = time.perf_counter()
        t_conc.append(t1 - t0)
        t_kin.append(t2 - t1)
        t_acr.append(t3 - t2)
        times.append(t3 - t0)

    module_of = {}
    for module_id, group in enumerate(result.modules, start=1):
        for i in group:
            module_of[i] = module_id if len(group) > 1 else -1
    for i in result.balanced:
        module_of[i] = 0

    Ycsc = Y.tocsc()
    # The empty complex has no composition; it exists so boundary reactions have a neighbour in
    # the complex graph, where terminality is judged. It is never a result, so it is left out of
    # what is reported -- as COCOA does -- and the complex counts then line up between the two.
    zero_complexes = {j for j in range(nc) if Ycsc[:, j].nnz == 0}
    complexes = []
    for j in range(nc):
        if j in zero_complexes:
            continue
        start, end = Ycsc.indptr[j], Ycsc.indptr[j + 1]
        comp = {crn["metabolites"][int(r)]: float(v)
                for r, v in zip(Ycsc.indices[start:end], Ycsc.data[start:end])}
        complexes.append({"id": f"complex_{j + 1}", "module": module_of[j], "composition": comp})

    succ, pred = complex_graph(A)
    sccs = _sccs(succ, set(range(nc)))
    terminal = sum(1 for c in sccs if not any(w not in c for v in c for w in succ[v]))

    # The supplement panels count free metabolites only: the reference's elementary-step split
    # appends enzyme species after E_1, and those are not part of the original network.
    mets = crn["metabolites"]
    free_upto = mets.index("E_1") if "E_1" in mets else len(mets)
    free = set(mets[:free_upto])

    for j, c in enumerate(complexes):
        c["kinetic_module"] = int(kin.module_of[j]) if kin.module_of[j] is not None else -1

    out = {"tool": "Kimonu.py", "seconds": min(times), "seconds_cold": times[0],
           "seconds_concordance": min(t_conc), "seconds_kinetic": min(t_kin),
           "seconds_acr": min(t_acr), "repeats": args.repeats,
           "n_reactions": nr, "n_metabolites": nm, "n_complexes": nc - len(zero_complexes),
           "n_blocked": n_blocked, "seconds_blocked": t_blocked,
           "n_balanced": result.n_balanced,
           "n_concordance_modules": sum(1 for g in result.modules if len(g) > 1),
           "n_concordant_pairs": result.n_concordant_pairs,
           "n_sccs": len(sccs), "largest_scc": max((len(c) for c in sccs), default=0),
           "n_terminal_sccs": terminal,
           "n_kinetic_modules": kin.n_modules, "giant_kinetic_module": kin.giant_module,
           "deficiency": kin.deficiency, "kinetic_iterations": kin.iterations,
           "acr": list(acr.acr), "acrr": [list(p) for p in acr.acrr],
           "acr_free": [m for m in acr.acr if m in free],
           "acrr_free": [list(p) for p in acr.acrr if p[0] in free and p[1] in free],
           "complexes": complexes}
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"{nc} complexes, {out['n_balanced']} balanced, {out['n_concordance_modules']} modules, "
          f"{out['n_concordant_pairs']} pairs, {out['n_kinetic_modules']} kinetic "
          f"(giant {out['giant_kinetic_module']}), delta={out['deficiency']}, "
          f"ACR {len(out['acr_free'])} free / {len(out['acr'])} all, "
          f"ACRR {len(out['acrr_free'])} free / {len(out['acrr'])} all, {min(times):.2f}s")


if __name__ == "__main__":
    main()
