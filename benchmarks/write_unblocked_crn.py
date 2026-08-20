#!/usr/bin/env python
"""Write a crn.json with the blocked reactions removed, so every tool can be given the same network.

Kimonu drops blocked reactions internally because its kernel test is only equivalent to the
feasible-set test on an unblocked network. Handing the other tools the same reduced network is what
makes a disagreement attributable to the method rather than to the input.
"""
import argparse, json, sys
import numpy as np
import cobra
cobra.Configuration().solver = "glpk"
import straindesign as sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("crn_in")
    ap.add_argument("crn_out")
    ap.add_argument("--solver", default="scip")
    args = ap.parse_args()

    crn = json.load(open(args.crn_in))
    nr, nm, nc = crn["n_reactions"], crn["n_metabolites"], crn["n_complexes"]

    model = cobra.Model("split")
    mets = [cobra.Metabolite(m) for m in crn["metabolites"]]
    model.add_metabolites(mets)
    terms = {}
    for i, j, v in crn["S"]:
        terms.setdefault(int(j) - 1, {})[int(i) - 1] = float(v)
    rxns = [cobra.Reaction(f"r{j}") for j in range(nr)]
    for j, r in enumerate(rxns):
        r.lower_bound = float(crn["lb"][j]); r.upper_bound = float(crn["ub"][j])
    model.add_reactions(rxns)
    for j, t in terms.items():
        rxns[j].add_metabolites({mets[i]: c for i, c in t.items()})

    ranges = sd.fva(model, solver=args.solver)
    keep = [j for j in range(nr)
            if not (abs(ranges.iat[j, 0]) <= 1e-9 and abs(ranges.iat[j, 1]) <= 1e-9)]
    pos = {j: k for k, j in enumerate(keep)}

    out = dict(crn)
    out["n_reactions"] = len(keep)
    out["reactions"] = [crn["reactions"][j] for j in keep]
    out["lb"] = [crn["lb"][j] for j in keep]
    out["ub"] = [crn["ub"][j] for j in keep]
    for key in ("A", "S"):
        out[key] = [[i, pos[int(j) - 1] + 1, v] for i, j, v in crn[key] if int(j) - 1 in pos]
    # complexes left with no reaction would make A have empty rows; keep them so complex indices
    # stay comparable across tools, which is what the whole comparison is keyed on.
    json.dump(out, open(args.crn_out, "w"))
    print(f"{args.crn_in}: {nr} -> {len(keep)} reactions ({nr - len(keep)} blocked removed), "
          f"{nc} complexes kept -> {args.crn_out}")


if __name__ == "__main__":
    main()
