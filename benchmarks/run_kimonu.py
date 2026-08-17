#!/usr/bin/env python
"""Stage 2b of the benchmark: run Kimonu's exact nullspace concordance on the preprocessed model.

The model comes from ``prepare_models.jl`` / ``run_cocoa.jl``, so Kimonu analyses exactly the
network COCOA analysed — no preprocessing difference can leak into the comparison. Results are
written with complexes keyed by composition, matching the COCOA exporter.

Usage:
    python benchmarks/run_kimonu.py [--profile P] [--scenario i] [--repeats N] [--model M]...
"""
import argparse
import json
import os
import time

import cobra

import kimonu
from kimonu.io import build_crn

# No LP is solved here (preprocess=False), but cobra instantiates its configured solver with every
# model. Pin the free one so the benchmark never occupies a commercial licence seat.
cobra.Configuration().solver = "glpk"

HERE = os.path.dirname(os.path.abspath(__file__))
PROCESSED = os.path.join(HERE, "processed")
RESULTS = os.path.join(HERE, "results")

DEFAULT_MODELS = ["EnvZ_OmpR", "e_coli_core", "e_coli_core_anaerobic", "iAB_RBC_283"]


def load_model(path):
    """Rebuild a cobra model from the JSON the Julia side wrote."""
    with open(path) as fh:
        data = json.load(fh)

    model = cobra.Model(os.path.splitext(os.path.basename(path))[0])
    mets = {m: cobra.Metabolite(m) for m in data["metabolites"]}
    model.add_metabolites(list(mets.values()))

    reactions = []
    for r in data["reactions"]:
        rxn = cobra.Reaction(r["id"])
        rxn.lower_bound = float(r["lb"])
        rxn.upper_bound = float(r["ub"])
        reactions.append((rxn, {mets[m]: float(c) for m, c in r["stoichiometry"].items()}))
    model.add_reactions([r for r, _ in reactions])
    for rxn, stoich in reactions:
        rxn.add_metabolites(stoich)
    return model


def composition_key(complex_tuple):
    """Canonical key for a Kimonu complex — a sorted (metabolite, coefficient) tuple."""
    return tuple(sorted((str(m), float(c)) for m, c in complex_tuple))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="concordance")
    ap.add_argument("--scenario", default="i", choices=["i"],
                    help="scenario (ii) needs the FVA-restricted kernel — not implemented yet")
    ap.add_argument("--model", action="append", dest="models", default=None)
    ap.add_argument("--repeats", type=int, default=2,
                    help="report the fastest run; the Julia side needs this to shed JIT warmup")
    args = ap.parse_args()
    models = args.models or DEFAULT_MODELS

    os.makedirs(RESULTS, exist_ok=True)
    for name in models:
        path = os.path.join(PROCESSED, args.profile, name + ".json")
        if not os.path.exists(path):
            print(f"missing preprocessed model, skipping: {path}")
            continue

        model = load_model(path)

        # preprocess=False: the model is already normalized, unblocked and split by the Julia
        # pipeline. Running Kimonu's own FVA here would analyse a different network than COCOA did.
        times = []
        for _ in range(args.repeats):
            t0 = time.perf_counter()
            result = kimonu.concordance(model, preprocess=False)
            times.append(time.perf_counter() - t0)
        seconds = min(times)

        crn = build_crn(model)
        module_of = {}
        for module_id, group in enumerate(result.modules, start=1):
            for i in group:
                module_of[i] = module_id if len(group) > 1 else -1
        for i in result.balanced:
            module_of[i] = 0

        out = {
            "tool": "Kimonu.py",
            "model": name,
            "profile": args.profile,
            "scenario": args.scenario,
            "seconds": seconds,
            "seconds_cold": times[0],
            "repeats": args.repeats,
            "n_reactions": len(model.reactions),
            "n_metabolites": len(model.metabolites),
            "complexes": [
                {"id": f"complex_{i + 1}",
                 "module": module_of[i],
                 "composition": {m: float(c) for m, c in comp}}
                for i, comp in enumerate(crn.complexes)
            ],
        }
        dest = os.path.join(RESULTS, f"kimonu_{args.profile}_{args.scenario}_{name}.json")
        with open(dest, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"{name}: {result.n_complexes} complexes, {result.n_balanced} balanced, "
              f"{result.n_modules} groups, {result.n_concordant_pairs} pairs, {seconds:.2f}s")


if __name__ == "__main__":
    main()
