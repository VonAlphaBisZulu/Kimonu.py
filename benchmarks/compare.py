#!/usr/bin/env python
"""Stage 3 of the benchmark: correctness and speed comparison across implementations.

Compares Kimonu against COCOA.jl and against the MATLAB reference implementation
(ankueken/Upstream_Algorithm, run under Octave — see ../matlab_reference).

Correctness is checked at the level of the *partition*, not of summary counts. Two runs can report
the same number of modules while disagreeing about which complexes belong together, so complexes
are matched by composition and the induced partitions are compared set-wise. When they disagree the
direction is reported: a tool that misses concordances yields a strict subset of the pairs, one
that invents them a strict superset.

Singletons are excluded from module and pair counts throughout. COCOA labels a complex concordant
with nothing as -1, Kimonu returns it as a one-element group; counting them differently is the
easiest way to manufacture a disagreement that is not there.

Usage:
    python benchmarks/compare.py [--profile P] [--scenario i|ii] [--max-detail N]
"""
import argparse
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")

TOOLS = [("kimonu", "Kimonu"), ("cocoa", "COCOA.jl"), ("reference", "MATLAB ref")]


def load(tool, profile, scenario, model):
    path = os.path.join(RESULTS, f"{tool}_{profile}_{scenario}_{model}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def key(composition):
    return tuple(sorted((m, float(c)) for m, c in composition.items()))


def partition(run):
    """Return (balanced, groups) keyed by complex composition; groups excludes singletons."""
    balanced, by_module = set(), {}
    for c in run["complexes"]:
        k = key(c["composition"])
        m = c["module"]
        if m == 0:
            balanced.add(k)
        elif m == -1:
            continue
        else:
            by_module.setdefault(m, set()).add(k)
    return balanced, {frozenset(g) for g in by_module.values() if len(g) > 1}


def concordant_pairs(groups):
    out = set()
    for g in groups:
        members = sorted(g)
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                out.add((members[i], members[j]))
    return out


def relation(a, b, name_a, name_b):
    pa, pb = concordant_pairs(a), concordant_pairs(b)
    if pa == pb:
        return "identical"
    if pa > pb:
        return f"{name_a} is a strict superset of {name_b} (+{len(pa - pb)} pairs)"
    if pa < pb:
        return f"{name_a} is a strict subset of {name_b} (-{len(pb - pa)} pairs)"
    return f"neither contains the other (+{len(pa - pb)} / -{len(pb - pa)} pairs)"


def describe(k):
    return " + ".join(f"{c:g} {m}" if c != 1 else m for m, c in k)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="concordance")
    ap.add_argument("--scenario", default="i")
    ap.add_argument("--max-detail", type=int, default=5)
    args = ap.parse_args()

    prefix = f"kimonu_{args.profile}_{args.scenario}_"
    models = sorted(os.path.basename(p)[len(prefix):-len(".json")]
                    for p in glob.glob(os.path.join(RESULTS, prefix + "*.json")))
    if not models:
        print(f"no Kimonu results for profile={args.profile} scenario={args.scenario}")
        return

    header = (f"{'model':<24} {'tool':<11} {'complexes':>9} {'balanced':>8} {'modules':>7} "
              f"{'pairs':>7} {'seconds':>9}  vs Kimonu")
    print(header)
    print("-" * len(header))

    disagreements = []
    for model in models:
        runs = [(label, tool, load(tool, args.profile, args.scenario, model))
                for tool, label in TOOLS]
        runs = [(label, tool, r) for label, tool, r in runs if r is not None]
        if not runs:
            continue
        base_label, _, base = runs[0][0], runs[0][1], runs[0][2]
        kb, kg = partition(base)

        for label, tool, run in runs:
            b, g = partition(run)
            if tool == "kimonu":
                verdict = "-"
            elif len(run["complexes"]) != len(base["complexes"]):
                verdict = f"COMPLEX COUNT DIFFERS ({len(run['complexes'])})"
            elif b == kb and g == kg:
                verdict = "exact match"
            else:
                bits = []
                if b != kb:
                    bits.append(f"balanced +{len(b - kb)}/-{len(kb - b)}")
                if g != kg:
                    bits.append(f"modules +{len(g - kg)}/-{len(kg - g)}")
                verdict = "MISMATCH: " + ", ".join(bits)
                disagreements.append((model, base_label, label, kb, b, kg, g))
            print(f"{model:<24} {label:<11} {len(run['complexes']):>9} {len(b):>8} {len(g):>7} "
                  f"{len(concordant_pairs(g)):>7} {run['seconds']:>9.2f}  {verdict}")
        print()

    print("modules and pairs exclude singletons\n")

    for model, na, nb, kb, b, kg, g in disagreements:
        print(f"--- {model}: {nb} vs {na} ---")
        print(f"concordant pairs: {relation(g, kg, nb, na)}")
        for lab, only in ((f"balanced only in {nb}", b - kb), (f"balanced only in {na}", kb - b)):
            if only:
                print(f"{lab} ({len(only)}):")
                for k in sorted(only)[:args.max_detail]:
                    print(f"    {describe(k)}")
        for lab, only in ((f"groups only in {nb}", g - kg), (f"groups only in {na}", kg - g)):
            if only:
                shown = sorted(only, key=lambda s: sorted(s))[:args.max_detail]
                print(f"{lab} ({len(only)}), first {len(shown)}:")
                for grp in shown:
                    print("    {" + "; ".join(describe(x) for x in sorted(grp)) + "}")
        print()


if __name__ == "__main__":
    main()
