#!/usr/bin/env python
"""Stage 3 of the benchmark: correctness and speed comparison of Kimonu against COCOA.jl.

Correctness is checked at the level of the *partition*, not of summary counts. Two runs can report
the same number of modules while disagreeing about which complexes belong together, so complexes
are matched by composition and the induced partitions are compared set-wise.

Singletons are reported separately throughout: COCOA labels a complex that is concordant with
nothing as -1, whereas Kimonu returns it as a one-element group. Counting them differently is the
easiest way to manufacture a disagreement that is not there.

Usage:
    python benchmarks/compare.py [--profile P] [--scenario i|ii]
"""
import argparse
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")


def load(tool, profile, scenario, model):
    path = os.path.join(RESULTS, f"{tool}_{profile}_{scenario}_{model}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def key(composition):
    return tuple(sorted((m, float(c)) for m, c in composition.items()))


def partition(run):
    """Return (balanced, groups) keyed by complex composition.

    ``groups`` is a set of frozensets and excludes singletons; ``balanced`` is a set of keys.
    """
    balanced, by_module = set(), {}
    for c in run["complexes"]:
        k = key(c["composition"])
        m = c["module"]
        if m == 0:
            balanced.add(k)
        elif m == -1:
            continue                       # singleton: concordant with nothing
        else:
            by_module.setdefault(m, set()).add(k)
    groups = {frozenset(g) for g in by_module.values() if len(g) > 1}
    return balanced, groups


def pairs(groups):
    return sum(len(g) * (len(g) - 1) // 2 for g in groups)


def concordant_pairs(groups):
    """The set of unordered concordant pairs induced by a grouping."""
    out = set()
    for g in groups:
        members = sorted(g)
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                out.add((members[i], members[j]))
    return out


def relation(a, b):
    """Describe how partition ``a`` relates to partition ``b`` at the level of pairs.

    Which direction a disagreement runs is the whole diagnosis: a method that misses concordances
    reports a strict subset of the pairs, whereas a method that invents them reports a strict
    superset. Reporting only "they differ" throws that away.
    """
    pa, pb = concordant_pairs(a), concordant_pairs(b)
    if pa == pb:
        return "identical"
    if pa > pb:
        return f"kimonu is a strict superset (+{len(pa - pb)} pairs)"
    if pa < pb:
        return f"kimonu is a strict subset (-{len(pb - pa)} pairs)"
    return f"neither contains the other (+{len(pa - pb)} / -{len(pb - pa)} pairs)"


def describe(k):
    """Readable form of a complex composition key."""
    return " + ".join(f"{c:g} {m}" if c != 1 else m for m, c in k)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="concordance")
    ap.add_argument("--scenario", default="i")
    ap.add_argument("--max-detail", type=int, default=5,
                    help="how many disagreeing complexes/groups to print per model")
    args = ap.parse_args()

    prefix = f"kimonu_{args.profile}_{args.scenario}_"
    pattern = os.path.join(RESULTS, prefix + "*.json")
    models = sorted(os.path.basename(p)[len(prefix):-len(".json")] for p in glob.glob(pattern))
    if not models:
        print(f"no Kimonu results for profile={args.profile} scenario={args.scenario}")
        return

    disagreements = []
    header = f"{'model':<24} {'complexes':>9} {'balanced':>16} {'modules':>14} {'pairs':>14} " \
             f"{'kimonu s':>9} {'cocoa s':>9} {'speedup':>8}  agreement"
    print(header)
    print("-" * len(header))

    for model in models:
        k = load("kimonu", args.profile, args.scenario, model)
        c = load("cocoa", args.profile, args.scenario, model)
        if c is None:
            print(f"{model:<24} (no COCOA result)")
            continue

        kb, kg = partition(k)
        cb, cg = partition(c)

        n_k, n_c = len(k["complexes"]), len(c["complexes"])
        if n_k != n_c:
            verdict = f"COMPLEX COUNT DIFFERS ({n_k} vs {n_c})"
        elif kb == cb and kg == cg:
            verdict = "exact match"
        else:
            bits = []
            if kb != cb:
                bits.append(f"balanced differ (+{len(kb - cb)}/-{len(cb - kb)})")
            if kg != cg:
                bits.append(f"modules differ (+{len(kg - cg)}/-{len(cg - kg)})")
            verdict = "MISMATCH: " + ", ".join(bits)

        speedup = c["seconds"] / k["seconds"] if k["seconds"] > 0 else float("inf")
        print(f"{model:<24} {f'{n_k}/{n_c}':>9} {f'{len(kb)}/{len(cb)}':>16} "
              f"{f'{len(kg)}/{len(cg)}':>14} {f'{pairs(kg)}/{pairs(cg)}':>14} "
              f"{k['seconds']:>9.2f} {c['seconds']:>9.2f} {speedup:>7.1f}x  {verdict}")

        if verdict != "exact match" and n_k == n_c:
            disagreements.append((model, kb, cb, kg, cg))

    print("\ncolumns are kimonu/cocoa; modules and pairs exclude singletons")

    for model, kb, cb, kg, cg in disagreements:
        print(f"\n--- {model} ---")
        print(f"concordant pairs: {relation(kg, cg)}")
        for label, only in (("balanced only in kimonu", kb - cb), ("balanced only in cocoa", cb - kb)):
            if only:
                print(f"{label} ({len(only)}):")
                for key_ in sorted(only)[:args.max_detail]:
                    print(f"    {describe(key_)}")
        merged = [g for g in kg - cg]
        if merged:
            print(f"groups only in kimonu ({len(merged)}), first {min(len(merged), args.max_detail)}:")
            for g in sorted(merged, key=lambda g: sorted(g))[:args.max_detail]:
                print("    {" + "; ".join(describe(x) for x in sorted(g)) + "}")
        split = [g for g in cg - kg]
        if split:
            print(f"groups only in cocoa ({len(split)}), first {min(len(split), args.max_detail)}:")
            for g in sorted(split, key=lambda g: sorted(g))[:args.max_detail]:
                print("    {" + "; ".join(describe(x) for x in sorted(g)) + "}")


if __name__ == "__main__":
    main()
