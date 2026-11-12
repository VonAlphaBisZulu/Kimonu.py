#!/usr/bin/env python
"""Validate a Kimonu sweep against the reference implementation's published results.

Three layers are compared, each in the form that makes the direction of a disagreement visible
rather than just its size.

  balanced   set equality, and which side holds the extras
  modules    whether each published class sits inside a single Kimonu module. Kimonu grouping two
             published classes together is a coarsening, which is what the reference's missed
             concordances look like; a published class split across Kimonu modules is a genuine
             conflict, since the two implementations then disagree about a specific pair.
  pairs      concordant pairs implied by each partition, and how many published pairs Kimonu also
             calls concordant

The published files store the partition in `class_with_balanced`, one index array per class, with
singleton classes omitted and the balanced complexes present as a single class. `group` is a sign
class (U/B/P/N), not a module id.

Usage:
    python benchmarks/validate_sweep.py --upstream <dir> --sweep <dir> [--variant fixed]
        [--out validation.json] [models ...]
"""
import argparse
import collections
import json
import os

import numpy as np
import scipy.io as sio


def published(path):
    """(balanced set, list of module sets) from a published concordance result, 1-based."""
    d = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    balanced = {int(x) for x in np.atleast_1d(d["B"]) if str(x).lstrip("-").isdigit()}
    classes = [set(np.atleast_1d(x).astype(int).tolist()) for x in np.atleast_1d(d["class_with_balanced"])]
    return balanced, [c for c in classes if c != balanced]


def kimonu(path):
    """(balanced set, {complex: module id}, list of module sets) from a Kimonu result, 1-based."""
    d = json.load(open(path))
    module_of = {i + 1: c["module"] for i, c in enumerate(d["complexes"])}
    balanced = {i for i, m in module_of.items() if m == 0}
    groups = collections.defaultdict(set)
    for i, m in module_of.items():
        if m > 0:
            groups[m].add(i)
    return balanced, module_of, list(groups.values()), d


def pairs_of(groups):
    return sum(len(g) * (len(g) - 1) // 2 for g in groups)


def compare(kim_path, pub_path):
    kb, kof, kmods, raw = kimonu(kim_path)
    pb, pmods = published(pub_path)

    inside = split = 0
    conflicts = []
    for c in pmods:
        labels = {kof.get(i) for i in c}
        if len(labels) == 1:
            inside += 1
        else:
            split += 1
            if len(conflicts) < 5:
                conflicts.append(sorted(c)[:6])

    ppairs = set()
    for c in pmods:
        s = sorted(c)
        for a in range(len(s)):
            for b in range(a + 1, len(s)):
                ppairs.add((s[a], s[b]))
    shared = sum(1 for a, b in ppairs
                 if kof.get(a) is not None and kof.get(a) == kof.get(b) and kof.get(a) != -1)

    return {
        "balanced_kimonu": len(kb), "balanced_published": len(pb),
        "balanced_identical": kb == pb,
        "balanced_kimonu_only": len(kb - pb), "balanced_published_only": len(pb - kb),
        "modules_kimonu": len(kmods), "modules_published": len(pmods),
        "published_classes_inside_one_kimonu_module": inside,
        "published_classes_split_by_kimonu": split,
        "example_conflicts": conflicts,
        "pairs_kimonu": pairs_of(kmods), "pairs_published": len(ppairs),
        "published_pairs_also_kimonu": shared,
        "n_blocked": raw.get("n_blocked"),
        "seconds": raw.get("seconds"), "seconds_blocked": raw.get("seconds_blocked"),
        "seconds_concordance": raw.get("seconds_concordance"),
        "seconds_kinetic": raw.get("seconds_kinetic"), "seconds_acr": raw.get("seconds_acr"),
        "n_complexes": raw.get("n_complexes"), "n_reactions": raw.get("n_reactions"),
        "kinetic_modules": raw.get("n_kinetic_modules"), "giant": raw.get("giant_kinetic_module"),
        "deficiency": raw.get("deficiency"),
        "acr": len(raw.get("acr", [])), "acr_free": len(raw.get("acr_free", [])),
        "acrr": len(raw.get("acrr", [])), "acrr_free": len(raw.get("acrr_free", [])),
        "sccs": raw.get("n_sccs"), "largest_scc": raw.get("largest_scc"),
        "terminal_sccs": raw.get("n_terminal_sccs"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", required=True)
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--variant", default="fixed")
    ap.add_argument("--out", default=None)
    ap.add_argument("models", nargs="*")
    args = ap.parse_args()

    pubdir = os.path.join(args.upstream, "Results", "concordant")
    names = args.models or sorted(
        f[len("kimonu_"):-len(f"_{args.variant}.json")]
        for f in os.listdir(args.sweep)
        if f.startswith("kimonu_") and f.endswith(f"_{args.variant}.json"))

    rows = {}
    for name in names:
        kim = os.path.join(args.sweep, f"kimonu_{name}_{args.variant}.json")
        pub = os.path.join(pubdir, f"{name}_concordant_{args.variant}.mat")
        if not (os.path.isfile(kim) and os.path.isfile(pub)):
            continue
        try:
            rows[name] = compare(kim, pub)
        except Exception as exc:
            rows[name] = {"error": f"{type(exc).__name__}: {exc}"}

    head = (f"{'model':28s} {'sec':>7} {'blk':>5} | {'bal K':>6} {'bal P':>6} {'=':>3} | "
            f"{'mod K':>6} {'mod P':>6} {'in':>6} {'split':>5} | {'pairs K':>9} {'pairs P':>8} {'shared':>7}")
    print(head)
    print("-" * len(head))
    ok = split_free = 0
    for name in sorted(rows, key=lambda n: rows[n].get("seconds") or 0):
        r = rows[name]
        if "error" in r:
            print(f"{name:28s} {r['error']}")
            continue
        ok += r["balanced_identical"]
        split_free += (r["published_classes_split_by_kimonu"] == 0)
        print(f"{name:28s} {r['seconds']:7.2f} {r['n_blocked'] or 0:5d} | "
              f"{r['balanced_kimonu']:6d} {r['balanced_published']:6d} "
              f"{'yes' if r['balanced_identical'] else 'NO':>3} | "
              f"{r['modules_kimonu']:6d} {r['modules_published']:6d} "
              f"{r['published_classes_inside_one_kimonu_module']:6d} "
              f"{r['published_classes_split_by_kimonu']:5d} | "
              f"{r['pairs_kimonu']:9d} {r['pairs_published']:8d} {r['published_pairs_also_kimonu']:7d}")
    n = sum(1 for r in rows.values() if "error" not in r)
    print(f"\n{n} models | balanced identical to published on {ok} | "
          f"no published class split by Kimonu on {split_free}")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(rows, fh, indent=1)
        print(f"written {args.out}")


if __name__ == "__main__":
    main()
