#!/usr/bin/env python
"""Run the three implementations in stages, fastest model first, with a wall-clock cap per cell.

Order matters for turnaround: every stage is run over the whole model set before the next starts,
and the model order for a stage is the previous stage's measured runtime, ascending. That way the
cheap models finish first at every stage and a Kimonu error surfaces before hours are spent on
COCOA or the reference.

  stage 1  Kimonu       all models, natural order          -> sorted by Kimonu seconds
  stage 2  COCOA        that order, one run per seed       -> sorted by COCOA seconds
  stage 3  reference    that order

Each cell gets its own timeout, so a single hung model cannot hold up the sweep — the iIT341
reference run sat for 16 hours while five other models waited behind it.

Results are appended to the output JSON after every cell, so intermediate state is always readable.

Usage:
    python benchmarks/staged_sweep.py --upstream <dir> --out results.json
        [--models m1 m2 ...] [--variant fixed] [--seeds 1234 42 ...]
        [--timeout 3600] [--jobs 1] [--stages kimonu,cocoa,reference]
"""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


def discover(upstream, variant):
    """Models that have both an elementary-step network and a published concordance result."""
    mdir = os.path.join(upstream, "Models", "models_with_elementary_steps")
    rdir = os.path.join(upstream, "Results", "concordant")
    suffix = f"_pre_balanced_{variant}.mat"
    names = []
    for f in sorted(os.listdir(mdir)):
        if not f.endswith(suffix):
            continue
        name = f[: -len(suffix)]
        if os.path.isfile(os.path.join(rdir, f"{name}_concordant_{variant}.mat")):
            names.append(name)
    return names


def run(cmd, timeout, env=None, cwd=REPO):
    """Run one cell. Returns (status, seconds, stdout tail)."""
    t = time.perf_counter()
    e = dict(os.environ)
    e.update(env or {})
    try:
        p = subprocess.run(cmd, timeout=timeout, capture_output=True, text=True, env=e, cwd=cwd)
        el = time.perf_counter() - t
        status = "ok" if p.returncode == 0 else f"exit {p.returncode}"
        out = (p.stdout or "")[-2000:] + (p.stderr or "")[-2000:]
    except subprocess.TimeoutExpired:
        el = time.perf_counter() - t
        status, out = "timeout", ""
    return status, el, out


def cell_paths(args, name, stage, seed=None):
    tag = f"{name}_{args.variant}" + (f"_seed{seed}" if seed is not None else "")
    return os.path.join(args.workdir, f"{stage}_{tag}.json")


def kimonu_cell(args, name):
    crn = os.path.join(args.workdir, f"{name}_{args.variant}.crn.json")
    out = cell_paths(args, name, "kimonu")
    if not os.path.isfile(crn):
        status, el, log = run([sys.executable, os.path.join(HERE, "export_split_model.py"),
                               args.upstream, name, "--variant", args.variant, "--out", crn],
                              args.timeout)
        if status != "ok":
            return {"status": f"export {status}", "seconds": el, "log": log}
    status, el, log = run([sys.executable, os.path.join(HERE, "run_kimonu_crn.py"), crn, out,
                           "--repeats", str(args.repeats)], args.timeout)
    rec = {"status": status, "seconds": el, "log": log if status != "ok" else ""}
    if status == "ok" and os.path.isfile(out):
        with open(out) as fh:
            rec["result"] = json.load(fh)
            rec["seconds"] = rec["result"].get("seconds", el)
    return rec


def cocoa_cell(args, name, seed):
    crn = os.path.join(args.workdir, f"{name}_{args.variant}.crn.json")
    out = cell_paths(args, name, "cocoa", seed)
    status, el, log = run(["julia", f"--project={os.path.join(HERE, 'julia_revision')}",
                           os.path.join(HERE, "run_cocoa_acr.jl"), crn, out, str(seed)],
                          args.timeout)
    rec = {"status": status, "seconds": el, "seed": seed, "log": log if status != "ok" else ""}
    if status == "ok" and os.path.isfile(out):
        with open(out) as fh:
            rec["result"] = json.load(fh)
            rec["seconds"] = rec["result"].get("seconds", el)
    return rec


def reference_cell(args, name):
    model = os.path.join(args.upstream, "Models", "models_with_elementary_steps",
                         f"{name}_pre_balanced_{args.variant}.mat")
    published = os.path.join(args.upstream, "Results", "concordant",
                             f"{name}_concordant_{args.variant}.mat")
    status, el, log = run(["octave", "--no-gui",
                           os.path.join(HERE, "matlab_reference", "validate_against_published.m"),
                           model, published], args.timeout)
    return {"status": status, "seconds": el, "log": log[-4000:]}


def save(args, state):
    tmp = args.out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(state, fh, indent=1)
    os.replace(tmp, args.out)


def stage(args, state, label, jobs, tasks):
    """Run one stage's cells, saving after each, and return per-model seconds."""
    seconds = {}
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(fn): key for key, fn in tasks}
        for fut in as_completed(futures):
            key = futures[fut]
            rec = fut.result()
            state.setdefault(label, {}).setdefault(key[0], []).append(rec)
            seconds.setdefault(key[0], []).append(rec["seconds"] if rec["status"] == "ok" else None)
            save(args, state)
            print(f"  [{label}] {key[0]}{'' if key[1] is None else ' seed ' + str(key[1])}: "
                  f"{rec['status']} in {rec['seconds']:.1f}s", flush=True)
    return seconds


def order_by(seconds, previous):
    """Ascending by measured time; cells that failed or timed out go last, in previous order."""
    def key(name):
        vals = [v for v in seconds.get(name, []) if v is not None]
        return (0, max(vals)) if vals else (1, previous.index(name))
    return sorted(previous, key=key)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", default=os.path.join(HERE, "sweep"))
    ap.add_argument("--models", nargs="*", default=None)
    ap.add_argument("--variant", default="fixed")
    ap.add_argument("--seeds", nargs="*", type=int, default=[1234, 42, 7, 999, 2026, 31337])
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--stages", default="kimonu,cocoa,reference")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    models = args.models or discover(args.upstream, args.variant)
    stages = args.stages.split(",")
    state = {"models": models, "variant": args.variant, "seeds": args.seeds,
             "timeout": args.timeout, "order": {}}
    if os.path.isfile(args.out):
        with open(args.out) as fh:
            state = json.load(fh)
            state["models"] = models
    print(f"{len(models)} models, variant {args.variant}, "
          f"{args.timeout}s cap per cell, {args.jobs} concurrent\n", flush=True)

    order = models
    if "kimonu" in stages:
        print("stage 1: Kimonu", flush=True)
        sec = stage(args, state, "kimonu", args.jobs, [((m, None), (lambda m=m: kimonu_cell(args, m))) for m in order])
        order = order_by(sec, order)
        state["order"]["after_kimonu"] = order
        save(args, state)
        print(f"\n  order after Kimonu: {', '.join(order)}\n", flush=True)

    if "cocoa" in stages:
        print("stage 2: COCOA", flush=True)
        tasks = [((m, s), (lambda m=m, s=s: cocoa_cell(args, m, s))) for m in order for s in args.seeds]
        sec = stage(args, state, "cocoa", args.jobs, tasks)
        order = order_by(sec, order)
        state["order"]["after_cocoa"] = order
        save(args, state)
        print(f"\n  order after COCOA: {', '.join(order)}\n", flush=True)

    if "reference" in stages:
        print("stage 3: reference", flush=True)
        stage(args, state, "reference", args.jobs, [((m, None), (lambda m=m: reference_cell(args, m))) for m in order])
        save(args, state)


if __name__ == "__main__":
    main()
