#!/usr/bin/env python
"""Run Kimonu, COCOA and the reference on the same networks and collect one comparable record.

Each cell is (tool, model, seed) and reports the same quantities plus its own stage timings, so a
disagreement can be traced to a stage rather than just a total. The three arms are deliberately not
given the same hardware: COCOA batches candidate pairs across Distributed workers, which is part of
its design, while the reference is serial by construction and has no parallel path at all. Core
counts are recorded per cell rather than equalised.

  kimonu     one process, SCIP for the blocked-reaction FVA, exact arithmetic elsewhere
  cocoa      HiGHS, `--procs` Distributed workers, one seed per cell
  reference  Octave (balanced + concordance) -> R (kinetic modules) -> R (ACR/ACRR tables)

The reference's ACR/ACRR tables additionally need a stoichiometric coupling matrix, which its
MATLAB workflow produces with F2C2. F2C2 does not run under Octave, so the matrices published with
the reference are used instead; that omits F2C2's runtime from the reference's total, which
understates it.
"""
import argparse
import json
import os
import shutil
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
HOME = os.path.expanduser("~")


def run(cmd, log_path, env=None, cwd=REPO, timeout=None):
    """Run one stage, keeping the full log. Returns (ok, seconds)."""
    e = dict(os.environ)
    e.update(env or {})
    t = time.perf_counter()
    with open(log_path, "w") as fh:
        try:
            p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=e, cwd=cwd,
                               timeout=timeout)
            ok = p.returncode == 0
        except subprocess.TimeoutExpired:
            ok = False
    return ok, time.perf_counter() - t


def kimonu_cell(args, model, outdir):
    out = os.path.join(outdir, f"kimonu_{model}.json")
    log = os.path.join(outdir, f"kimonu_{model}.log")
    ok, wall = run([args.python, os.path.join(HERE, "run_kimonu_crn.py"),
                    os.path.join(args.crn_dir, f"{model}_fixed.crn.json"), out,
                    "--repeats", "1", "--solver", args.kimonu_solver, "--assume-unblocked"],
                   log, env={"PYTHONPATH": args.straindesign}, timeout=args.timeout)
    rec = {"tool": "kimonu", "model": model, "ok": ok, "wall": wall, "cpus": 1,
           "solver": args.kimonu_solver, "log": log}
    if ok and os.path.isfile(out):
        d = json.load(open(out))
        rec.update(
            n_complexes=d.get("n_complexes"), n_reactions=d.get("n_reactions"),
            n_blocked=d.get("n_blocked"), balanced=d.get("n_balanced"),
            modules=d.get("n_concordance_modules"), pairs=d.get("n_concordant_pairs"),
            kinetic_modules=d.get("n_kinetic_modules"), giant=d.get("giant_kinetic_module"),
            deficiency=d.get("deficiency"),
            acr=len(d.get("acr", [])), acrr=len(d.get("acrr", [])),
            acr_free=len(d.get("acr_free", [])), acrr_free=len(d.get("acrr_free", [])),
            stages={"blocked_fva": d.get("seconds_blocked"),
                    "concordance": d.get("seconds_concordance"),
                    "kinetic": d.get("seconds_kinetic"), "acr": d.get("seconds_acr")})
    return rec


def cocoa_cell(args, model, seed, outdir):
    out = os.path.join(outdir, f"cocoa_{model}_seed{seed}.json")
    log = os.path.join(outdir, f"cocoa_{model}_seed{seed}.log")
    ok, wall = run(["julia", f"--project={os.path.join(HERE, 'julia_revision')}",
                    os.path.join(HERE, "run_cocoa_acr.jl"),
                    os.path.join(args.crn_dir, f"{model}_fixed.crn.json"), out, str(seed)],
                   log, env={"COCOA_SOLVER": args.cocoa_solver,
                             "COCOA_PROCS": str(args.cocoa_procs),
                             "JULIA_DEPOT_PATH": args.julia_depot},
                   timeout=args.timeout)
    rec = {"tool": "cocoa", "model": model, "seed": seed, "ok": ok, "wall": wall,
           "cpus": max(args.cocoa_procs, 1), "solver": args.cocoa_solver, "log": log}
    if ok and os.path.isfile(out) and os.path.getsize(out) > 0:
        d = json.load(open(out))
        rec.update(n_complexes=d.get("n_complexes"),
                   acr=len(d.get("acr", [])), acrr=len(d.get("acrr", [])),
                   giant=d.get("giant_module"),
                   stages={"concordance": d.get("seconds")})
    return rec


def reference_cell(args, model, outdir):
    """Octave balanced+concordance, then the reference's own R kinetic and table stages."""
    work = os.path.join(outdir, f"ref_{model}")
    for sub in ("Results/concordant", "Results/stoichiometric_coupling", "Results/Overview",
                "Results/MetSingle", "Results/MetDouble", "Results/Reactions_Giant", "Code"):
        os.makedirs(os.path.join(work, sub), exist_ok=True)
    shutil.copytree(os.path.join(args.upstream, "Code", "kinetic_modules"),
                    os.path.join(work, "Code", "kinetic_modules"), dirs_exist_ok=True)
    coupling = os.path.join(args.upstream, "Results", "stoichiometric_coupling",
                            f"{model}_stoichiometric_coupling_fixed.mat")
    if os.path.isfile(coupling):
        shutil.copy(coupling, os.path.join(work, "Results", "stoichiometric_coupling"))

    stem = f"{model}_concordant_fixed"
    mat = os.path.join(work, "Results", "concordant", f"{stem}.mat")
    jsn = os.path.join(outdir, f"reference_{model}.json")
    rec = {"tool": "reference", "model": model, "cpus": 1, "solver": "glpk", "stages": {}}

    # One wall budget, exhaustive throughout. The reference's own prefilter is stochastic and only
    # ever removes candidates, so a completed sampled run is a subset of the true concordance just
    # as a capped exhaustive run is -- same direction, without the nondeterminism. A capped result
    # is marked *: its concordance is a subset, so modules come out finer and pairs undercounted,
    # while the balanced set is unaffected.
    ok1, t1 = run([args.octave, "--no-gui",
                   os.path.join(HERE, "matlab_reference", "run_reference.m"),
                   os.path.join(args.crn_dir, f"{model}_fixed.crn.json"), jsn],
                  os.path.join(outdir, f"reference_{model}_octave.log"),
                  env={"REFERENCE_MAT_OUT": mat, "UPSTREAM_ALGORITHM": args.upstream,
                       "REFERENCE_MODE": "exhaustive",
                       "REFERENCE_TIME_BUDGET": str(args.budget),
                       "PATH": f"{os.path.dirname(args.octave)}:{os.environ['PATH']}"},
                  timeout=args.budget + 1800)
    d = json.load(open(jsn)) if (ok1 and os.path.isfile(jsn)) else {}
    rec["stages"]["octave_balanced_concordance"] = t1
    rec["completed"] = d.get("completed")
    rec["mode"] = "exhaustive" if d.get("completed") else "exhaustive*"
    rec["pairs_tested"] = d.get("pairs_tested")
    rec["candidate_pairs"] = d.get("candidate_pairs")
    if d:
        rec.update(balanced=d.get("n_balanced"), modules=d.get("n_modules"),
                   pairs=d.get("n_concordant_pairs"), n_complexes=d.get("n_complexes"))
        rec["stages"]["balanced"] = d.get("seconds_balanced")
        rec["stages"]["concordance"] = d.get("seconds_concordant")

    kdir = os.path.join(work, "Code", "kinetic_modules")
    ok2, t2 = run([args.rscript, "code_kineticModule_analysis.R", f"{stem}.mat"],
                  os.path.join(outdir, f"reference_{model}_rkinetic.log"), cwd=kdir,
                  timeout=args.timeout)
    rec["stages"]["r_kinetic"] = t2
    ok3, t3 = run([args.rscript, "Get_excel_tables_kinetic_modules.R", f"{stem}.RData"],
                  os.path.join(outdir, f"reference_{model}_rtables.log"), cwd=kdir,
                  timeout=args.timeout)
    rec["stages"]["r_tables"] = t3

    def count_rows(path):
        if not os.path.isfile(path):
            return None
        with open(path) as fh:
            return max(sum(1 for _ in fh) - 1, 0)   # drop the header

    rec["acr"] = count_rows(os.path.join(work, "Results", "MetSingle", f"MetSingle_{stem}.csv"))
    rec["acrr"] = count_rows(os.path.join(work, "Results", "MetDouble", f"MetDouble_{stem}.csv"))
    rec["ok"] = ok1 and ok2 and ok3
    rec["wall"] = t1 + t2 + t3
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("models", nargs="+")
    ap.add_argument("--out", default=os.path.join(HERE, "three_way"))
    ap.add_argument("--crn-dir", default=os.path.join(HERE, "unblocked"),
                    help="networks handed to every tool; blocked reactions already removed, "
                         "so the three arms differ by method rather than by input")
    ap.add_argument("--upstream", default=os.path.join(HOME, "pipelines/2026-08_Kimonu/Upstream_Algorithm"))
    ap.add_argument("--straindesign", default=os.path.join(HOME, "sd-scip-pricing"))
    ap.add_argument("--python", default=os.path.join(HOME, "miniconda3/envs/straindesign/bin/python"))
    ap.add_argument("--octave", default=os.path.join(HOME, "miniconda3/envs/octave/bin/octave"))
    ap.add_argument("--rscript", default=os.path.join(HOME, "miniconda3/envs/R/bin/Rscript"))
    ap.add_argument("--julia-depot", default="/scratch/sadreyev_lab/pschneider/.julia")
    ap.add_argument("--kimonu-solver", default="scip")
    ap.add_argument("--cocoa-solver", default="HiGHS")
    ap.add_argument("--cocoa-procs", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--tools", default="kimonu,cocoa,reference")
    ap.add_argument("--timeout", type=int, default=7200)
    ap.add_argument("--budget", type=int, default=86400,
                    help="wall seconds per stage; a capped reference run is marked *")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    tools = args.tools.split(",")
    records = []
    for model in args.models:
        for tool in tools:
            print(f"--- {tool} / {model}", flush=True)
            if tool == "kimonu":
                r = kimonu_cell(args, model, args.out)
            elif tool == "cocoa":
                r = cocoa_cell(args, model, args.seed, args.out)
            else:
                r = reference_cell(args, model, args.out)
            records.append(r)
            print(f"    ok={r['ok']} wall={r['wall']:.1f}s stages={r.get('stages')}", flush=True)
            with open(os.path.join(args.out, "three_way.json"), "w") as fh:
                json.dump(records, fh, indent=1)

    print(f"\n{'tool':10s} {'model':24s} {'cpus':>4} {'wall':>9} | {'bal':>6} {'mods':>6} "
          f"{'pairs':>9} {'kmod':>6} {'giant':>6} | {'ACR':>5} {'ACRR':>6} {'mode':>16}")
    for r in records:
        f = lambda k: "-" if r.get(k) is None else str(r.get(k))
        print(f"{r['tool']:10s} {r['model']:24s} {r['cpus']:>4} {r['wall']:>8.1f}s | "
              f"{f('balanced'):>6} {f('modules'):>6} {f('pairs'):>9} {f('kinetic_modules'):>6} "
              f"{f('giant'):>6} | {f('acr'):>5} {f('acrr'):>6} {f('mode'):>16}")


if __name__ == "__main__":
    main()
