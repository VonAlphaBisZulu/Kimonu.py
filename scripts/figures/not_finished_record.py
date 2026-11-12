"""Record of the runs behind the "not finished" cells of Fig. 4A, fig. S4 and Data S1, and of the node memory, for
the Methods ("Benchmark tools, hardware and software"). Writes, under KIMONU_OUT (records.py),

  results/cocoa_jl/runs_not_finished/runs_not_finished.tsv
      every run of COCOA.jl on a network and setting without a result (no seed finished), one row per seed, from the
      job logs v110_<array job>_<task>.log of results/cocoa_jl/efficient_false and efficient_true and from sacct
      (state, start, end, elapsed and time limit of the array task). The slurmstepd line of the log tells a run
      stopped at its time limit ("DUE TO TIME LIMIT") from one cancelled by hand. A seed without a log has no run.
  results/reference/trial_runs_2026-08/trial_runs.tsv
      the trial runs of the reference of August and September 2026, with the code before the adaptations of the
      fork (pull request 1 of the reference repository) and a budget of 24 h per network (REFERENCE_TIME_BUDGET
      86400 s, exhaustive SLURM job 521953, sampled with seed 1234 jobs 527897, 538873 and 546030), one row per
      network and mode, from the per-network logs copied into exhaustive/ and sampled/ of that folder. Pairs tested
      in the budget, candidate pairs, seconds, and a linear projection of the time to test every candidate pair
      (balanced stage plus concordance seconds times candidate pairs over pairs tested, in days, without the R
      stage). For the 11 networks of the paper, the wall time of the run reported (results/reference/exhaustive/
      and sampled/, the median of the three seeds for the sampled mode) stands beside it. The job of a sampled log
      is the job whose run time holds the file time of the log; a log without a concordance line was killed by the
      job's out-of-memory handler (the "Killed" lines of the job's own log).
  results/hardware/node_memory.tsv
      CPUs and memory of the nodes as SLURM configures them (scontrol show node, RealMemory in MiB, that is GB as
      SLURM counts them, 2^20 bytes per MB). memcheck_572322.log in the same folder holds the memory online on node02
      by the memory blocks of /sys/devices/system/memory (256.0 GiB) and the MemTotal of the kernel.

It runs on the cluster only, because sacct and scontrol are not part of the archive.

    python figures/not_finished_record.py        # run in scripts/
"""
import glob, json, os, re, statistics, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import records

FIELDS = "JobID,State,Start,End,Elapsed,Timelimit,NodeList"
SEEDS = (1234, 2026, 31337)
TRIAL = {"exhaustive": ["521953"], "sampled": ["527897", "538873", "546030"]}


def sacct(jobs):
    """{JobID: row} of the jobs and their array tasks, allocations only"""
    txt = subprocess.run(["sacct", "-X", "-j", ",".join(sorted(set(jobs))), "-o", FIELDS, "-P", "-n"],
                         capture_output=True, text=True, check=True).stdout
    return {p[0]: dict(zip(FIELDS.split(","), p)) for p in (l.split("|") for l in txt.splitlines())}


def days(s):
    d, _, t = s.rpartition("-")
    h, m, sec = (int(x) for x in t.split(":"))
    return round(int(d or 0) + (h * 3600 + m * 60 + sec) / 86400, 2)


def write_tsv(path, cols, rows):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join("" if r.get(c) is None else str(r[c]) for c in cols) + "\n")
    print(f"wrote {path}, {len(rows)} rows")


def cocoa():
    rows = []
    for eff, part, prefix in ((False, records.COCOA_FALSE, "cocoa_v110"), (True, records.COCOA_TRUE, "cocoa_eff")):
        d = records.read(part)
        done = {re.match(rf"{prefix}_(.+)_HiGHS_p64_seed\d+\.json", os.path.basename(f)).group(1)
                for f in glob.glob(f"{d}/{prefix}_*_HiGHS_p64_seed*.json")}
        logs = {}                                                    # (model, seed) -> [(array task, header, text)]
        for f in sorted(glob.glob(f"{d}/v110_*_*.log")):
            t = open(f).read()
            h = re.match(r"=== (\S+) seed (\d+) \| .*\| (\S+)$", t.split("\n", 1)[0])
            if h:
                task = "_".join(os.path.basename(f)[5:-4].split("_"))
                logs.setdefault((h.group(1), int(h.group(2))), []).append((task, os.path.basename(f), t))
        missing = [m for m in records.models() if m not in done]
        S = sacct([task for m in missing for s in SEEDS for task, *_ in logs.get((m, s), [])] or ["0"])
        for m in missing:
            for s in SEEDS:
                if (m, s) not in logs:
                    rows.append(dict(tool="COCOA.jl v1.1.0", efficient=str(eff).lower(), network=m, seed=s, stopped_by="no run"))
                for task, name, t in logs.get((m, s), []):
                    a = S[task]
                    c = re.search(r"JOB (\d+) ON (\S+) CANCELLED AT (\S+)( DUE TO TIME LIMIT)?", t)
                    rows.append(dict(tool="COCOA.jl v1.1.0", efficient=str(eff).lower(), network=m, seed=s, log=f"{part}/{name}",
                                     slurm_array_task=task, slurm_job=c.group(1) if c else None, node=a["NodeList"],
                                     state=a["State"], start=a["Start"], end=a["End"], elapsed=a["Elapsed"],
                                     elapsed_days=days(a["Elapsed"]), time_limit=a["Timelimit"],
                                     stopped_by=("time limit" if c and c.group(4) else "cancelled by hand" if c else "finished")))
    cols = ["tool", "efficient", "network", "seed", "stopped_by", "elapsed_days", "elapsed", "time_limit", "state", "start",
            "end", "node", "slurm_array_task", "slurm_job", "log"]
    write_tsv(records.write(records.COCOA_NOT_FINISHED, "runs_not_finished.tsv"), cols, rows)


def reference():
    S = sacct([j for js in TRIAL.values() for j in js])
    stamp = lambda s: subprocess.run(["date", "-d", s, "+%s"], capture_output=True, text=True, check=True).stdout.strip()
    window = {j: (int(stamp(S[j]["Start"])), int(stamp(S[j]["End"])) + 60) for js in TRIAL.values() for j in js}
    rows = []
    for mode, jobs in TRIAL.items():
        d = records.read(records.REF_TRIAL, mode)
        killed = {m for j in jobs if os.path.exists(f"{d}/slurm_{j}.log")
                  for m in re.findall(r"Killed.*\n.*done (\S+) rc=137", open(f"{d}/slurm_{j}.log").read())}
        for f in sorted(glob.glob(f"{d}/reference_*.log")):
            m, t, mt = os.path.basename(f)[len("reference_"):-4], open(f).read(), os.path.getmtime(f)
            job = next(j for j in jobs if window[j][0] <= mt <= window[j][1])
            cx = re.search(r"complexes (\d+)", t)
            b = re.search(r"balanced: \d+ complexes in ([\d.]+) s", t)
            c = re.search(r"concordance \(\w+\): (\d+) of (\d+) pairs in ([\d.]+) s, completed=(\d)", t)
            r = dict(network=m, mode=mode, complexes=int(cx.group(1)) if cx else None, slurm_job=job, job_state=S[job]["State"],
                     job_elapsed=S[job]["Elapsed"], log=f"{records.REF_TRIAL}/{mode}/{os.path.basename(f)}",
                     balanced_seconds=float(b.group(1)) if b else None)
            if c:
                tested, cand, sec, comp = int(c.group(1)), int(c.group(2)), float(c.group(3)), c.group(4) == "1"
                r.update(pairs_tested=tested, candidate_pairs=cand, concordance_seconds=sec, completed=str(comp).lower(),
                         fraction_tested=round(tested / cand, 4),
                         projected_days=round(((r["balanced_seconds"] or 0) + sec * cand / tested) / 86400, 1))
            else:
                r.update(completed="false", note="killed by the out-of-memory handler" if m in killed else "no concordance line")
            if mode == "exhaustive":
                p = records.read(records.REF_EXHAUSTIVE, f"reference_{m}.json")
                if os.path.exists(p):
                    r["paper_run_days"] = round(json.load(open(p))["seconds"] / 86400, 2)
            else:
                sec = [json.load(open(p))["seconds"] for p in glob.glob(records.read(records.REF_SAMPLED, f"reference_{m}_seed*.json"))]
                if sec:
                    r["paper_run_days"] = round(statistics.median(sec) / 86400, 2)
            rows.append(r)
    cols = ["network", "mode", "complexes", "completed", "pairs_tested", "candidate_pairs", "fraction_tested", "balanced_seconds",
            "concordance_seconds", "projected_days", "paper_run_days", "note", "slurm_job", "job_state", "job_elapsed", "log"]
    write_tsv(records.write(records.REF_TRIAL, "trial_runs.tsv"), cols, sorted(rows, key=lambda r: (r["mode"], r["complexes"] or 0, r["network"])))


def node_memory():
    rows = []
    for n in sorted({f"node{i:02d}" for i in range(1, 17)}):
        t = subprocess.run(["scontrol", "show", "node", n], capture_output=True, text=True, check=True).stdout
        mib = int(re.search(r"RealMemory=(\d+)", t).group(1))
        rows.append(dict(node=n, cpus=re.search(r"CPUTot=(\d+)", t).group(1), real_memory_mib=mib, real_memory_gib=round(mib / 1024, 1)))
    write_tsv(records.write(records.HARDWARE, "node_memory.tsv"), ["node", "cpus", "real_memory_mib", "real_memory_gib"], rows)


if __name__ == "__main__":
    cocoa()
    reference()
    node_memory()
