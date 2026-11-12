"""Hardware, memory, thread and software record of the benchmark runs, for the Methods ("Benchmark tools, hardware and
software", "Floating-point comparison"). Writes results/hardware/hardware_record.json (records.py, KIMONU_OUT) from
SLURM accounting (sacct), the job logs and the result files, so that every such number in the text has a source. It
runs on the cluster only, because sacct and the SLURM logs of the jobs (CLUSTER below) are not part of the archive.

    python figures/hardware_record.py            # run in scripts/

Sources (records of KIMONU_RECORDS unless named otherwise)
  COCOA.jl   results/cocoa_jl/efficient_false and efficient_true: the 183 records gather_figdata.py reads; each is
             matched to the job log (v110_<job>_<task>.log) with the same model and seed whose end time is closest to
             the file time of the record, and the wall time of that log must agree with the "seconds" of the record.
  reference  job ids from the log names in results/reference/exhaustive/ and sampled/
  Kimonu     results/kimonu/timing (hardware_*.txt, provenance_*.txt) and the SLURM logs of its jobs (CLUSTER
             kimonu_rerun_1111/logs, lines of /usr/bin/time -v)
  node speed CLUSTER kimonu_timing (job 570234, node16) against kimonu_timing_4314 (job 570464, node08), same commit 323fe40
  threads    CLUSTER kimonu_float_four/logs and results/kimonu/float_comparison/parts: the job that wrote each part and
             its allocated CPUs
A job that sacct no longer holds is recorded with "sacct": null. MaxRSS of sacct is in units of 1024 bytes and is reported in GB as SLURM counts them (2^30 bytes); it is
sampled every 30 s and is missing (a few hundred KiB) for many batch steps, which the record counts.
CPU models were read from /proc/cpuinfo of every node on 2026-11-12 (CPU_MODEL below); --probe reads them again by ssh.
"""
import collections, datetime, glob, json, os, re, statistics, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import records
CLUSTER = "/scratch/sadreyev_lab/pschneider"                    # SLURM logs and runs outside the archive
COCOA = {"cocoa_v110": records.read(records.COCOA_FALSE), "cocoa_eff": records.read(records.COCOA_TRUE)}
REF = {"reference_final": records.read(records.REF_EXHAUSTIVE), "reference_final_sampled": records.read(records.REF_SAMPLED)}
CPU_MODEL = {f"node{i:02d}": ("Xeon Silver 4314" if i <= 8 else "Xeon Silver 4514Y") for i in range(1, 17)}
FIELDS = "JobID,JobName,NodeList,ReqMem,ReqCPUS,AllocCPUS,MaxRSS,Elapsed,Start,End,State"

def sacct(jobs):
    """{JobID: row} for the jobs, their array tasks and steps"""
    out = {}
    jobs = sorted(set(jobs))
    for i in range(0, len(jobs), 50):
        txt = subprocess.run(["sacct", "-j", ",".join(jobs[i:i + 50]), "-o", FIELDS, "-P", "-n"], capture_output=True, text=True).stdout
        for line in txt.splitlines():
            p = line.split("|")
            out[p[0]] = dict(zip(FIELDS.split(","), p))
    return out

def gib(s):
    if not s: return None
    unit = {"K": 2**-20, "M": 2**-10, "G": 1.0}.get(s[-1])
    return round(float(s[:-1]) * unit, 2) if unit else round(float(s) * 2**-30, 2)

def hms(s):
    p = [float(x) for x in s.split(":")]
    return p[0] * 3600 + p[1] * 60 + p[2] if len(p) == 3 else p[0] * 60 + p[1]

def job_row(S, jid):
    """allocation of the job and the memory peak of its batch step"""
    a, b = S.get(jid), S.get(jid + ".batch")
    if a is None: return dict(job=jid, sacct=None)
    return dict(job=jid, name=a["JobName"], node=a["NodeList"], cpu=CPU_MODEL.get(a["NodeList"]), req_mem=a["ReqMem"], req_cpus=int(a["ReqCPUS"]),
                alloc_cpus=int(a["AllocCPUS"]), max_rss_gib=gib(b["MaxRSS"]) if b else None, elapsed=a["Elapsed"], state=a["State"].split()[0])

def cocoa():
    jobs = {re.match(r"v110_(\d+)_", os.path.basename(f)).group(1) for sub in ("cocoa_v110", "cocoa_eff") for f in glob.glob(f"{COCOA[sub]}/v110_*_*.log")}
    S = sacct(jobs)
    runs = []
    for sub, eff in (("cocoa_v110", False), ("cocoa_eff", True)):
        logs = collections.defaultdict(list)
        for f in glob.glob(f"{COCOA[sub]}/v110_*_*.log"):
            for blk in open(f).read().split("=== ")[1:]:
                h = re.match(r"(\S+) seed (\d+) \|.*\| (\d+) workers \| (\S+)\n", blk)
                el, rss, dn = re.search(r"Elapsed.*: (\S+)", blk), re.search(r"Maximum resident set size \(kbytes\): (\d+)", blk), re.search(r"done (\S+) ---", blk)
                if h and el and dn:
                    logs[h.group(1), h.group(2)].append(dict(job=os.path.basename(f)[5:-4], workers=int(h.group(3)), wall=hms(el.group(1)),
                                                             done=datetime.datetime.fromisoformat(dn.group(1)).timestamp(), main_rss_gib=round(int(rss.group(1)) * 2**-20, 2)))
        for j in sorted(glob.glob(f"{COCOA[sub]}/{sub}_*_HiGHS_p64_seed*.json")):
            m, seed = re.match(rf"{sub}_(.+)_HiGHS_p64_seed(\d+)\.json", os.path.basename(j)).groups()
            sec, mt = json.load(open(j))["seconds"], os.path.getmtime(j)
            x = min(logs[m, seed], key=lambda x: abs(x["done"] - mt))
            assert abs(x["done"] - mt) < 600 and sec <= x["wall"] <= 1.3 * sec + 120, (j, x)       # the run that wrote the record
            runs.append(dict(model=m, seed=int(seed), efficient=eff, seconds=round(sec, 1), workers=x["workers"], runs_in_logs=len(logs[m, seed]),
                             main_process_rss_gib=x["main_rss_gib"], **job_row(S, x["job"])))
    known = [r for r in runs if r.get("node")]
    new = [r for r in known if r["cpu"].endswith("4514Y")]
    rss = [r["max_rss_gib"] for r in known if r["max_rss_gib"] and r["max_rss_gib"] > 0.01]
    # every job of the sweeps, the runs that hit the 7-day limit included
    allpeak = max(((gib(v["MaxRSS"]) or 0, k, v["State"].split()[0], v["Elapsed"]) for k, v in S.items() if k.endswith(".batch")), default=None)
    both = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in known: both[r["efficient"], r["model"]][r["cpu"][-5:].strip()].append(r["seconds"])
    ratios = sorted(statistics.mean(v["4314"]) / statistics.mean(v["4514Y"]) for v in both.values() if len(v) == 2)
    return dict(n_runs=len(runs), n_without_sacct=len(runs) - len(known), n_on_4514Y=len(new), n_on_4314=len(known) - len(new),
                n_on_4514Y_by_efficient={str(e): sum(1 for r in new if r["efficient"] == e) for e in (False, True)},
                nodes=dict(collections.Counter(r["node"] for r in known)), req_mem=dict(collections.Counter(r["req_mem"] for r in known)),
                alloc_cpus=dict(collections.Counter(r["alloc_cpus"] for r in known)), workers=dict(collections.Counter(r["workers"] for r in runs)),
                n_with_max_rss=len(rss), peak_max_rss_gib_of_the_183=max(rss, default=None),
                peak_max_rss_gib_all_sweep_jobs=dict(gib=allpeak[0], step=allpeak[1], state=allpeak[2], elapsed=allpeak[3]) if allpeak else None,
                wall_ratio_4314_over_4514Y_same_model_other_seed=dict(n=len(ratios), median=round(statistics.median(ratios), 3), min=round(ratios[0], 3),
                                                                      max=round(ratios[-1], 3), note="different seeds, not the same computation") if ratios else None,
                runs=runs)

def reference():
    out = {}
    for d in ("reference_final", "reference_final_sampled"):
        jobs = {m.group(1) for f in glob.glob(f"{REF[d]}/*.log") for m in [re.search(r"_(\d{6})(?:_\d+)?\.log$", f)] if m}
        S = sacct(jobs)
        rows = [job_row(S, k) for k in sorted(S) if "." not in k]
        missing = sorted(j for j in jobs if not any(k.split("_")[0] == j for k in S))
        rss = [r["max_rss_gib"] for r in rows if r["max_rss_gib"]]
        out[d] = dict(jobs=sorted(jobs), jobs_without_sacct=missing, n_tasks=len(rows), nodes=dict(collections.Counter(r["node"] for r in rows)),
                      cpus=dict(collections.Counter(CPU_MODEL.get(r["node"]) for r in rows)), alloc_cpus=dict(collections.Counter(r["alloc_cpus"] for r in rows)),
                      req_mem_by_state={f"{k[0]} {k[1]}": v for k, v in sorted(collections.Counter((r["req_mem"], r["state"]) for r in rows).items())},
                      peak_max_rss_gib=max(rss, default=None))
    return out

def kimonu():
    out = {}
    T = records.read(records.TIMING)
    S = sacct({os.path.basename(f)[9:-4] for f in glob.glob(f"{T}/*/hardware_*.txt")})
    for f in sorted(glob.glob(f"{T}/*/hardware_*.txt")):
        arm, jid = f.split("/")[-2], os.path.basename(f)[9:-4]
        txt = open(f).read()
        prov = dict(l.split(" ", 1) for l in open(f.replace("hardware_", "provenance_")).read().splitlines() if " " in l)
        peak = 0
        for lg in glob.glob(f"{CLUSTER}/kimonu_rerun_1111/logs/*_{jid}.log"):                       # lines "<time> rc=0 <wall> <max RSS in KiB> <file>" of /usr/bin/time
            for l in open(lg):
                p = l.split()
                if len(p) >= 5 and p[1].startswith("rc=") and p[3].isdigit(): peak = max(peak, int(p[3]))
        out[arm] = dict(host=txt.split()[0], cpu_lscpu=re.search(r"Model name:\s+(.+)", txt).group(1), provenance=prov,
                        peak_rss_gib_time_v_whole_job=round(peak * 2**-20, 2) or None, **job_row(S, jid))
    return out

def node_speed():
    """the same Kimonu computation (commit 323fe40, 1 core, exclusive node) on node16 (4514Y) and node08 (4314)"""
    out = dict(jobs={k: job_row(sacct(["570234", "570236", "570464", "570492"]), k) for k in ("570234", "570236", "570464", "570492")})
    for sub in ("det", "att", "ref", "tanatt"):
        r = []
        for f in sorted(glob.glob(f"{CLUSTER}/kimonu_timing_4314/{sub}/*.json")):
            g = f.replace("kimonu_timing_4314", "kimonu_timing")
            if os.path.exists(g): r.append(json.load(open(f))["seconds"] / json.load(open(g))["seconds"])
        if r:
            r.sort()
            out[sub] = dict(n=len(r), median=round(statistics.median(r), 3), q10=round(r[int(.1 * (len(r) - 1))], 3), q90=round(r[int(.9 * (len(r) - 1))], 3),
                            min=round(r[0], 3), max=round(r[-1], 3))
    return out

def float_threads():
    """threads = allocated CPUs of the job that wrote each part (OMP, OPENBLAS and MKL threads are set to SLURM_CPUS_PER_TASK)"""
    L = f"{CLUSTER}/kimonu_float_four"
    wrote = collections.defaultdict(list)
    for f in glob.glob(f"{L}/logs/k[fx]_*.log"):
        for l in open(f, errors="ignore"):
            m = re.search(r"wrote (\S+)", l)
            if m: wrote[os.path.basename(m.group(1))].append((os.path.getmtime(f), os.path.basename(f)[3:-4]))
    S = sacct({j.split("_")[0] for v in wrote.values() for _, j in v})
    arms, models = collections.defaultdict(collections.Counter), collections.defaultdict(lambda: collections.defaultdict(set))
    for p in sorted(os.path.basename(x) for x in glob.glob(records.read(records.FLOAT_PARTS, "*.json"))):
        model, arm = p[:-5].split("__"); arm = re.sub(r"_t?\d+$", "", arm)
        jid = max(wrote[p])[1] if p in wrote else None
        c = S[jid]["AllocCPUS"] if jid in S else "unknown"
        arms[arm][c] += 1; models[arm][c].add(model)
    return dict(parts_by_arm_and_threads={a: dict(v) for a, v in sorted(arms.items())},
                models_with_more_than_2_threads={a: {c: sorted(ms) for c, ms in v.items() if c != "2"} for a, v in sorted(models.items())},
                n_models_all_parts_on_2_threads={a: len(v["2"] - set().union(*[ms for c, ms in v.items() if c != "2"] or [set()])) for a, v in sorted(models.items())})

def versions():
    home = os.path.expanduser("~")
    run = lambda cmd: (lambda r: (r.stdout + r.stderr).strip().splitlines()[0] if (r.stdout + r.stderr).strip() else None)(subprocess.run(cmd, capture_output=True, text=True, shell=True))
    py = run(f"{home}/miniconda3/envs/straindesign/bin/python -c \"import sys,json,numpy,scipy,cobra,flint,sympy,pyscipopt,cplex; "
             "print(json.dumps(dict(python=sys.version.split()[0],numpy=numpy.__version__,scipy=scipy.__version__,cobra=cobra.__version__,"
             "python_flint=flint.__version__,flint=flint.__FLINT_VERSION__,sympy=sympy.__version__,pyscipopt=pyscipopt.__version__,"
             "scip=str(pyscipopt.Model().version()),cplex=cplex.Cplex().get_version())))\" 2>/dev/null")
    return dict(straindesign_env=json.loads(py) if py else None, note="the environment as it is now; the timing jobs record only Python, Kimonu and StrainDesign (provenance)",
                julia_manifest=re.search(r'julia_version = "(.+)"', open("julia_v110/Manifest.toml").read()).group(1),
                octave=run(f"{home}/miniconda3/envs/octave/bin/octave --version"), R=run(f"{home}/miniconda3/envs/R/bin/R --version"))

if __name__ == "__main__":
    os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    if "--probe" in sys.argv:
        for n in CPU_MODEL:
            r = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=6", n, "grep -m1 'model name' /proc/cpuinfo"], capture_output=True, text=True)
            CPU_MODEL[n] = r.stdout.split(":")[-1].strip() or CPU_MODEL[n]
    D = dict(generated=datetime.datetime.now().isoformat(timespec="seconds"), cpu_model=CPU_MODEL, cocoa=cocoa(), reference=reference(), kimonu_timing=kimonu(),
             node_speed_4314_over_4514Y=node_speed(), float_threads=float_threads(), versions=versions())
    json.dump(D, open(records.write(records.HARDWARE, "hardware_record.json"), "w"), indent=1)
    c = D["cocoa"]
    print(f"COCOA.jl {c['n_runs']} runs, {c['n_on_4514Y']} on 4514Y {c['n_on_4514Y_by_efficient']}, ReqMem {c['req_mem']}, CPUs {c['alloc_cpus']}, "
          f"peak of the 183 {c['peak_max_rss_gib_of_the_183']} GB ({c['n_with_max_rss']} with MaxRSS), all jobs {c['peak_max_rss_gib_all_sweep_jobs']}")
    for k, v in D["reference"].items(): print(k, v["cpus"], v["alloc_cpus"], v["req_mem_by_state"], "peak", v["peak_max_rss_gib"], "no sacct", v["jobs_without_sacct"])
    for k, v in D["kimonu_timing"].items(): print(k, v["host"], v["cpu_lscpu"], v.get("req_mem"), v.get("max_rss_gib"), v["peak_rss_gib_time_v_whole_job"])
    print({k: v for k, v in D["node_speed_4314_over_4514Y"].items() if k != "jobs"})
    print(D["float_threads"]["parts_by_arm_and_threads"], D["float_threads"]["n_models_all_parts_on_2_threads"])
    print(D["versions"])
