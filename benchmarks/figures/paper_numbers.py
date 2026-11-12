"""Every number of the Results, Materials and Methods and Supplementary Text, recomputed from the result files, each with
its source.

One line per number, in the order of the text: the section, a short key, what the number is, the value, and the file
or directory it comes from. After a rerun, point the variables below at the new result directories, run this script
and update the text where a value changed. Numbers that no result file holds (software versions, SLURM settings of the
other tools) are listed at the end with where they come from.

    python figures/paper_numbers.py        # from benchmarks/

Result directories (environment variable, default; RERUN = /scratch/sadreyev_lab/pschneider/kimonu_rerun_1111, the
runs of the paper on the final code):
    KIMONU_TIMING                     RERUN/timing   (det, att, ref, tan, tanatt, small; hardware_*.txt, slurm_jobs.txt)
    KIMONU_DET_DIR, KIMONU_REF_DIR    KIMONU_TIMING/det, KIMONU_TIMING/ref
    KIMONU_LACR                       /scratch/sadreyev_lab/pschneider/kimonu_local_acr   (symbolic baselines)
    KF_OUT                            /scratch/sadreyev_lab/pschneider/kimonu_float_four   (Fig. 3E, F, figs. S5, S6)
    GB_OUT                            RERUN/genericity_bound
    GAP_OUT                           RERUN/gap_roots
    SZ_DIR                            RERUN/sz_flux_33
    SCALE_DIR                         RERUN/scale_check   (S8, the rounding of the LP point)
    FVA_OUT                           /scratch/sadreyev_lab/pschneider/kimonu_timing_seeds/fva   (preprocessing time)
    UPTAKE_LISTS                      RERUN/uptake_lists
Repository files: figdata.json, figures/pair_counts.json, figures/pair_sets.json, figures/seed_spread.json,
kinetic_variants/uptake_modes/, member_counts/, fb_vs_ma/, second_order_results/, injectivity_results/,
descartes_*.json, one_lp_results.log, deficiency/deficiency_33.json, unblocked/, reference_final/, sweep/.

Complexes are counted with the zero complex and without the complexes that occur in no reaction (figdata.json
n_complexes, table S2). ACRR is counted as independent ratios, k - 1 per group of k metabolites (k free metabolites for
the free counts), from the group sizes of kinetic_variants/member_counts/. Coupled complex pairs stay pairs.
"""
import glob, json, os, re, statistics, sys

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(B)
sys.path.insert(0, os.path.join(B, "figures"))
SCR = "/scratch/sadreyev_lab/pschneider"
RERUN = f"{SCR}/kimonu_rerun_1111"
env = lambda k, d: os.environ.get(k, d)
TIMING = env("KIMONU_TIMING", f"{RERUN}/timing")
REF_DIR = env("KIMONU_REF_DIR", f"{TIMING}/ref")
DET_DIR = env("KIMONU_DET_DIR", f"{TIMING}/det")
LACR = env("KIMONU_LACR", f"{SCR}/kimonu_local_acr")
GB = env("GB_OUT", f"{RERUN}/genericity_bound")
GAP = env("GAP_OUT", f"{RERUN}/gap_roots")
SZ = env("SZ_DIR", f"{RERUN}/sz_flux_33")
SCALE = env("SCALE_DIR", f"{RERUN}/scale_check")
FVA = env("FVA_OUT", f"{SCR}/kimonu_timing_seeds/fva")
MODELS = [l.strip() for l in open("kinetic_variants/models33.txt") if l.strip()]
GAP_ORDER = ["iCN718", "iNJ661", "iHN637"]                     # the order of the text
med = statistics.median


def out(section, key, what, value, source):
    print(f"{section:8s} | {key} | {what} | {value} | {source}")


def pct(a, b):
    return f"{100 * a / b:.1f}%"


def maxrss_gb(paths):
    """largest 'Maximum resident set size' of /usr/bin/time -v logs, in GB (10^6 kB, as SLURM reports memory)"""
    kb = [int(m.group(1)) for p in paths for m in [re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(p).read())] if m]
    return max(kb) / 1e6 if kb else None


def times(arm):
    """{network: [seconds of each replicate]} of a timing arm"""
    t = {}
    for f in glob.glob(f"{TIMING}/{arm}/*_rep*.json"):
        r = json.load(open(f)); t.setdefault(r.get("network") or r.get("model") or re.sub(r"^(kimonu|tangent)_|_rep\d+\.json$", "", os.path.basename(f)), []).append(r["seconds"])
    return t


ratios = lambda sizes: sum(k - 1 for k in sizes)               # k - 1 independent ratios per ACRR group of k
pairs_of = lambda sizes: sum(k * (k - 1) // 2 for k in sizes)

# ---- inputs shared by several numbers
D = json.load(open("figdata.json"))
nC = D["n_complexes"]
UM = {}
for f in sorted(glob.glob("kinetic_variants/uptake_modes/*.json")):
    d = json.load(open(f))
    if "detached" in d:
        UM.setdefault(d["model"], d)
assert sorted(UM) == sorted(MODELS), len(UM)
det = {m: UM[m]["detached"] for m in MODELS}
att = {m: UM[m]["attached"] for m in MODELS}
it = {m: det[m]["it"] for m in MODELS}
gap = sorted(m for m in MODELS if det[m]["gap"])
MCd = {d["model"]: d for d in (json.load(open(f"kinetic_variants/member_counts/{m}.json")) for m in MODELS)}
MCa = {m: json.load(open(f"kinetic_variants/member_counts/{m}.attached.json")) for m in MODELS}
for m in MODELS:                                                 # the group sizes reproduce the pair counts
    assert pairs_of(MCd[m]["it"]["acrr_group_sizes"]) == it[m]["acrr_pairs"], m
    assert pairs_of(MCd[m]["it"]["acrr_free_group_sizes"]) == it[m]["acrr_pairs_free"], m
    assert MCd[m]["n_complexes"] - D["complexes_in_no_reaction"][m] == nC[m], m

# ======================================================================================== Results
R = "Results"
out(R, "R1", "networks", len(MODELS), "kinetic_variants/models33.txt")
out(R, "R2", "complexes per network (zero complex counted, complexes of no reaction left out, as table S2)",
    f"{min(nC.values()):,} to {max(nC.values()):,}", "figdata.json n_complexes")
out(R, "R2b", "complexes that occur in no reaction, left out", ", ".join(f"{m} {k}" for m, k in D["complexes_in_no_reaction"].items() if k),
    "figdata.json complexes_in_no_reaction")
out(R, "R3", "networks on which the routes agree", f"{len(MODELS) - len(gap)} of {len(MODELS)} (gap: {', '.join(gap)})", "uptake_modes/ detached gap")
v = sorted(x["pairs"] for x in it.values())
out(R, "R4", "coupled complex pairs per network, iterated route", f"{v[0]:,} to {v[-1]:,}, median {med(v):,}", "uptake_modes/ detached it.pairs")
v = [x["modules"] for x in it.values()]
out(R, "R5", "kinetic modules per network", f"{min(v):,} to {max(v):,}", "uptake_modes/ detached it.modules")
v = [x["giant"] for x in it.values()]
out(R, "R6", "complexes in the largest module", f"{min(v):,} to {max(v):,}", "uptake_modes/ detached it.giant")
out(R, "R7", "ACR metabolites, all 33", f"{sum(x['acr'] for x in it.values()):,}", "uptake_modes/ detached it.acr")
out(R, "R8", "ACRR ratios (k - 1 per group), all 33", f"{sum(ratios(MCd[m]['it']['acrr_group_sizes']) for m in MODELS):,} "
    f"(pairs {sum(x['acrr_pairs'] for x in it.values()):,}, groups {sum(len(MCd[m]['it']['acrr_group_sizes']) for m in MODELS):,})",
    "member_counts/ it.acrr_group_sizes")
out(R, "R9", "free metabolites with ACR, networks", f"{sum(x['acr_free'] for x in it.values())} on {sum(1 for x in it.values() if x['acr_free'])}",
    "uptake_modes/ detached it.acr_free")
out(R, "R10", "free ACRR ratios (k - 1 per group of k free metabolites)", f"{sum(ratios(MCd[m]['it']['acrr_free_group_sizes']) for m in MODELS):,} "
    f"(pairs {sum(x['acrr_pairs_free'] for x in it.values()):,}, groups {sum(len(MCd[m]['it']['acrr_free_group_sizes']) for m in MODELS):,})",
    "member_counts/ it.acrr_free_group_sizes")
ex = [MCd[m] for m in MODELS if not MCd[m]["gap"]]
cx = sorted(100 * d["it"]["complexes_in_modules"] / nC[d["model"]] for d in ex)
fr = sorted(100 * d["it"]["fixed_members_free"] / d["n_free"] for d in ex)
out(R, "R11", f"complexes in a kinetic module, {len(ex)} exact networks (of the complexes of table S2)",
    f"{cx[0]:.1f} to {cx[-1]:.1f}%, median {med(cx):.1f}%", "member_counts/ it.complexes_in_modules, figdata.json n_complexes")
out(R, "R12", f"free metabolites with ACR or in an ACRR group, {len(ex)} exact networks", f"{fr[0]:.1f} to {fr[-1]:.1f}%, median {med(fr):.1f}%",
    "member_counts/ it.fixed_members_free")

FM = [json.load(open(f"kinetic_variants/fb_vs_ma/{m}.json")) for m in MODELS]
S = lambda get: sum(get(d) for d in FM)
nr = S(lambda d: d["n_reactions"])
fb, sub = S(lambda d: d["flux_coupling"]["members"]), S(lambda d: d["flux_coupling_and_shared_substrate"]["members"])
ma, ma_t = S(lambda d: d["rates_iterated"]["members"]), S(lambda d: (d["rates_tangent"] or d["rates_iterated"])["members"])
out(R, "R13", "reactions fully coupled under flux balance", pct(fb, nr), "fb_vs_ma/ flux_coupling.members")
out(R, "R14", "reactions in a rate class under mass action", pct(ma, nr), "fb_vs_ma/ rates_iterated.members")
out(R, "R15", "of these, points from flux balance, from a shared substrate complex, from kinetic modules",
    f"{100 * fb / nr:.1f}, {100 * (sub - fb) / nr:.1f}, {100 * (ma - sub) / nr:.1f}", "fb_vs_ma/ flux_coupling_and_shared_substrate.members")
out(R, "R15b", "tangent route, reactions in a rate class (not in the text)", pct(ma_t, nr), "fb_vs_ma/ rates_tangent.members")

import make_data_s1_s2
n2, tot = make_data_s1_s2.data_s2(os.devnull)
out(R, "R16", "more free ACR metabolites with attached uptakes", tot["lost"], "UPTAKE_LISTS (uptake_lists.py), as Data S2")
out(R, "R17", "of these extracellular or periplasmic", f"{tot['extracellular'] + tot['periplasmic']} = {pct(tot['extracellular'] + tot['periplasmic'], tot['lost'])}",
    "UPTAKE_LISTS")

pairs = [det[m]["tan"]["pairs"] - det[m]["it"]["pairs"] for m in GAP_ORDER]
fixed_free = lambda r: r["acr_free"] + ratios(r["acrr_free_group_sizes"])
rel = [fixed_free(MCd[m]["tan"]) - fixed_free(MCd[m]["it"]) for m in GAP_ORDER]
rel_pairs = [det[m]["tan"]["acr_free"] + det[m]["tan"]["acrr_pairs_free"] - det[m]["it"]["acr_free"] - det[m]["it"]["acrr_pairs_free"] for m in GAP_ORDER]
out(R, "R18", f"tangent-only coupled complex pairs on {', '.join(GAP_ORDER)}", ", ".join(f"{x:,}" for x in pairs), "uptake_modes/ detached tan - it")
out(R, "R19", "tangent-only free ACR and ACRR relations (ACRR as k - 1 ratios), same networks", ", ".join(map(str, rel)) +
    f" (as pairs {', '.join(map(str, rel_pairs))})", "member_counts/ tan - it")
cases = [json.load(open(f)) for f in sorted(glob.glob(f"{GAP}/case_*.json"))]
roots = []
for m in GAP_ORDER:
    ok = [len(c["roots"]) for c in cases if c["model"] == m and c["detached"] and c["roots"] != ["none"] and c["closed"] == c["n_relations"]]
    roots.append(min(ok) if ok else "none closes all")
out(R, "R20", "open questions (fewest roots found that close every open relation), same networks", ", ".join(map(str, roots)), f"{GAP}/case_*.json")

# floating-point comparison (Fig. 3E, F)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import make_float_four as ff
_, axs = plt.subplots(1, 2)
cert = ff.draw_certificate(axs, ["exactcmp", "svdabs"])
plt.close("all")
k14 = min(range(len(ff.TAUS)), key=lambda k: abs(ff.TAUS[k] - 1e-14))
e, s = cert["exactcmp"], cert["svdabs"]
out(R, "R21", f"exact nullspace, tolerance {ff.TAUS[k14]:g}: networks on which the routes differ", f"{e['gap'][k14] + e['float_differ'][k14]} of {len(e['models'])}",
    f"{ff.KF}/parts (make_float_four.draw_certificate)")
out(R, "R22", "exact nullspace, most networks certified at a tolerance above 0", max(e["certified"]), f"{ff.KF}/parts")
out(R, "R23", "SVD nullspace, most networks certified at any tolerance", max(s["certified"] + [s["at_zero"][0]]), f"{ff.KF}/parts")
out(R, "R24", "SVD nullspace at tolerance 0: false certificates", f"{s['at_zero'][3]} of {len(s['models'])}", f"{ff.KF}/parts")

G = json.load(open(f"{GB}/summary.json"))
for key, name in (("bound_tangent", "tangent route"), ("bound_iterated", "iterated route")):
    out(R, "R25", f"Schwartz-Zippel bound of the uniform r* (bits {G['bits']}), {name}, over {G['networks']} networks",
        f"min {G[key]['min']:.2g}, median {G[key]['median']:.2g}, max {G[key]['max']:.2g} ({G[key]['argmax']})", f"{GB}/summary.json")

# timing and comparison with the other tools (Fig. 4)
kd = D["kimonu_det"]
out(R, "R26", "Kimonu default, largest network (median of three seeds)", f"{kd['iAF1260b']['seconds']:.1f} s (iAF1260b); slowest of all "
    f"{max(v['seconds'] for v in kd.values()):.1f} s ({max(kd, key=lambda m: kd[m]['seconds'])})", "figdata.json kimonu_det")
sec = lambda v: med(r["seconds"] for r in v) if isinstance(v, list) else v["seconds"]
for arm, name in (("ref_exh", "exhaustive reference"), ("ref_sampled", "sampled reference"), ("cocoa_false", "COCOA.jl efficient=false"),
                  ("cocoa_true", "COCOA.jl efficient=true")):
    su = [sec(v) / kd[m]["seconds"] for m, v in D[arm].items()]
    out(R, "R27", f"median speed-up against the {name}", f"{med(su):.2g} over {len(su)} networks (range {min(su):.2g} to {max(su):.2g})", "figdata.json")
tan, tanatt, dett = times("tan"), times("tanatt"), times("det")
m_tan = max(tan, key=lambda m: med(tan[m]))
out(R, "R28", "tangent route, longest run (median of three, detached)", f"{med(tan[m_tan]) / 60:.1f} min ({m_tan}, {med(tan[m_tan]):.0f} s)", f"{TIMING}/tan")

PS = json.load(open("figures/pair_sets.json"))
keys = ["balanced", "conc", "n_kinetic", "giant", "acr", "acrr", "acr_free", "acrr_free"]
same = sum(D["kimonu_ref"][m][k] == D["ref_exh"][m][k] for m in D["ref_exh"] for k in keys)
out(R, "R29", "reference-equivalent mode equals the exhaustive reference", f"{same} of {len(keys) * len(D['ref_exh'])} values on {len(D['ref_exh'])} networks",
    "figdata.json kimonu_ref, ref_exh")
r0 = {m: PS[m]["ref_exh"][0] for m in PS}
out(R, "R30", "reference pairs not in Kimonu's default", sum(r["pairs_not_in_kimonu"] for r in r0.values()), "figures/pair_sets.json ref_exh")
miss = {m: len(r["acr_not_in_kimonu"]) + len(r["acrr_not_in_kimonu"]) for m, r in r0.items()}
out(R, "R31", "reference ACR and ACRR relations not in Kimonu's default", ", ".join(f"{m} {k}" for m, k in miss.items() if k) or "none", "figures/pair_sets.json")
out(R, "R32", "networks of the exhaustive reference on which the routes agree", f"{sum(1 for m in D['ref_exh'] if not det[m]['gap'])} of {len(D['ref_exh'])}",
    "uptake_modes/ detached gap")


def partition(f):
    d = json.load(open(f))
    groups = {}
    for c in d["complexes"]:
        groups.setdefault(c["kinetic_module"], set()).add(c["id"])
    return {frozenset(g) for g in groups.values()}, set(d["acr"]), {frozenset(p) for p in d["acrr"]}


same = 0
for m in MODELS:
    reps = [partition(f) for f in sorted(glob.glob(f"{DET_DIR}/kimonu_{m}_rep*.json"))]
    same += len(reps) == 3 and all(r == reps[0] for r in reps)
out(R, "R33", "networks with identical modules, ACR and ACRR under three seeds of r* (timing runs, seeds 0, 1, 2)", f"{same} of {len(MODELS)}", DET_DIR)
n = nC["iLJ478"]
out(R, "R34", "Fig. 4F, iLJ478 complexes and complex pairs", f"{n:,} complexes, {n * (n - 1) / 2e6:.2f} M pairs", "figdata.json n_complexes")

# ======================================================================================== Materials and Methods
M = "Methods"
LP = [json.loads(l) for l in open("kinetic_variants/one_lp_results.log") if l.startswith("{")]
out(M, "M1", "one LP finds exactly the blocked reactions of FVA", f"{sum(r['lp_wrong'] == 0 and r['lp_blocked'] == r['fva_blocked'] for r in LP)} of {len(LP)}",
    "kinetic_variants/one_lp_results.log")
out(M, "M2", "time of the one LP", f"{min(r['lp_seconds'] for r in LP):.3f} to {max(r['lp_seconds'] for r in LP):.3f} s", "kinetic_variants/one_lp_results.log")
lb, ub = set(), set()
for m in MODELS:
    c = json.load(open(f"unblocked/{m}_fixed.crn.json")); lb |= set(c["lb"]); ub |= set(c["ub"])
out(M, "M3", "bounds on every reaction", f"lower {sorted(lb)}, upper {sorted(ub)}", "unblocked/*.crn.json")
w = [float(x) for x in re.findall(r"iAF1260b.*\n.*wall ([\d.]+) s", open(sorted(glob.glob(f"{FVA}/preprocess_*.log"))[-1]).read())]
out(M, "M4", "preprocessing iAF1260b, load to factorization", f"{med(w):.1f} s (runs {', '.join(f'{x:.1f}' for x in w)})", f"{FVA}/preprocess_*.log")
Z = [json.load(open(f)) for f in sorted(glob.glob(f"{SZ}/*.json"))]
lp = lambda z: z["lp"] if "lp" in z else z["default"]         # records before 2026-11-11 call the LP arm "default"
ri = [z["uniform"]["iterated"]["seconds"] / lp(z)["iterated"]["seconds"] for z in Z]
rt = [z["uniform"]["tangent"]["seconds"] / lp(z)["tangent"]["seconds"] for z in Z]
out(M, "M5", f"uniform 2^64 flux against the LP point, runtime ratio ({len(Z)} networks, one run each, not exclusive)",
    f"iterated {min(ri):.1f} to {max(ri):.1f} (median {med(ri):.1f}), tangent {min(rt):.1f} to {max(rt):.1f} (median {med(rt):.1f}); identical results "
    f"{sum(z['iterated_agree'] and z['tangent_agree'] for z in Z)} of {len(Z)}", SZ)
for key, name in (("bound_tangent", "tangent"), ("bound_iterated", "iterated")):
    out(M, "M6", f"Schwartz-Zippel bound per network, uniform r*, {name} route (largest)", f"{G[key]['max']:.2g} ({G[key]['argmax']})", f"{GB}/summary.json")
rows = [json.load(open(f)) for f in sorted(glob.glob(f"{GB}/*.json")) if not f.endswith("summary.json")]
out(M, "M6b", "rounds of the iterated route at the uniform r*, tests per network", f"rounds {sorted({r['rounds'] for r in rows})}, tests "
    f"{min(r['tests'] for r in rows):.2g} to {max(r['tests'] for r in rows):.2g}", f"{GB}/<model>.json")
for f in ("kinetic_variants/float_exactcmp.py", "kinetic_variants/float_four.py"):
    r = re.search(r"rnd >= (\d+)", open(f).read())
    out(M, "M7", "rounds of the iterated route in the floating-point comparison, at most", r.group(1) if r else "?", f)
for arm in ("cocoa_false", "cocoa_true", "ref_sampled"):
    out(M, "M8", f"seeds of {arm}", sorted({len(v) for v in D[arm].values()}), "figdata.json")
for m in ("iSB619", "iNF517"):
    out(M, "M9", f"exhaustive reference run time, {m}", f"{D['ref_exh'][m]['seconds'] / 86400:.1f} days", "figdata.json ref_exh")
hw = sorted(glob.glob(f"{TIMING}/*/hardware_*.txt"))
cpus = {re.search(r"Model name:\s*(.+)", open(f).read()).group(1).strip() for f in hw}
hosts = {open(f).readline().strip() for f in hw}
mem = {re.search(r"Mem:\s+(\d+)", open(f).read()).group(1) for f in hw}
cores = {re.search(r"^CPU\(s\):\s+(\d+)", open(f).read(), re.M).group(1) for f in hw}
out(M, "M10", "hardware of the timed Kimonu runs", f"{', '.join(sorted(cpus))}; {', '.join(sorted(cores))} logical cores; nodes {', '.join(sorted(hosts))}; "
    f"{', '.join(sorted(mem))} GB (free -g)", f"{TIMING}/*/hardware_*.txt")
req = sorted({x for f in glob.glob(f"{TIMING}/*/provenance_*.txt") for x in re.findall(r"^kimonu (\w+)", open(f).read(), re.M)})
out(M, "M11", "commit of the timed runs", ", ".join(c[:7] for c in req), f"{TIMING}/*/provenance_*.txt")
jobs = f"{TIMING}/slurm_jobs.txt"                              # sacct of the timing jobs, written when they finished
if os.path.exists(jobs):
    J = [l.split("|") for l in open(jobs) if l.strip() and not l.startswith("JobID")]
    out(M, "M12", "memory reserved for the timed jobs, cores per run", "; ".join(f"{j[1]} {j[2]}" for j in J) + "; one process pinned to one core", jobs)
logs = lambda arms: [p for a in arms for p in glob.glob(f"{TIMING}/{a}/*.log")]
out(M, "M13", "peak memory, iterated route and reference-equivalent mode", f"{maxrss_gb(logs(['det', 'att', 'ref'])):.2f} GB", f"{TIMING}/{{det,att,ref}}/*.log (/usr/bin/time -v)")
out(M, "M14", "peak memory, tangent route", f"{maxrss_gb(logs(['tan', 'tanatt'])):.2f} GB", f"{TIMING}/{{tan,tanatt}}/*.log")
import numpy as np
fits = {}
for arm in ("ref_exh", "ref_sampled", "cocoa_false", "cocoa_true", "kimonu_det", "kimonu_ref"):
    ms = [m for m in D["models"] if m in D[arm]]
    fits[arm] = (np.polyfit(np.log10([nC[m] for m in ms]), np.log10([sec(D[arm][m]) for m in ms]), 1)[0], len(ms))
    out(M, "M15", f"fit exponent of wall time on complexes, {arm} (slope of Fig. 4A)", f"{fits[arm][0]:.2f} over {fits[arm][1]} networks",
        "figdata.json (as make_figure3.py)")
pw = [fits[a][0] for a in ("ref_exh", "ref_sampled", "cocoa_false", "cocoa_true")]
out(M, "M16", "exponents of the pairwise tests; networks per fit", f"{min(pw):.1f} to {max(pw):.1f}; {min(v[1] for v in fits.values())} to "
    f"{max(v[1] for v in fits.values())}", "figdata.json")
out(M, "M17", "tangent route, longest run", f"{med(tan[m_tan]) / 60:.1f} min on {m_tan}", f"{TIMING}/tan")
LR = json.load(open(f"{LACR}/small.json")) + [json.load(open(f)) for f in sorted(glob.glob(f"{LACR}/ladder_*.json"))]
copies = sorted(int(r["network"][4:]) for r in LR if re.fullmatch(r"envz\d+", r["network"]))
small_ = [r["network"] for r in LR if not r["network"].startswith(("envz", "e_coli"))]
out(M, "M18", "symbolic baselines: small networks, EnvZ-OmpR copies, total", f"{len(small_)}, copies {copies}, total {len(LR)}", LACR)
for key, name in (("rref", "elimination"), ("minors", "minors")):
    done = [r for r in LR if r.get(key)]
    agree = all(r.get(f"agree_{key}") for r in done)
    out(M, "M19", f"symbolic {name} finished, agrees with the tangent route", f"{len(done)} of {len(LR)}, agree {agree}", LACR)
timeouts = sorted(int(r["network"][4:]) for r in LR if re.fullmatch(r"envz\d+", r["network"]) and r.get("minors_timeout") is not None)
out(M, "M20", "copies on which the minors reached the limit", timeouts, LACR)
kr = D["kimonu_ref"]
out(M, "M21", "reference-equivalent mode, longest run; default mode on the same networks", f"{max(v['seconds'] for v in kr.values()):.1f} s "
    f"({max(kr, key=lambda m: kr[m]['seconds'])}); default {max(kd[m]['seconds'] for m in kr):.1f} s", "figdata.json kimonu_ref, kimonu_det")

# ======================================================================================== Supplementary Text
X = "Supp"
# S3. Making the tools comparable
ref = json.load(open("reference_final/reference_iLJ478.json"))
log = open("sweep/cocoa_eff/cell_iLJ478_seed1234.log").read()
total = int(re.search(r"total_possible_pairs = (\d+)", log).group(1))
cand, trans = map(int, re.search(r"Processing: (\d+) candidates, (\d+) filtered by transitivity", log).groups())
skip_bal = total - ref["candidate_pairs"]
out(X, "S3a", "iLJ478: pairs the reference tests, of all complex pairs (the 3439 complexes the LP tools test, without the zero complex)",
    f"{ref['candidate_pairs'] / 1e6:.2f} M of {total / 1e6:.2f} M", "reference_final/reference_iLJ478.json, sweep/cocoa_eff/cell_iLJ478_seed1234.log")
out(X, "S3b", "iLJ478: COCOA.jl pairs left for LPs; skipped as balanced, by the spread, by transitivity (efficient=true, seed 1234)",
    f"{cand / 1e3:.0f} k; {skip_bal / 1e6:.2f} M, {(total - skip_bal - cand - trans) / 1e6:.2f} M, {trans / 1e6:.2f} M", "sweep/cocoa_eff/cell_iLJ478_seed1234.log")
k = D["kimonu_det"]["iLJ478"]
out(X, "S3c", "iLJ478 kinetic stage, Kimonu (iterated route, kinetic and ACR stages)", f"{k['sec_kin'] + k['sec_acr']:.1f} s", "figdata.json kimonu_det sec_kin + sec_acr")
out(X, "S3d", "iLJ478, reference-equivalent mode in all", f"{D['kimonu_ref']['iLJ478']['seconds']:.1f} s", "figdata.json kimonu_ref")
for arm in ("cocoa_false", "cocoa_true"):
    by_seed = ", ".join(f"seed {r['seed']} {r['conc']}" for r in D[arm]["iLJ478"])
    out(X, "S3e", f"iLJ478 concordance modules, {arm}, against the exact partition", f"{by_seed} against {D['kimonu_det']['iLJ478']['conc']}",
        "figdata.json")
sub = {arm: sum(all(r["pairs_not_in_kimonu"] == 0 and not r["acr_not_in_kimonu"] and not r["acrr_not_in_kimonu"] for r in PS[m][arm]) for m in PS)
       for arm in ("ref_exh", "cocoa_false", "cocoa_true")}
out(X, "S3f", "networks on which a tool's relations are a subset of Kimonu's default (of 11)", ", ".join(f"{a} {v}" for a, v in sub.items()), "figures/pair_sets.json")
SS = json.load(open("figures/seed_spread.json"))
for arm in ("cocoa_false", "cocoa_true", "ref_sampled", "kimonu_det"):
    v = SS[arm]
    out(X, "S3g", f"networks on which the coupled pairs change with the seed, {arm}", f"{sum(1 for x in v.values() if x['pairs'][0])} of {len(v)}",
        "figures/seed_spread.json pairs[0]")
# S4. Run records
it_max = max(dett, key=lambda m: med(dett[m]))
out(X, "S4a", "tangent route, longest (median of three, detached; attached in brackets)", f"{med(tan[m_tan]) / 60:.1f} min ({m_tan}); attached "
    f"{max(med(v) for v in tanatt.values()) / 60:.1f} min", f"{TIMING}/tan, tanatt")
out(X, "S4b", "iterated route, longest (median of three, detached)", f"{med(dett[it_max]):.1f} s ({it_max})", f"{TIMING}/det")
ratio = sorted(med(tan[m]) / med(dett[m]) for m in tan if m in dett)
out(X, "S4c", "tangent route over iterated route, per network (medians)", f"{ratio[0]:.1f} to {ratio[-1]:.0f}-fold, median {med(ratio):.0f}", f"{TIMING}/tan, det")
# S5. Attached and detached uptakes
out(X, "S5a", "free ACR metabolites attached, detached", f"{tot['attached']}, {tot['both'] + tot['gained']}", "UPTAKE_LISTS (make_data_s1_s2.data_s2)")
out(X, "S5b", "attached only; extracellular, periplasmic; in a detached ACRR group; gained by detaching",
    f"{tot['lost']}; {tot['extracellular']}, {tot['periplasmic']}; {tot['lost_in_acrr']}; {tot['gained']}", "UPTAKE_LISTS")
FA = [json.load(open(f"kinetic_variants/fb_vs_ma/{m}.attached.json")) for m in MODELS]
out(X, "S5c", "networks with complexes of fixed concentration, detached, attached (attached, the zero complex is fixed and not counted)",
    f"{sum(1 for d in FM if d['fixed_complexes_iterated'])} and {sum(1 for d in FA if d['fixed_complexes_iterated'] > 1)} of {len(MODELS)}",
    "fb_vs_ma/<model>[.attached].json fixed_complexes_iterated")
out(X, "S5d", "gap networks attached", ", ".join(sorted(m for m in MODELS if att[m]["gap"])), "uptake_modes/ attached gap")
# S6. The gap between the routes
out(X, "S6a", "open relations (tangent route only) on iCN718, iNJ661, iHN637, detached; attached",
    f"{', '.join(str(det[m]['n_relations']) for m in GAP_ORDER)}; {', '.join(str(att[m]['n_relations']) for m in GAP_ORDER)}", "uptake_modes/ n_relations")
out(X, "S6b", "dimensions they add to the iterated route's fixed combinations", ", ".join(str(det[m].get("g_rel", att[m].get("g_rel"))) for m in GAP_ORDER),
    "uptake_modes/ g_rel")
acr_open = {m: [r["c"] and list(r["c"])[0] for r in det[m]["relations"] if r["kind"] == "ACR"] for m in GAP_ORDER}
out(X, "S6c", "open ACR metabolites", "; ".join(f"{m}: {', '.join(v) or 'none'}" for m, v in acr_open.items()), "uptake_modes/ detached relations")
for c in cases:
    out(X, "S6d", f"roots {c['model']} {'detached' if c['detached'] else 'attached'}", f"{', '.join(c['roots'])}: closes {c['closed']} of {c['n_relations']}",
        f"{GAP}/case_*.json")
out(X, "S6e", "open relations in all, root questions", f"{sum(det[m]['n_relations'] for m in GAP_ORDER)}, {sum(r for r in roots if isinstance(r, int))}",
    "uptake_modes/, gap_roots/")
for m in GAP_ORDER:
    d = json.load(open(f"kinetic_variants/descartes_{m}.json"))
    out(X, "S6f", f"Descartes step, {m}: settled relations, seconds (detached; attached)",
        f"{len(d['detached']['firings'])}, {d['detached']['seconds_descartes']:.0f} s; {d['attached']['seconds_descartes']:.0f} s", f"kinetic_variants/descartes_{m}.json")
for m in GAP_ORDER:
    i, q = (json.load(open(f"kinetic_variants/injectivity_results/{m}_{x}.json")) for x in ("inject", "quotient"))
    out(X, "S6g", f"injectivity test and sign program for all relations, {m}", f"{i['status']} after {i['seconds'] / 3600:.1f} h; {q['status']} after "
        f"{q['seconds'] / 3600:.1f} h", "kinetic_variants/injectivity_results/")
AB = [json.loads(l) for l in open("kinetic_variants/autonomous_blocks_33.jsonl")]
out(X, "S6h", "autonomous blocks: relations proven on the gap networks", ", ".join(f"{r['model']} {'detached' if r['detach'] else 'attached'} "
    f"{r['tangent_only_relations_proven']}" for r in AB if r["model"] in GAP_ORDER), "kinetic_variants/autonomous_blocks_33.jsonl")
# S7. Nondegeneracy
for m in GAP_ORDER:
    d = json.load(open(f"kinetic_variants/injectivity_results/{m}_dims.json"))
    out(X, "S7a", f"tangent dimension at four interior fluxes, and of the sum of two, {m} (attached)", f"{d['dims_T']}, {d['dim_T0_plus_Tk']}",
        "kinetic_variants/injectivity_results/<model>_dims.json")
# S8. Genericity
for f in sorted(f for f in glob.glob(f"{SCALE}/*.json") if not f.endswith("_1e4.json")):
    d = json.load(open(f))
    for route in ("tangent", "iterated"):
        a, b, c = (d[f"{route}_{x}"] for x in ("scale_1e6", "scale_100", "unrounded"))
        out(X, "S8a", f"{d['model']}, {route} route, LP point at scale 10^6, at scale 100, unrounded (seconds, largest bits; same results as at 10^6)",
            f"{a['seconds']:.1f} s ({a['max_bits']} bits), {b['seconds']:.1f} s ({b['max_bits']} bits, same {b['same_as_1e6']}), "
            f"{c['seconds']:.1f} s ({c['max_bits']} bits, same {c['same_as_1e6']})", f"{f} (one run each, not exclusive)")
for m in sorted({os.path.basename(f).split("_rep")[0] for f in glob.glob(f"{SCALE}/exclusive/*_rep*.json")}):
    reps = [json.load(open(f)) for f in sorted(glob.glob(f"{SCALE}/exclusive/{m}_rep*.json"))]
    med = lambda key: statistics.median(r[key]["seconds"] for r in reps)
    out(X, "S8a", f"{m}, tangent route, LP point at scale 10^6 against unrounded, median of {len(reps)} (exclusive node, 1 core)",
        f"{med('tangent_scale_1e6'):.1f} s against {med('tangent_unrounded'):.1f} s (iterated {med('iterated_scale_1e6'):.1f} s against "
        f"{med('iterated_unrounded'):.1f} s); same results {all(r['tangent_unrounded']['same_as_1e6'] for r in reps)}",
        f"{SCALE}/exclusive/{m}_rep*.json (job 572320)")
for f in sorted(glob.glob(f"{SCALE}/*_1e4.json")):
    d = json.load(open(f))
    out(X, "S8b", f"{d['model']}, same results as at scale 10^6 at scales 10^4 and 10^5 (tangent, iterated)",
        ", ".join(f"{k.split('_')[0]} {k.split('_')[-1]}: {v['same_as_1e6']}" for k, v in d.items() if k != "model"), f)
# S10. Deficiency
F = json.load(open("deficiency/deficiency_33.json"))
out(X, "S10a", "deficiency", f"{min(r['delta'] for r in F)} to {max(r['delta'] for r in F)}", "deficiency/deficiency_33.json")
# S11. Symbolic and floating-point baselines
env_ = [r for r in LR if re.fullmatch(r"envz\d+", r["network"])]
eco = [r for r in LR if r["network"].startswith("e_coli")][0]
out(X, "S11a", "EnvZ-OmpR copies, species; e_coli_core species, complexes", f"{min(r['species'] for r in env_)} to {max(r['species'] for r in env_)}; "
    f"{eco['species']}, {eco['complexes']}", LACR)
mi = [r for r in env_ if r.get("minors")]
el = [r for r in env_ if r.get("rref")]
out(X, "S11b", "minors finish on the copies (species, seconds)", ", ".join(f"{r['species']} in {r['minors']['seconds']:.1f} s" for r in mi), LACR)
big = max(el, key=lambda r: r["species"])
out(X, "S11c", "elimination, largest copy finished; first copy not finished", f"{big['species']} species in {big['rref']['seconds']:.0f} s; "
    f"{min((r['species'] for r in env_ if not r.get('rref')), default='none')} species; e_coli_core finished {bool(eco.get('rref'))}", LACR)
sm = times("small")
lad = [m for m in sm if m.startswith("envz") or m.startswith("e_coli")]
out(X, "S11d", "tangent route on the copies and e_coli_core, longest (median of three, 1 core)", f"{max(med(sm[m]) for m in lad):.2f} s", f"{TIMING}/small")
k12 = min(range(len(ff.TAUS)), key=lambda k: abs(ff.TAUS[k] - 1e-12))
out(X, "S11e", f"exact nullspace, routes differ at tolerance {ff.TAUS[k12]:g}", f"{e['gap'][k12] + e['float_differ'][k12]} of {len(e['models'])}", f"{ff.KF}/parts")
out(X, "S11f", "tolerance 0, exact rows: certified, gap", f"{e['at_zero'][0]}, {e['at_zero'][1]}", f"{ff.KF}/parts")
# S12. The exact second-order test
for m in ("iCN718", "iHN637", "iNJ661"):
    d = json.load(open(f"kinetic_variants/second_order_results/{m}_modp.json"))
    c = d["candidates"]
    ma_ = MCa[m]
    acrr_k1 = ratios(ma_["tan"]["acrr_group_sizes"]) - ratios(ma_["it"]["acrr_group_sizes"])
    out(X, "S12a", f"{m}: tested ACR, ACRR pairs, complex pairs (attached); nonzero", f"{c['acr']['n']}, {c['acrr']['n']} (as k - 1 ratios {acrr_k1}), "
        f"{c['complex_pairs']['n']}; nonzero {c['acr']['nonzero'] + c['acrr']['nonzero'] + c['complex_pairs']['nonzero']}",
        "kinetic_variants/second_order_results/, member_counts/<model>.attached.json")
    out(X, "S12b", f"{m}: controls (proven relations nonzero; random row-space vectors nonzero), consistent, ranks N, W, J",
        f"{d['control_proven']['nonzero']} of {d['control_proven']['n']}; {d['control_random_rowspace']['nonzero']} of {d['control_random_rowspace']['n']}; "
        f"{d['second_order_consistent']}; {d['rank_N']}, {d['rank_W']}, {d['rank_J']}", "kinetic_variants/second_order_results/")

# ---- numbers that no result file holds
N = "not here"
out(N, "N1", "COCOA.jl runs on the 4514Y nodes (53 of 183), 125 GB reserved, peak 86 GB; reference 2 cores, 4 to 26 GB, peak 20 GB",
    "-", "SLURM accounting (sacct) of the COCOA.jl and reference jobs, and the sbatch files")
out(N, "N2", "threads of the SVD (2 on most networks, up to 16)", "-", "kinetic_variants/float_four.sbatch, float_exactcmp.sbatch")
out(N, "N3", "software versions", "-", "the environment (python -c 'import numpy, scipy, cobra, flint')")
out(N, "N4", "kinetic stage on iLJ478, COCOA.jl 4.5 s and reference 3.5 min", "-", "constants of figures/make_fig3_methods.py, from "
    "scratch cocoa_kin_timing_2026-11-03 and ref_rstage_timing_2026-11-03")
