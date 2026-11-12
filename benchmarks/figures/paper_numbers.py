"""Every number of the Results and Materials and Methods text, recomputed from the result files, each with its source.

One line per number, in the order of the text: the section, what the number is, the value, and the file or directory it
comes from. After a rerun, point the variables below at the new result directories, run this script and update the
text where a value changed. Numbers that no result file holds (hardware, software versions, SLURM settings) are listed
at the end with where they come from.

    python figures/paper_numbers.py        # from benchmarks/

Result directories (environment variable, default):
    KIMONU_DET_DIR, KIMONU_ATT_DIR    /scratch/sadreyev_lab/pschneider/kimonu_timing_seeds/{det,att}   (figdata.json)
    KIMONU_REF_DIR                    /scratch/sadreyev_lab/pschneider/kimonu_timing_final/ref
    KIMONU_TIMING                     /scratch/sadreyev_lab/pschneider/kimonu_timing_seeds   (tan, small)
    KIMONU_LACR                       /scratch/sadreyev_lab/pschneider/kimonu_local_acr
    KF_OUT                            /scratch/sadreyev_lab/pschneider/kimonu_float_four
    CERT_OUT                          /scratch/sadreyev_lab/pschneider/kimonu_cert_modp
    GAP_OUT                           /scratch/sadreyev_lab/pschneider/kimonu_gaproots
    SZ_DIR                            /scratch/sadreyev_lab/pschneider/kimonu_sz_flux (8 networks; kimonu_sz_flux_33 for all)
    FVA_OUT                           /scratch/sadreyev_lab/pschneider/kimonu_timing_seeds/fva
    UPTAKE_LISTS                      /scratch/sadreyev_lab/pschneider/kimonu_uptake_modes/lists
Repository files: figdata.json, figures/pair_counts.json, figures/pair_sets.json, kinetic_variants/uptake_modes/,
kinetic_variants/member_counts/, kinetic_variants/fb_vs_ma/, kinetic_variants/one_lp_results.log, unblocked/.
"""
import glob, json, os, re, statistics, sys

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(B)
sys.path.insert(0, os.path.join(B, "figures"))
SCR = "/scratch/sadreyev_lab/pschneider"
env = lambda k, d: os.environ.get(k, d)
TIMING = env("KIMONU_TIMING", f"{SCR}/kimonu_timing_seeds")
REF_DIR = env("KIMONU_REF_DIR", f"{SCR}/kimonu_timing_final/ref")
DET_DIR = env("KIMONU_DET_DIR", f"{TIMING}/det")
LACR = env("KIMONU_LACR", f"{SCR}/kimonu_local_acr")
CERT = env("CERT_OUT", f"{SCR}/kimonu_cert_modp")
GAP = env("GAP_OUT", f"{SCR}/kimonu_gaproots")
SZ = env("SZ_DIR", f"{SCR}/kimonu_sz_flux")
FVA = env("FVA_OUT", f"{TIMING}/fva")
MODELS = [l.strip() for l in open("kinetic_variants/models33.txt") if l.strip()]
GAP_ORDER = ["iCN718", "iNJ661", "iHN637"]                     # the order of the text
med = statistics.median


def out(section, what, value, source):
    print(f"{section:8s} | {what} | {value} | {source}")


def pct(a, b):
    return f"{100 * a / b:.1f}%"


def maxrss_gb(paths):
    """largest 'Maximum resident set size' of /usr/bin/time -v logs, in GB (10^6 kB, as SLURM reports memory)"""
    kb = [int(m.group(1)) for p in paths for m in [re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(p).read())] if m]
    return max(kb) / 1e6 if kb else None


# ---- inputs shared by several numbers
D = json.load(open("figdata.json"))
nC = D["n_complexes"]
UM = {}
for f in sorted(glob.glob("kinetic_variants/uptake_modes/*.json")):      # <model>.detached.json repeats the detached block
    d = json.load(open(f))
    if "detached" in d:
        UM.setdefault(d["model"], d)
assert sorted(UM) == sorted(MODELS), len(UM)
det = {m: UM[m]["detached"] for m in MODELS}
it = {m: det[m]["it"] for m in MODELS}
gap = sorted(m for m in MODELS if det[m]["gap"])

# ---- Results: genome-scale networks
R = "Results"
out(R, "networks", len(MODELS), "kinetic_variants/models33.txt")
out(R, "complexes per network, zero complex included (text)", f"{min(UM[m]['n_complexes'] for m in MODELS):,} to {max(UM[m]['n_complexes'] for m in MODELS):,}",
    "kinetic_variants/uptake_modes/ n_complexes")
out(R, "complexes per network without the zero complex (table S2)", f"{min(nC.values()):,} to {max(nC.values()):,}", "figdata.json n_complexes")
out(R, "networks on which the routes agree", f"{len(MODELS) - len(gap)} of {len(MODELS)} (gap: {', '.join(gap)})", "uptake_modes/ detached gap")
v = sorted(x["pairs"] for x in it.values())
out(R, "coupled complex pairs per network, iterated route", f"{v[0]:,} to {v[-1]:,}, median {med(v):,}", "uptake_modes/ detached it.pairs")
v = [x["modules"] for x in it.values()]
out(R, "kinetic modules per network", f"{min(v):,} to {max(v):,}", "uptake_modes/ detached it.modules")
v = [x["giant"] for x in it.values()]
out(R, "complexes in the largest module", f"{min(v):,} to {max(v):,}", "uptake_modes/ detached it.giant")
out(R, "ACR metabolites, all 33", f"{sum(x['acr'] for x in it.values()):,}", "uptake_modes/ detached it.acr")
out(R, "ACRR pairs, all 33", f"{sum(x['acrr_pairs'] for x in it.values()):,}", "uptake_modes/ detached it.acrr_pairs")
out(R, "free metabolites with ACR, networks", f"{sum(x['acr_free'] for x in it.values())} on {sum(1 for x in it.values() if x['acr_free'])}",
    "uptake_modes/ detached it.acr_free")
out(R, "free ACRR pairs", f"{sum(x['acrr_pairs_free'] for x in it.values()):,}", "uptake_modes/ detached it.acrr_pairs_free")

MC = [json.load(open(f)) for f in sorted(glob.glob("kinetic_variants/member_counts/*.json")) if not f.endswith(".attached.json")]
assert len(MC) == 33
ex = [d for d in MC if not d["gap"]]
cx = sorted(100 * d["it"]["complexes_in_modules"] / d["n_complexes"] for d in ex)
fr = sorted(100 * d["it"]["fixed_members_free"] / d["n_free"] for d in ex)
out(R, f"complexes in a kinetic module, {len(ex)} exact networks", f"{cx[0]:.1f} to {cx[-1]:.1f}%, median {med(cx):.1f}%", "kinetic_variants/member_counts/")
out(R, f"free metabolites with ACR or in an ACRR group, {len(ex)} exact networks", f"{fr[0]:.1f} to {fr[-1]:.1f}%, median {med(fr):.1f}%",
    "kinetic_variants/member_counts/")

FM = [json.load(open(f"kinetic_variants/fb_vs_ma/{m}.json")) for m in MODELS]
S = lambda get: sum(get(d) for d in FM)
nr = S(lambda d: d["n_reactions"])
fb, sub = S(lambda d: d["flux_coupling"]["members"]), S(lambda d: d["flux_coupling_and_shared_substrate"]["members"])
ma, ma_t = S(lambda d: d["rates_iterated"]["members"]), S(lambda d: (d["rates_tangent"] or d["rates_iterated"])["members"])
out(R, "reactions fully coupled under flux balance", pct(fb, nr), "kinetic_variants/fb_vs_ma/ flux_coupling.members")
out(R, "reactions in a rate class under mass action", pct(ma, nr), "fb_vs_ma/ rates_iterated.members")
out(R, "of these, points from flux balance, from a shared substrate complex, from kinetic modules",
    f"{100 * fb / nr:.1f}, {100 * (sub - fb) / nr:.1f}, {100 * (ma - sub) / nr:.1f}", "fb_vs_ma/ flux_coupling_and_shared_substrate.members")
out(R, "tangent route, reactions in a rate class", pct(ma_t, nr), "fb_vs_ma/ rates_tangent.members")

import make_data_s1_s2
n, tot = make_data_s1_s2.data_s2(os.devnull)
out(R, "more free ACR metabolites with attached uptakes", tot["lost"], "UPTAKE_LISTS (uptake_lists.py), as Data S2")
out(R, "of these extracellular or periplasmic", f"{tot['extracellular'] + tot['periplasmic']} = {pct(tot['extracellular'] + tot['periplasmic'], tot['lost'])}"
    f" ({tot['extracellular']} e, {tot['periplasmic']} p); in a detached ACRR group {tot['lost_in_acrr']}; gained by detaching {tot['gained']}", "UPTAKE_LISTS")

pairs = [det[m]["tan"]["pairs"] - det[m]["it"]["pairs"] for m in GAP_ORDER]
rel = [det[m]["tan"]["acr_free"] + det[m]["tan"]["acrr_pairs_free"] - det[m]["it"]["acr_free"] - det[m]["it"]["acrr_pairs_free"] for m in GAP_ORDER]
out(R, f"tangent-only coupled complex pairs on {', '.join(GAP_ORDER)}", ", ".join(f"{x:,}" for x in pairs), "uptake_modes/ detached tan - it")
out(R, "tangent-only free ACR and ACRR relations, same networks", ", ".join(map(str, rel)), "uptake_modes/ detached tan - it")
cases = [json.load(open(f)) for f in sorted(glob.glob(f"{GAP}/case_*.json"))]
roots = []
for m in GAP_ORDER:
    ok = [len(c["roots"]) for c in cases if c["model"] == m and c["detached"] and c["roots"] != ["none"] and c["closed"] == c["n_relations"]]
    roots.append(min(ok) if ok else "none closes all")
out(R, "open questions (fewest roots that close every open relation), same networks", ", ".join(map(str, roots)), f"{GAP}/case_*.json")

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
out(R, f"exact nullspace, tolerance {ff.TAUS[k14]:g}: networks on which the routes differ", f"{e['gap'][k14] + e['float_differ'][k14]} of {len(e['models'])}",
    f"{ff.KF}/parts (make_float_four.draw_certificate)")
out(R, "exact nullspace, most networks certified at a tolerance above 0", max(e["certified"]), f"{ff.KF}/parts")
out(R, "SVD nullspace, most networks certified at any tolerance", max(s["certified"] + [s["at_zero"][0]]), f"{ff.KF}/parts")
out(R, "SVD nullspace at tolerance 0: false certificates", f"{s['at_zero'][3]} of {len(s['models'])}", f"{ff.KF}/parts")

C = [json.load(open(f)) for f in sorted(glob.glob(f"{CERT}/*.json"))]
out(R, "genericity bound per network (largest over both routes)", f"{max(max(c['bound_tangent'], c['bound_iterated']) for c in C):.2g} ({len(C)} networks)",
    f"{CERT}/*.json bound_*")

# timing and comparison with the other tools (Fig. 4)
kd = D["kimonu_det"]
out(R, "Kimonu default, largest network", f"{kd['iAF1260b']['seconds']:.1f} s (iAF1260b); slowest of all {max(v['seconds'] for v in kd.values()):.1f} s",
    "figdata.json kimonu_det (median of three seeds)")
sec = lambda v: med(r["seconds"] for r in v) if isinstance(v, list) else v["seconds"]
for arm, name in (("ref_exh", "exhaustive reference"), ("ref_sampled", "sampled reference"), ("cocoa_false", "COCOA.jl efficient=false"),
                  ("cocoa_true", "COCOA.jl efficient=true")):
    su = [sec(v) / kd[m]["seconds"] for m, v in D[arm].items()]
    out(R, f"median speed-up against the {name}", f"{med(su):.2g} over {len(su)} networks", "figdata.json")
tan = {}
for f in glob.glob(f"{TIMING}/tan/tangent_*_rep*.json"):
    r = json.load(open(f)); tan.setdefault(r["network"], []).append(r["seconds"])
m_tan = max(tan, key=lambda m: med(tan[m]))
out(R, "tangent route, longest run (median of three)", f"{med(tan[m_tan]) / 60:.1f} min ({m_tan})", f"{TIMING}/tan")

PS = json.load(open("figures/pair_sets.json"))
keys = ["balanced", "conc", "n_kinetic", "giant", "acr", "acrr", "acr_free", "acrr_free"]
same = sum(D["kimonu_ref"][m][k] == D["ref_exh"][m][k] for m in D["ref_exh"] for k in keys)
out(R, "reference-equivalent mode equals the exhaustive reference", f"{same} of {len(keys) * len(D['ref_exh'])} values on {len(D['ref_exh'])} networks",
    "figdata.json kimonu_ref, ref_exh")
r0 = {m: PS[m]["ref_exh"][0] for m in PS}
out(R, "reference pairs not in Kimonu's default", sum(r["pairs_not_in_kimonu"] for r in r0.values()), "figures/pair_sets.json ref_exh")
miss = {m: len(r["acr_not_in_kimonu"]) + len(r["acrr_not_in_kimonu"]) for m, r in r0.items()}
out(R, "reference ACR and ACRR relations not in Kimonu's default", ", ".join(f"{m} {k}" for m, k in miss.items() if k) or "none", "figures/pair_sets.json")
out(R, "networks of the exhaustive reference on which the routes agree", f"{sum(1 for m in D['ref_exh'] if not det[m]['gap'])} of {len(D['ref_exh'])}",
    "uptake_modes/ detached gap")
SEEDS = [json.load(open(f)) for f in sorted(glob.glob(f"{env('SEED_OUT', SCR + '/kimonu_seeds')}/*.json"))]
out(R, "networks with identical modules, ACR and ACRR under three seeds", f"{sum(all(d['identical'].values()) for d in SEEDS)} of {len(SEEDS)}",
    "SEED_OUT (seed_invariance.py)")


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
out(R, "same check on the timing runs (seeds 0, 1, 2)", f"{same} of {len(MODELS)}", DET_DIR)

# ---- Materials and Methods
M = "Methods"
LP = [json.loads(l) for l in open("kinetic_variants/one_lp_results.log") if l.startswith("{")]
out(M, "one LP finds exactly the blocked reactions of FVA", f"{sum(r['lp_wrong'] == 0 and r['lp_blocked'] == r['fva_blocked'] for r in LP)} of {len(LP)}",
    "kinetic_variants/one_lp_results.log")
out(M, "time of the one LP", f"{min(r['lp_seconds'] for r in LP):.3f} to {max(r['lp_seconds'] for r in LP):.3f} s", "kinetic_variants/one_lp_results.log")
lb, ub = set(), set()
for m in MODELS:
    c = json.load(open(f"unblocked/{m}_fixed.crn.json")); lb |= set(c["lb"]); ub |= set(c["ub"])
out(M, "bounds on every reaction", f"lower {sorted(lb)}, upper {sorted(ub)}", "unblocked/*.crn.json")
w = [float(x) for x in re.findall(r"iAF1260b.*\n.*wall ([\d.]+) s", open(sorted(glob.glob(f"{FVA}/preprocess_*.log"))[-1]).read())]
out(M, "preprocessing iAF1260b, load to factorization", f"{med(w):.1f} s (runs {', '.join(f'{x:.1f}' for x in w)})", f"{FVA}/preprocess_*.log")
out(M, "second evaluation (mod p) reproduces both routes", f"{sum(c['tangent_match'] and c['iterated_match'] for c in C)} of {len(C)}", f"{CERT}/*.json")
Z = [json.load(open(f)) for f in sorted(glob.glob(f"{SZ}/*.json"))]
lp = lambda z: z["lp"] if "lp" in z else z["default"]         # records before 2026-11-11 call the LP arm "default"
ri = [z["uniform"]["iterated"]["seconds"] / lp(z)["iterated"]["seconds"] for z in Z]
rt = [z["uniform"]["tangent"]["seconds"] / lp(z)["tangent"]["seconds"] for z in Z]
out(M, f"uniform 2^64 flux against the LP point, runtime ratio ({len(Z)} networks)",
    f"iterated {min(ri):.1f} to {max(ri):.1f}, tangent {min(rt):.1f} to {max(rt):.1f}; identical results "
    f"{sum(z['iterated_agree'] and z['tangent_agree'] for z in Z)} of {len(Z)}", SZ)
for f in ("kinetic_variants/float_exactcmp.py", "kinetic_variants/float_four.py"):
    r = re.search(r"rnd >= (\d+)", open(f).read())
    out(M, "rounds of the iterated route in the floating-point comparison, at most", r.group(1) if r else "?", f)
for arm in ("cocoa_false", "cocoa_true", "ref_sampled"):
    out(M, f"seeds of {arm}", sorted({len(v) for v in D[arm].values()}), "figdata.json")
for m in ("iSB619", "iNF517"):
    out(M, f"exhaustive reference run time, {m}", f"{D['ref_exh'][m]['seconds'] / 86400:.1f} days", "figdata.json ref_exh")
logs = lambda arms, d: [p for a in arms for p in glob.glob(f"{d}/{a}/*.log")]
out(M, "Kimonu peak memory, iterated route and reference-equivalent mode",
    f"{max(maxrss_gb(logs(['det', 'att'], TIMING)), maxrss_gb(glob.glob(f'{REF_DIR}/*.log'))):.2f} GB", f"{TIMING}/{{det,att}}, {REF_DIR} (/usr/bin/time -v)")
out(M, "Kimonu peak memory, tangent route", f"{maxrss_gb(logs(['tan', 'tanatt'], TIMING)):.2f} GB", f"{TIMING}/{{tan,tanatt}}")
import numpy as np
for arm in ("ref_exh", "ref_sampled", "cocoa_false", "cocoa_true", "kimonu_det", "kimonu_ref"):
    ms = [m for m in D["models"] if m in D[arm]]
    b = np.polyfit(np.log10([nC[m] for m in ms]), np.log10([sec(D[arm][m]) for m in ms]), 1)[0]
    out(M, f"fit exponent of wall time on complexes, {arm}", f"{b:.2f} over {len(ms)} networks", "figdata.json (as make_figure3.py, Fig. 4A)")
out(M, "tangent route, longest run", f"{med(tan[m_tan]) / 60:.1f} min on {m_tan}", f"{TIMING}/tan")
rows = json.load(open(f"{LACR}/small.json")) + [json.load(open(f)) for f in sorted(glob.glob(f"{LACR}/ladder_*.json"))]
copies = sorted(int(r["network"][4:]) for r in rows if re.fullmatch(r"envz\d+", r["network"]))
small_ = [r["network"] for r in rows if not r["network"].startswith(("envz", "e_coli"))]
out(M, "symbolic baselines: small networks, EnvZ-OmpR copies, total", f"{len(small_)}, copies {copies}, total {len(rows)}", f"{LACR}")
for key, name in (("rref", "elimination"), ("minors", "minors")):
    done = [r for r in rows if r.get(key)]
    agree = all(r.get(f"agree_{key}") for r in done)
    out(M, f"symbolic {name} finished, agrees with the tangent route", f"{len(done)} of {len(rows)}, agree {agree}", LACR)
timeouts = sorted(int(r["network"][4:]) for r in rows if re.fullmatch(r"envz\d+", r["network"]) and r.get("minors_timeout") is not None)
out(M, "copies on which the minors reached the limit", timeouts, LACR)
out(M, "reference-equivalent mode, longest run", f"{max(v['seconds'] for v in D['kimonu_ref'].values()):.1f} s", "figdata.json kimonu_ref")

# ---- numbers that no result file holds
N = "not here"
out(N, "COCOA.jl runs on the 4514Y nodes (53 of 183), 125 GB reserved, peak 86 GB; reference 2 cores, 4 to 26 GB, peak 20 GB",
    "-", "SLURM accounting (sacct) of the COCOA.jl and reference jobs, and the sbatch files")
out(N, "Kimonu 1 GB and tangent route 3 GB reserved, 1 core, exclusive node", "-", "kinetic_variants/timing_exclusive.sbatch")
out(N, "threads of the SVD (2 on most networks, up to 16)", "-", "kinetic_variants/float_four.sbatch, float_exactcmp.sbatch")
out(N, "software versions", "-", "the environment (python -c 'import numpy, scipy, cobra, flint')")
