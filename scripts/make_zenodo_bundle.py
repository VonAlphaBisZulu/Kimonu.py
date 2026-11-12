#!/usr/bin/env python3
"""Assemble the Zenodo archive of the paper, the records that its figures, tables, Data S1 and S2 and numbers read.

    python make_zenodo_bundle.py <out dir> <figures dir> [recorded lists ...]      # from scripts/
    bash zenodo_bundle.sh [out dir]          # all steps in one command, on the login node (records, build, rebuild test)

PARTS maps each folder of the archive to its source (a path under scripts/ or on the cluster scratch space), to
where the scripts of scripts/ look for it (a path under scripts/, or an environment variable), to the paper
items that read it and to the script or SLURM job script that produced it. <figures dir> holds the figures as published,
exactly the files of PUBLISHED (Figure1.svg to Figure4.svg, supp/FigS1_stage_timing.svg to supp/FigS6_float_errors.svg,
supp/DataS1_values.csv, supp/DataS2_uptake_acr.csv). The recorded lists are those of figures/record_inputs.py, one per figure or number
script, run from scripts/ of a copy of the checkout, because the scripts rewrite figdata.json and figures/*.json:

    python figures/record_inputs.py lists/<script>.txt figures/<script>.py [args]

for gather_figdata, pair_sets, seed_spread, dataS1_S2, tableS2_S3 and paper_numbers (Python 3.12 environment
of the runs) and fig4_figS1_S2_benchmark, figS3_S4_uptakes_seeds, figS5_S6_float, fig3AB_hierarchy and fig3CD_coverage (matplotlib
3.8), as zenodo_bundle.sh does. Every file a script read must be in the archive (the published SVGs of <figures dir>
included), or in SNAPSHOT_ONLY, or the script stops. Code is not copied into the data folders. The source snapshot (git archive of $KIMONU_REV, default HEAD) holds it. ZENODO_DOI,
ARTICLE_DOI, ARCHIVE_VERSION (the version of this record) and RELEASE_DATE fill the placeholders of README.md and CITATION.cff when
set (for a rebuild once the DOI is reserved). Every archive path and every repository script that README.md names
in backticks must exist, or the script stops. Writes into <out dir>

    bundle/   README.md (ZENODO_README.md), LICENSE (ZENODO_LICENSE), CITATION.cff, link_records.sh,
              Kimonu.py-<commit>.zip, networks/, results/, figures/, SHA256SUMS (every other file)
    upload/   Kimonu_records.zip (bundle/ as one deflate zip, top folder Kimonu_records/) and README.md
    zenodo_metadata.json, description.html, UPLOAD.md (ZENODO_UPLOAD.md), when $MANUSCRIPT names the manuscript
"""
import datetime, fnmatch, glob, hashlib, json, os, re, shutil, subprocess, sys, zipfile

SCRATCH = "/scratch/sadreyev_lab/pschneider/"
RERUN = SCRATCH + "kimonu_rerun_1111/"
CODE = (".py", ".pyc", ".sbatch", ".sh", ".R", ".m", ".jl")
RUN_COMMIT = "451a748"          # Kimonu commit of the final runs (provenance_*.txt); later src gives the same results
PLACED, ENV = "place", "env"

# archive path, source (a directory, taken whole without code, or a glob), (PLACED, path under scripts/) or
# (ENV, variable), the paper items it feeds, content, and what produced it (paths under scripts/, SLURM job ids)
TE = "scripts/kinetic_variants/timing_exclusive.sbatch"
PARTS = [
    ("networks/unblocked", "unblocked/*_fixed.crn.json", (PLACED, "unblocked"), "every tool, Fig. 3C, D, Tables S2, S3",
     "the 33 networks with blocked reactions removed, the input of every tool (crn.json, triplets of A, Y and S)",
     "`scripts/write_unblocked_crn.py`, run by `scripts/make_unblocked.sbatch` (FVA with SCIP), from `networks/with_blocked_reactions/`"),
    ("networks/with_blocked_reactions", "sweep/*_fixed.crn.json", (PLACED, "sweep"), "Materials and Methods",
     "the same networks before blocked-reaction removal, converted from the elementary-step models of the reference",
     "`scripts/export_split_model.py`, from the elementary-step models (.mat) of Langary et al. (2025)"),
    ("networks/reference_coupling_matrices", "coupling_exact/*.mat", (PLACED, "coupling_exact"), "input of the reference",
     "exact coupling matrices that the reference's R stage read", "`scripts/matlab_reference/write_coupling_exact.py`"),
    ("networks/deficiency", "deficiency/deficiency_33.json", (PLACED, "deficiency"), "Table S2, Supplementary Text S10",
     "deficiency, linkage classes and ranks of the 33 networks", "`scripts/deficiency/deficiency_split.py`"),
    ("results/kimonu/timing", RERUN + "timing", (ENV, "KIMONU_TIMING"), "Fig. 4A to D, F, figs. S1 to S4, Data S1, Materials and Methods",
     "timed runs, each one process on one core of an exclusive node, replicates rep1 to rep3 at seeds 0, 1, 2 of r*. "
     "Arms det and att, iterated route (the default) with detached and attached uptakes. Arm ref, reference-equivalent "
     "mode (no r*, one computation repeated). Arms tan and tanatt, tangent route. Arm small, tangent route on the small "
     "networks and the EnvZ-OmpR ladder. Each run has its JSON (modules, ACR, ACRR, seconds) and /usr/bin/time log, "
     "each arm provenance_*.txt and hardware_*.txt, and slurm_jobs.txt lists the SLURM jobs",
     f"`{TE}`, which runs `scripts/run_kimonu_crn.py` (det, att), `scripts/run_kimonu_refmode.py` (ref) and "
     "`scripts/kinetic_variants/tangent_cell.py` (tan, tanatt, small), SLURM jobs 572311 (det, att, ref, small), 572312 (tan), 572313 (tanatt)"),
    ("results/kimonu/uptake_modes", "kinetic_variants/uptake_modes", (PLACED, "kinetic_variants/uptake_modes"),
     "Fig. 3A, B (through hierarchy_detached), Data S1, Results", "both routes with detached and attached uptakes, per network",
     "`scripts/kinetic_variants/uptake_modes.py`, run by `scripts/kinetic_variants/uptake_modes.sbatch`, SLURM job 572113"),
    ("results/kimonu/hierarchy_detached", "kinetic_variants/hierarchy_detached", (PLACED, "kinetic_variants/hierarchy_detached"),
     "Fig. 3A, B", "per network, the counts of both routes that Fig. 3A, B shows (detached uptakes)",
     "`scripts/kinetic_variants/hierarchy_detached.py`, from `results/kimonu/uptake_modes/` and `results/kimonu/member_counts/`"),
    ("results/kimonu/member_counts", "kinetic_variants/member_counts", (PLACED, "kinetic_variants/member_counts"),
     "Fig. 3C, D, Tables S2, S3, Data S1, Results", "complexes in kinetic modules, metabolites with ACR or in an ACRR group, group sizes",
     "`scripts/kinetic_variants/member_counts.py`, run by `scripts/kinetic_variants/member_counts.sbatch`, SLURM jobs 572115 (detached), 572116 (attached)"),
    ("results/kimonu/flux_balance_vs_mass_action", "kinetic_variants/fb_vs_ma", (PLACED, "kinetic_variants/fb_vs_ma"), "Results",
     "reactions coupled under flux balance against rate classes under mass action",
     "`scripts/kinetic_variants/fb_vs_ma.py`, run by `scripts/kinetic_variants/fb_vs_ma.sbatch` (SLURM job 572117), summary by `scripts/kinetic_variants/fb_vs_ma_summary.py`"),
    ("results/kimonu/acr_lists_by_uptake_mode", RERUN + "uptake_lists", (ENV, "UPTAKE_LISTS"), "Data S2, Results",
     "ACR and ACRR lists of the iterated route per uptake mode, metabolites as indices into the network's own list",
     "`scripts/kinetic_variants/uptake_lists.py`, run by `scripts/kinetic_variants/uptake_lists.sbatch`, SLURM job 572114"),
    ("results/kimonu/gap_roots", RERUN + "gap_roots", (ENV, "GAP_OUT"), "Results, Supplementary Text S6",
     "the open questions the three gap networks come down to",
     "`scripts/kinetic_variants/gap_roots.py`, run by `scripts/kinetic_variants/gap_roots.sbatch` with `scripts/kinetic_variants/gap_roots_cases.txt`, SLURM job 572118"),
    ("results/kimonu/second_order", "kinetic_variants/second_order_results", (PLACED, "kinetic_variants/second_order_results"),
     "Supplementary Text S12", "exact second-order test modulo a prime on the gap networks",
     "`scripts/kinetic_variants/second_order_modp.py`, run by `scripts/kinetic_variants/second_order_modp.sbatch`, SLURM job 572119"),
    ("results/kimonu/gap_experiments/injectivity", "kinetic_variants/injectivity_results", (PLACED, "kinetic_variants/injectivity_results"),
     "Supplementary Text S6, S7", "sign condition, injectivity and tangent dimensions on the gap networks",
     "`scripts/kinetic_variants/injectivity_gap.py`, run by `scripts/kinetic_variants/injectivity.sbatch`"),
    ("results/kimonu/gap_experiments/descartes", "kinetic_variants/descartes_*.json", (PLACED, "kinetic_variants"),
     "Supplementary Text S6", "Descartes step on the gap networks",
     "`scripts/kinetic_variants/descartes_vs_tangent.py`, run by `scripts/kinetic_variants/descartes_array.sbatch`"),
    ("results/kimonu/gap_experiments/autonomous_blocks", "kinetic_variants/autonomous_blocks_33.jsonl", (PLACED, "kinetic_variants"),
     "Supplementary Text S6", "relations proven by autonomous blocks",
     "`scripts/kinetic_variants/autonomous_blocks.py`, run by `scripts/kinetic_variants/autonomous_blocks.sbatch` (one line per network)"),
    ("results/kimonu/uniform_flux", RERUN + "sz_flux_33", (ENV, "SZ_DIR"), "Supplementary Text S8, Materials and Methods",
     "both routes at a uniformly drawn r* (64 bits) against the LP point, 33 networks",
     "`scripts/kinetic_variants/sz_flux.py`, run by `scripts/kinetic_variants/sz_flux33.sbatch`, SLURM job 572120"),
    ("results/kimonu/genericity_bound", RERUN + "genericity_bound", (ENV, "GB_OUT"), "Results, Supplementary Text S8",
     "Schwartz-Zippel bound of the uniform r* per network and route (summary.json)",
     "`scripts/kinetic_variants/genericity_bound.py`, from `results/kimonu/uniform_flux/`, SLURM job 572316"),
    ("results/kimonu/lp_point_rounding", RERUN + "scale_check", (ENV, "SCALE_DIR"), "Supplementary Text S8",
     "rounding of the LP point (scales 100, 10^4, 10^5, 10^6, unrounded) on iNJ661 and iHN637, and the same timed on "
     "an exclusive node (exclusive/, three replicates per network)",
     "`scripts/kinetic_variants/scale_check.py` (SLURM jobs 572317, iNJ661, and 572318, iHN637), `scripts/kinetic_variants/scale_check_1e4.py` "
     "(iHN637_1e4.json, job 572319) and `scripts/kinetic_variants/scale_check_exclusive.sbatch` (exclusive/, job 572320)"),
    ("results/kimonu/float_comparison/parts", SCRATCH + "kimonu_float_four/parts", (ENV, "KF_OUT"),
     "Fig. 3E, F, figs. S5, S6, Results", "floating-point comparison, per network, nullspace and tolerance",
     "`scripts/kinetic_variants/float_exactcmp.py` (the parts *_exactcmp*), run by `scripts/kinetic_variants/float_exactcmp.sbatch`, and "
     "`scripts/kinetic_variants/float_four.py` (the other parts), run by `scripts/kinetic_variants/float_four.sbatch` with `scripts/kinetic_variants/float_four_parts.txt`"),
    ("results/kimonu/symbolic_local_acr", SCRATCH + "kimonu_local_acr", (ENV, "KIMONU_LACR"),
     "Fig. 4B, Supplementary Text S11", "symbolic local-ACR criteria (minors, elimination) on the small networks and the ladder",
     "`scripts/kinetic_variants/local_acr_symbolic.py`, run by `scripts/kinetic_variants/local_acr_ladder.sbatch`"),
    ("results/kimonu/preprocessing_time", SCRATCH + "kimonu_timing_seeds/fva/preprocess_*.log", (ENV, "FVA_OUT"),
     "Materials and Methods", "preprocessing time of iAF1260b, load to factorization",
     "`scripts/kinetic_variants/preprocess_time.py`, run by `scripts/kinetic_variants/preprocess_time.sbatch`"),
    ("results/kimonu/one_lp_vs_fva", "kinetic_variants/one_lp_results.log", (PLACED, "kinetic_variants"), "Materials and Methods",
     "one LP against FVA, blocked reactions and times", "`scripts/kinetic_variants/one_lp.py`, run by `scripts/kinetic_variants/one_lp.sbatch`"),
    ("results/cocoa_jl/efficient_false", "sweep/cocoa_v110", (PLACED, "sweep/cocoa_v110"), "Fig. 4A, C, D, figs. S1, S4, Data S1",
     "COCOA.jl v1.1.0, efficient=false, three seeds (cell logs and module files)",
     "`scripts/run_cocoa_acr.jl` with the Julia project `scripts/julia_v110/`, run by `scripts/cocoa_v110.sbatch`, `scripts/cocoa_v110_seeds.sbatch`, "
     "`scripts/cocoa_v110_backfill.sbatch`, `scripts/cocoa_gap_false.sbatch` and `scripts/cocoa_iaf_false.sbatch` (cell lists `scripts/*_cells.txt`)"),
    ("results/cocoa_jl/efficient_true", "sweep/cocoa_eff", (PLACED, "sweep/cocoa_eff"), "Fig. 4A, C, D, F, figs. S1, S4, Data S1",
     "COCOA.jl v1.1.0, efficient=true, three seeds",
     "`scripts/run_cocoa_acr.jl` with `scripts/julia_v110/`, run by `scripts/cocoa_eff.sbatch`, `scripts/cocoa_eff_big.sbatch`, `scripts/cocoa_eff_all.sbatch` "
     "and `scripts/cocoa_gap_true.sbatch`"),
    ("results/cocoa_jl/kinetic_stage_iLJ478", SCRATCH + "cocoa_kin_timing_2026-11-03", None, "Fig. 4F",
     "COCOA.jl's kinetic stage on iLJ478 from saved concordance labels, three replicates (the 4.5 s of Fig. 4F)",
     "`run.sh` and `kin_only.jl` in this folder (SLURM job 568621), not in the repository"),
    ("results/reference/exhaustive", "reference_final", (PLACED, "reference_final"), "Fig. 4A, C, D, F, figs. S1, S2, Data S1",
     "the reference (Langary et al.), exhaustive, with concordance results, R-stage logs and results, retested partitions",
     "`scripts/matlab_reference/run_reference.m` under Octave, run by `scripts/ref_final.sbatch` (`scripts/final_models.txt`); R stage "
     "`scripts/reference_final/acr/run_chain.sh` (`scripts/reference_final/rerun_chain.sbatch`); retests by `scripts/ref_finish.sbatch` with "
     "`scripts/matlab_reference/retest_pairs.py` and `scripts/matlab_reference/retest_rejected.m`"),
    ("results/reference/sampled", "reference_final_sampled", (PLACED, "reference_final_sampled"), "Fig. 4A, C, D, fig. S4, Data S1",
     "the reference, sampled, three seeds",
     "`scripts/matlab_reference/run_reference.m`, run by `scripts/ref_final_sampled.sbatch` (`scripts/final_sampled_cells.txt`)"),
    ("results/reference/kinetic_stage_iLJ478", SCRATCH + "ref_rstage_timing_2026-11-03", None, "Fig. 4F",
     "the reference's R stages on iLJ478, three replicates (the 3.5 min of Fig. 4F)",
     "`run.sh` in this folder (SLURM job 568620), not in the repository"),
    ("results/hardware", "hardware_record.json", (PLACED, "."), "Materials and Methods",
     "nodes, CPUs, memory reserved and used, threads and software versions of the benchmark runs",
     "`scripts/figures/hardware_record.py`, from sacct, the job logs and the result files"),
    ("figures/source_data", "figdata.json", (PLACED, "."), "Fig. 4A, B, F, figs. S1 to S3, Tables S2, S3, Data S1",
     "the values the figures show, per network, tool and mode",
     "`scripts/figures/gather_figdata.py`, from `results/kimonu/timing/` (det, att, ref), `results/cocoa_jl/`, `results/reference/`"),
    ("figures/source_data", "figures/pair_counts.json", (PLACED, "figures"), "Fig. 4C, D, Data S1", "relation counts per tool",
     "`scripts/figures/pair_sets.py`, from `results/kimonu/timing/ref` and the module files of `results/cocoa_jl/` and `results/reference/`"),
    ("figures/source_data", "figures/pair_sets.json", (PLACED, "figures"), "Fig. 4C, D", "set containment against Kimonu",
     "`scripts/figures/pair_sets.py`, as above"),
    ("figures/source_data", "figures/seed_spread.json", (PLACED, "figures"), "fig. S4", "relations that change with the seed",
     "`scripts/figures/seed_spread.py`, from `results/kimonu/timing/det`, `results/cocoa_jl/`, `results/reference/sampled/`"),
    ("figures/source_data", "figures/reference_rstage_end.json", (PLACED, "figures"), "fig. S1",
     "the end of each R stage of the exhaustive reference, the modification time of its ACR log (a zip keeps file times "
     "to 2 s only)", "`scripts/figures/reference_rstage_end.py`, from `results/reference/exhaustive/acr/*.acr.log`"),
]
# the published figures, Data S1 and S2 as <figures dir> must hold them, each with the script that draws or writes it
PUBLISHED = {
    "Figure1.svg": ("Fig. 1", "Illustrator master; `scripts/docs/examples/network/fig1CDE_matrices.py`, then the `fig1_*.py` scripts there"),
    "Figure2.svg": ("Fig. 2", "Illustrator master; `scripts/docs/examples/network/fig2_kinase.py`, then the `fig2*_update_*.py` scripts there"),
    "Figure3.svg": ("Fig. 3", "Illustrator master; `scripts/figures/fig3AB_hierarchy.py`, `scripts/figures/fig3CD_coverage.py`, `scripts/figures/fig3EF_float.py`, then the `fig3*_update_*.py` scripts in `scripts/figures/`"),
    "Figure4.svg": ("Fig. 4", "Illustrator master; `scripts/figures/fig4_figS1_S2_benchmark.py` with `scripts/figures/fig4EF_methods_iLJ478.py`, then `scripts/figures/fig4ABF_update.py`, `scripts/figures/fig3_fig4_update_terms.py`"),
    "supp/FigS1_stage_timing.svg": ("fig. S1", "`scripts/figures/fig4_figS1_S2_benchmark.py`"),
    "supp/FigS2_reference_agreement.svg": ("fig. S2", "`scripts/figures/fig4_figS1_S2_benchmark.py`"),
    "supp/FigS3_uptakes.svg": ("fig. S3", "`scripts/figures/figS3_S4_uptakes_seeds.py`"),
    "supp/FigS4_seeds.svg": ("fig. S4", "`scripts/figures/figS3_S4_uptakes_seeds.py`"),
    "supp/FigS5_certificate.svg": ("fig. S5", "`scripts/figures/figS5_S6_float.py`"),
    "supp/FigS6_float_errors.svg": ("fig. S6", "`scripts/figures/figS5_S6_float.py`"),
    "supp/DataS1_values.csv": ("Data S1", "`scripts/figures/dataS1_S2.py`"),
    "supp/DataS2_uptake_acr.csv": ("Data S2", "`scripts/figures/dataS1_S2.py`"),
}
KEEP_CODE = ("results/cocoa_jl/kinetic_stage_iLJ478", "results/reference/kinetic_stage_iLJ478")   # wrappers not in the repository
SKIP = ("cocoa_kin_timing_2026-11-03/kin_debug.log",            # 23 MB of debug output, no number comes from it
        "Code/kinetic_modules/Get_excel_tables_kinetic_modules.R", "Code/kinetic_modules/code_kineticModule_analysis.R")  # the reference's own code
# files the scripts read that stay in the source snapshot only (tracked in git), with the reason
SNAPSHOT_ONLY = {
    "kinetic_variants/models33.txt": "the list of the 33 networks (array index)",
    "kinetic_variants/s210_*.json": "superseded upstream-set arm of September, read by gather_figdata.py, drawn by no figure",
    "sweep_shared/*": "same arm",
    "sweep_tentative/*": "same arm",
}


def files_of(src):
    """{source path: path relative to the part} for a directory (recursive, without code) or a glob"""
    if os.path.isdir(src):
        return {os.path.join(d, f): os.path.relpath(os.path.join(d, f), src) for d, _, fs in os.walk(src) for f in fs}
    return {p: os.path.basename(p) for p in glob.glob(src) if os.path.isfile(p)}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


out, figdir, listings = sys.argv[1], sys.argv[2], sys.argv[3:]
HERE = os.getcwd()
assert os.path.basename(HERE) == "scripts", "run from scripts/"
git = lambda *a: subprocess.run(["git", *a], capture_output=True, text=True, check=True).stdout.strip()
REV = os.environ.get("KIMONU_REV", "HEAD")
commit, commit_full = git("rev-parse", "--short=7", REV), git("rev-parse", REV)
meta = {k: os.environ.get(k, d) for k, d in dict(ZENODO_DOI="10.5281/zenodo.XXXXXXX", ARTICLE_DOI="10.1126/sciadv.XXXXXXX",
        ARCHIVE_VERSION="X.Y", RELEASE_DATE="YYYY-MM-DD").items()}
bundle, upload = os.path.join(out, "bundle"), os.path.join(out, "upload")
assert not os.path.exists(bundle), f"{bundle} exists"

# ---- collect: archive path -> source, and source -> archive path for the check
plan, where = {}, {}
for dest, src, look, items, what, made in PARTS:
    hit = files_of(src)
    assert hit, f"no files for {dest} ({src})"
    for p, rel in hit.items():
        if p.endswith(CODE) and dest not in KEEP_CODE or p.endswith(SKIP):
            continue
        a = f"{dest}/{rel}"
        assert a not in plan, f"two sources for {a}"
        plan[a] = p
        where[os.path.relpath(p, HERE) if not p.startswith(SCRATCH) else p] = a
have = sorted(os.path.relpath(f, figdir) for f in glob.glob(os.path.join(figdir, "*.svg")) + glob.glob(os.path.join(figdir, "supp", "*"))
              if f.endswith((".svg", ".csv")))
assert have == sorted(PUBLISHED), f"{figdir} must hold exactly {sorted(PUBLISHED)}, not {have}"
for f in have:
    plan[f"figures/{f}"] = os.path.join(figdir, f)

# tracked sources must equal REV, which the snapshot holds
paths = sorted({src for _, src, *_ in PARTS if not src.startswith("/")})
dirty = git("diff", "--name-only", REV, "--", *paths) + git("ls-files", "--others", "--exclude-standard", "--", *paths)
assert not dirty, f"uncommitted changes in the archived files:\n{dirty}"

# ---- check the recorded lists
recorded = {l.strip() for f in listings for l in open(f) if l.strip()}
FIG = os.path.abspath(figdir) + os.sep                          # the in-place figure scripts read the published SVGs
missing = sorted(p for p in recorded if p not in where and not p.startswith(FIG)
                 and not any(fnmatch.fnmatch(p, k) for k in SNAPSHOT_ONLY))
snap = sorted(p for p in recorded if p not in where and p not in missing and not p.startswith(FIG))
assert not missing, "read by a script but not in the archive:\n" + "\n".join(missing[:50])
unread = {}                                                    # archived files no recorded script read, per part
for p, a in where.items():
    if p not in recorded:
        part = next(d for d, *_ in PARTS if a.startswith(d + "/"))
        unread[part] = unread.get(part, 0) + 1

# ---- copy and write
for a, p in sorted(plan.items()):
    d = os.path.join(bundle, a)
    os.makedirs(os.path.dirname(d), exist_ok=True)
    shutil.copy2(p, d)
snapshot = f"Kimonu.py-{commit}.zip"
subprocess.run(["git", "archive", "--format=zip", f"--prefix=Kimonu.py-{commit}/", "-o", os.path.join(bundle, snapshot),
                REV], cwd=git("rev-parse", "--show-toplevel"), check=True)     # from the top, the whole repository


# the variables link_records.sh sets: archive folder, default of the scripts (the paper's runs on the cluster), readers
envs, default = {}, {}
for dest, src, look, *_ in PARTS:
    if look is not None and look[0] == ENV:
        envs[look[1]] = dest[:-len("/parts")] if dest.endswith("/parts") else dest
        default[look[1]] = os.path.dirname(src) if "*" in src or dest.endswith("/parts") else src
for k, arm in (("KIMONU_DET_DIR", "det"), ("KIMONU_ATT_DIR", "att"), ("KIMONU_REF_DIR", "ref")):
    envs[k], default[k] = envs["KIMONU_TIMING"] + "/" + arm, default["KIMONU_TIMING"] + "/" + arm


def readers(var):
    hits = subprocess.run(["git", "grep", "-l", "-E", f"(environ\\.get\\(|env\\()[\"']{var}[\"']", REV, "--", "figures/*.py"],
                          capture_output=True, text=True).stdout.split()
    return ", ".join(f"`{os.path.basename(h.split(':', 1)[1])}`" for h in hits) or "none"



def loc(src, look):
    """where the scripts of scripts/ read a part"""
    if look is None:
        return "no script reads it (Fig. 4F quotes the values as constants of `scripts/figures/fig4EF_methods_iLJ478.py`)"
    if look[0] == ENV:
        return f"variable `{look[1]}`, which `link_records.sh` sets to `{envs[look[1]]}/`"
    return f"`scripts/{look[1]}/`" if os.path.isdir(src) else f"`scripts/{os.path.normpath(os.path.join(look[1], os.path.basename(src)))}`"


def tally(pred):
    hit = [a for a in plan if pred(a)]
    return f"{len(hit):,} ({sum(os.path.getsize(os.path.join(bundle, a)) for a in hit) / 1e6:.2f} MB)"


layout = "| archive folder | content | feeds | produced by | files |\n|---|---|---|---|---|\n"
seen = set()
for dest, src, look, items, what, made in PARTS:
    if dest == "figures/source_data":
        what = f"`{os.path.basename(src)}`, {what}"
    layout += f"| `{dest}/` | {what} | {items} | {made} | {tally(lambda a: a.startswith(dest + '/')) if dest not in seen else '(same folder)'} |\n"
    seen.add(dest)
layout += "".join(f"| `figures/{f}` | the published file | {item} | {made} | 1 |\n" for f, (item, made) in PUBLISHED.items())
reading = "| archive folder or file | where the scripts of `scripts/` read it |\n|---|---|\n" + "\n".join(
    f"| `{dest}/{os.path.basename(src) if dest == 'figures/source_data' else ''}` | {loc(src, look)} |" for dest, src, look, *_ in PARTS)
variables = "| variable | `link_records.sh` sets it to | default in the scripts (the run of the paper on the cluster) | read by (`scripts/figures/`) |\n|---|---|---|---|\n" + "\n".join(
    f"| `{k}` | `{v}/`{' (the scripts read its `parts/`)' if k == 'KF_OUT' else ''} | {default[k]} | {readers(k)} |" for k, v in sorted(envs.items()))
fill = {"@COMMIT@": commit, "@COMMIT_FULL@": commit_full, "@RUN_COMMIT@": RUN_COMMIT, "@DATE@": datetime.date.today().isoformat(),
        "@LAYOUT@": layout.rstrip(), "@READING@": reading, "@VARIABLES@": variables, "@SNAPSHOT@": snapshot, "@YEAR@": meta["RELEASE_DATE"][:4]}
fill.update({f"@{k}@": v for k, v in meta.items()})
PLACEHOLDER = ("XXXXXXX", "X.Y", "YYYY")


def render(src, dst, t=None):
    t = open(src, encoding="utf-8").read() if t is None else t
    for k, v in fill.items():
        t = t.replace(k, v)
    assert not re.search(r"@[A-Z_]+@", t), f"placeholder left in {dst}"
    open(dst, "w", encoding="utf-8").write(t)


render("ZENODO_README.md", os.path.join(bundle, "README.md"))


def named_paths_exist(md):
    """every backticked path of the README must exist: in the archive, in the repository at REV, or (a bare name or
    pattern) as the name of a file in either"""
    tracked = git("ls-tree", "-r", "--name-only", REV, "--full-tree").split("\n")
    names = {os.path.basename(f) for f in tracked} | {n for _, ds, fs in os.walk(bundle) for n in ds + fs}
    targets = {os.path.normpath(f"scripts/{look[1]}") for _, _, look, *_ in PARTS if look and look[0] == PLACED}   # link_records.sh makes them
    bad = []
    for tok in sorted(set(re.findall(r"`([^`\s$=()]+)`", md))):
        if "/" not in tok and not re.search(r"\.(py|sbatch|sh|jl|m|R|txt|json|svg|csv|log|zip)$", tok):
            continue                                                       # a word or value, not a path
        parts = re.match(r"(.*)\{([^}]*)\}(.*)", tok)
        for e in ([parts.group(1) + x + parts.group(3) for x in parts.group(2).split(",")] if parts else [tok]):
            e = e.rstrip("/")
            ok = (glob.glob(os.path.join(bundle, e)) or e in targets or any(fnmatch.fnmatch(f, e) or f.startswith(e + "/") for f in tracked)
                  or "/" not in e and any(fnmatch.fnmatch(n, e) for n in names))
            if not ok:
                bad.append(e)
    assert not bad, "README.md names paths that exist neither in the archive nor in the repository: " + ", ".join(bad)


render("ZENODO_LICENSE", os.path.join(bundle, "LICENSE"))
render(None, os.path.join(bundle, "CITATION.cff"), f"""cff-version: 1.2.0
message: "If you use these records, please cite the article and this archive."
type: dataset
title: 'Run records and figure data for "Exact nullspaces find and certify kinetic modules and other mass-action invariants in metabolic networks"'
authors:
  - family-names: Schneider
    given-names: Philipp C.
    affiliation: "Department of Molecular Biology, Massachusetts General Hospital, Boston, MA, USA"
    orcid: "https://orcid.org/0000-0001-8858-428X"
version: "@ARCHIVE_VERSION@"
date-released: "@RELEASE_DATE@"
doi: "@ZENODO_DOI@"
license: CC-BY-4.0
repository-code: "https://github.com/VonAlphaBisZulu/Kimonu.py"
commit: "{commit_full}"
keywords:
  - kinetic modules
  - absolute concentration robustness
  - mass-action kinetics
  - chemical reaction network theory
  - constraint-based modeling
preferred-citation:
  type: article
  title: "Exact nullspaces find and certify kinetic modules and other mass-action invariants in metabolic networks"
  authors:
    - family-names: Schneider
      given-names: Philipp C.
      orcid: "https://orcid.org/0000-0001-8858-428X"
  journal: "Science Advances"
  year: @YEAR@
  doi: "@ARTICLE_DOI@"
""")
# link_records.sh: links the archive into a checkout, where the scripts of scripts/ look for each part
link = ["#!/usr/bin/env bash",
        "# Links this archive into a checkout of Kimonu.py (the snapshot Kimonu.py-*.zip, unpacked, or a clone at the same",
        "# commit), where the scripts of scripts/ look for each part, and prints the environment variables for the rest.",
        "#   eval \"$(bash link_records.sh <checkout>)\"      then run the scripts from <checkout>/scripts",
        "# Files the checkout holds already (tracked in git) are kept and compared with the archive copy.",
        "set -eu", 'A=$(cd "$(dirname "$0")" && pwd); B=$(cd "$1" && pwd)/scripts; test -d "$B"',
        'one() { if [ -e "$2" ]; then cmp -s "$1" "$2" || echo "differs from the archive: $2" >&2; '
        'else mkdir -p "$(dirname "$2")"; ln -s "$1" "$2"; fi; }',
        'put() { if [ -d "$A/$1" ]; then (cd "$A/$1" && find . -type f | sed "s|^\\./||") | while read -r f; do one "$A/$1/$f" "$B/$2/$f"; done; '
        'else one "$A/$1" "$B/$2"; fi; }']
for dest, src, look, *_ in PARTS:
    if look is not None and look[0] == PLACED:
        if os.path.isdir(src):
            link.append(f"put {dest} {look[1]}")
        else:
            for p, rel in sorted(files_of(src).items()):
                link.append(f"put {dest}/{rel} {'' if look[1] == '.' else look[1] + '/'}{rel}")
link.append("cat <<EOF")
link += [f"export {k}=$A/{v}" for k, v in sorted(envs.items())]
link += ["EOF", 'echo "linked into $B" >&2', ""]
open(os.path.join(bundle, "link_records.sh"), "w").write("\n".join(link))
os.chmod(os.path.join(bundle, "link_records.sh"), 0o755)
named_paths_exist(open(os.path.join(bundle, "README.md"), encoding="utf-8").read())

every = sorted(os.path.relpath(os.path.join(d, f), bundle) for d, _, fs in os.walk(bundle) for f in fs)
with open(os.path.join(bundle, "SHA256SUMS"), "w") as f:
    for a in every:
        f.write(f"{sha256(os.path.join(bundle, a))}  {a}\n")

# ---- zip, then reopen and verify every member against SHA256SUMS
os.makedirs(upload, exist_ok=True)
z = os.path.join(upload, "Kimonu_records.zip")
with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
    for a in every + ["SHA256SUMS"]:
        zf.write(os.path.join(bundle, a), f"Kimonu_records/{a}")
shutil.copy2(os.path.join(bundle, "README.md"), os.path.join(upload, "README.md"))
sums = dict(reversed(l.rstrip("\n").split("  ", 1)) for l in open(os.path.join(bundle, "SHA256SUMS")))
with zipfile.ZipFile(z) as zf:
    assert zf.testzip() is None
    names = [n[len("Kimonu_records/"):] for n in zf.namelist()]
    assert sorted(names) == sorted(every + ["SHA256SUMS"])
    for a in every:
        assert hashlib.sha256(zf.read(f"Kimonu_records/{a}")).hexdigest() == sums[a], a
print(f"commit {commit}; {len(every) + 1} files in bundle/, {sum(os.path.getsize(os.path.join(bundle, a)) for a in every) / 1e6:.0f} MB; "
      f"{z} {os.path.getsize(z) / 1e6:.0f} MB, verified")
print(f"recorded inputs {len(recorded)}: in the archive {len(recorded) - len(snap)}, in the snapshot only {len(snap)}")
print("lp_point_rounding:", ", ".join(a[len("results/kimonu/lp_point_rounding/"):] for a in sorted(plan) if a.startswith("results/kimonu/lp_point_rounding/")))
print("archived files that no recorded script read, per part:", ", ".join(f"{d} {n}" for d, n in unread.items()))
left = [f"{k}={v}" for k, v in meta.items() if any(s in v for s in PLACEHOLDER)]
print("placeholders left in README.md and CITATION.cff:", ", ".join(left) or "none")

# ---- Zenodo metadata (API format), the description and the upload steps, next to upload/
MS = os.environ.get("MANUSCRIPT")                              # the manuscript (Markdown), for title and abstract
if not MS:
    print("MANUSCRIPT not set, so no zenodo_metadata.json, description.html, UPLOAD.md")
    sys.exit(0)
ms = open(MS, encoding="utf-8").read()
abstract = ms[ms.index("# Abstract"):].split("\n")[2].strip()
paper = re.search(r"\*\*Title:\*\* (.+)", ms).group(1).strip()
assert len(abstract.split()) <= 150 and not re.search(r"~~|\*\*|==", abstract + paper), "title or abstract with markup, or too long"
REPO = "https://github.com/VonAlphaBisZulu/Kimonu.py"
TITLE = f'Run records and figure data for "{paper}"'
KEYWORDS = ["kinetic modules", "absolute concentration robustness", "absolute concentration ratio robustness",
            "mass-action kinetics", "chemical reaction network theory", "constraint-based metabolic modeling",
            "genome-scale metabolic networks", "concordant complexes", "flux coupling", "exact arithmetic", "nullspace",
            "benchmark"]
REFERENCES = [
    "Langary, D., Küken, A., Nikoloski, Z. (2025). Kinetic modules are sources of concentration robustness in biochemical "
    "networks. Science Advances 11, eads7269. https://doi.org/10.1126/sciadv.ads7269",
    "Schaffranke, A., Küken, A., Nikoloski, Z. (2026). COCOA.jl, a Julia package for high-performance analysis of "
    "concordance and kinetic modules in biochemical networks. Bioinformatics 42, btag724. "
    "https://doi.org/10.1093/bioinformatics/btag724",
    "Schneider, P., Bekiaris, P. S., von Kamp, A., Klamt, S. (2022). StrainDesign. Bioinformatics 38, 4981-4983. "
    "https://doi.org/10.1093/bioinformatics/btac632"]
LAYOUT = [
    ("networks/", "the 33 genome-scale networks as every tool received them (blocked reactions removed), the same "
                  "networks before removal, the exact coupling matrices the reference read, and their deficiencies"),
    ("results/kimonu/", "the Kimonu runs of 2026-11-11 on the final code. These are the timed runs on exclusive nodes "
                        "(iterated route with detached and attached uptakes, reference-equivalent mode, tangent route, "
                        "small networks), both routes in both uptake modes, coverage counts, ACR and ACRR lists, the gap "
                        "networks, the uniformly drawn flux with its Schwartz-Zippel bound, the rounding check, the "
                        "floating-point comparison and the symbolic baselines"),
    ("results/cocoa_jl/", "COCOA.jl v1.1.0 with efficient=false and efficient=true, three seeds each (frozen records)"),
    ("results/reference/", "the reference implementation of Langary et al. (2025), exhaustive and sampled (frozen records)"),
    ("results/hardware/", "nodes, CPUs, memory, threads and software versions of every benchmark run"),
    ("figures/", "the published figures (SVG), Data S1 and S2 (CSV) and the source data the figures are drawn from"),
    (snapshot, f"a snapshot of the source code and the benchmark scripts at commit {commit}")]
esc = lambda s: s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
desc = (f"<p>Run records and figure data of the article \"{esc(paper)}\" by Philipp C. Schneider (Science Advances, "
        f"DOI {meta['ARTICLE_DOI']}). The archive holds every file that the figures, tables, Data S1 and S2 and the "
        f"numbers of the article are computed from, and its README.md maps each file to the paper item it supports and "
        f"to the script that reads it.</p><p><strong>Abstract of the article.</strong> {esc(abstract)}</p>"
        f"<p><strong>Contents of Kimonu_records.zip.</strong></p><ul>"
        + "".join(f"<li><code>{esc(a)}</code> {esc(b)}</li>" for a, b in LAYOUT)
        + f"</ul><p>SHA256SUMS lists every file with its checksum, and link_records.sh links the records into a checkout "
          f"of the repository, so that the scripts in scripts/ rebuild the figures, tables and numbers. The Kimonu "
          f"results were computed with commit {RUN_COMMIT} of <a href=\"{REPO}\">{REPO}</a>, and later commits up to the "
          f"snapshot restructure the package without changing its results (scripts/kinetic_variants/"
          f"refactor_equivalence.py). Data, records and figures are licensed under CC BY 4.0, and "
          f"the code in the snapshot under the Apache License 2.0.</p>")
record = {"metadata": {
    "upload_type": "dataset", "title": TITLE,
    "creators": [{"name": "Schneider, Philipp C.", "orcid": "0000-0001-8858-428X",
                  "affiliation": "Department of Molecular Biology, Massachusetts General Hospital, Boston, MA, USA"}],
    "description": desc, "access_right": "open", "license": "cc-by-4.0", "keywords": KEYWORDS,
    "related_identifiers": [
        {"identifier": REPO, "relation": "isSupplementTo", "resource_type": "software", "scheme": "url"},
        {"identifier": f"{REPO}/tree/{commit_full}", "relation": "isCompiledBy", "resource_type": "software", "scheme": "url"},
        {"identifier": meta["ARTICLE_DOI"], "relation": "isSupplementTo", "resource_type": "publication-article", "scheme": "doi"}],
    "references": REFERENCES, "version": meta["ARCHIVE_VERSION"], "publication_date": meta["RELEASE_DATE"],
    "language": "eng",
    "notes": f"Built by scripts/make_zenodo_bundle.py of {REPO} at commit {commit}."}}
json.dump(record, open(os.path.join(out, "zenodo_metadata.json"), "w"), indent=2, ensure_ascii=False)
open(os.path.join(out, "description.html"), "w").write(f"<!doctype html><meta charset=utf-8><title>Zenodo description</title>\n{desc}\n")
fill.update({"@ZIPMB@": f"{os.path.getsize(z) / 1e6:.0f}", "@OUT@": os.path.abspath(out), "@TITLE@": TITLE,
             "@KEYWORDS@": ", ".join(f"`{k}`" for k in KEYWORDS)})
render("ZENODO_UPLOAD.md", os.path.join(out, "UPLOAD.md"))
print("wrote zenodo_metadata.json, description.html, UPLOAD.md")
