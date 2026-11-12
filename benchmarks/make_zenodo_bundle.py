#!/usr/bin/env python3
"""Assemble the Zenodo archive of the paper, the records that its figures, tables, Data S1 and S2 and numbers read.

    python make_zenodo_bundle.py <out dir> <figures dir> [recorded lists ...]      # from benchmarks/
    bash zenodo_bundle.sh [out dir]          # all steps in one command, on the login node (records, build, rebuild test)

PARTS maps each folder of the archive to its source (a path under benchmarks/ or on the cluster scratch space), to
where the scripts of benchmarks/ look for it (a path under benchmarks/, or an environment variable) and to the paper
items that read it. <figures dir> holds the figures as published (Figure*.svg, supp/*.svg, supp/DataS1_values.csv,
supp/DataS2_uptake_acr.csv). The recorded lists are those of figures/record_inputs.py, one per figure or number
script, run from benchmarks/ of a copy of the checkout, because the scripts rewrite figdata.json and figures/*.json:

    python figures/record_inputs.py lists/<script>.txt figures/<script>.py [args]

for gather_figdata, pair_sets, seed_spread, dataS1_S2, tableS2_S3 and paper_numbers (Python 3.12 environment
of the runs) and fig4_figS1_S2_benchmark, figS3_S4_uptakes_seeds, figS5_S6_float, fig3AB_hierarchy and fig3CD_coverage (matplotlib
3.8), as zenodo_bundle.sh does. Every file a script read must be in the archive (the published SVGs of <figures dir>
included), or in SNAPSHOT_ONLY, or the script stops. Code is not copied into the data folders. The source snapshot (git archive of $KIMONU_REV, default HEAD) holds it. ZENODO_DOI,
ARTICLE_DOI, RELEASE_TAG, ARCHIVE_VERSION and RELEASE_DATE fill the placeholders of README.md and CITATION.cff when
set (for a rebuild once the DOI is reserved). Writes into <out dir>

    bundle/   README.md (ZENODO_README.md), LICENSE (ZENODO_LICENSE), CITATION.cff, link_records.sh,
              Kimonu.py-<commit>.zip, networks/, results/, figures/, SHA256SUMS (every other file)
    upload/   Kimonu_records.zip (bundle/ as one deflate zip, top folder Kimonu_records/) and README.md
    zenodo_metadata.json, description.html, UPLOAD.md (ZENODO_UPLOAD.md), when $MANUSCRIPT names the manuscript
"""
import datetime, fnmatch, glob, hashlib, json, os, re, shutil, subprocess, sys, zipfile

SCRATCH = "/scratch/sadreyev_lab/pschneider/"
RERUN = SCRATCH + "kimonu_rerun_1111/"
CODE = (".py", ".pyc", ".sbatch", ".sh", ".R", ".m", ".jl")
RUN_COMMIT = "451a748"          # Kimonu commit of the final runs (provenance_*.txt); src changed only in docstrings since
PLACED, ENV = "place", "env"

# archive path, source (a directory, taken whole without code, or a glob), (PLACED, path under benchmarks/) or
# (ENV, variable), paper items, content
PARTS = [
    ("networks/unblocked", "unblocked/*_fixed.crn.json", (PLACED, "unblocked"), "every tool, Table S2",
     "the 33 networks with blocked reactions removed, the input of every tool (crn.json, triplets of A, Y and S)"),
    ("networks/original", "sweep/*_fixed.crn.json", (PLACED, "sweep"), "Materials and Methods",
     "the same networks before blocked-reaction removal, converted from the elementary-step models of the reference"),
    ("networks/coupling_exact", "coupling_exact/*.mat", (PLACED, "coupling_exact"), "reference (input)",
     "exact coupling matrices that the reference's R stage read"),
    ("networks/deficiency", "deficiency/deficiency_33.json", (PLACED, "deficiency"), "Table S2, Supplementary Text S10",
     "deficiency, linkage classes and ranks of the 33 networks"),
    ("results/kimonu/timing", RERUN + "timing", (ENV, "KIMONU_TIMING"), "Fig. 4A, B, figs. S1 to S4, Data S1, Methods",
     "timed runs, each one process on one core of an exclusive node, replicates rep1 to rep3 at seeds 0, 1, 2 of r*. "
     "Arms det and att, iterated route (the default) with detached and attached uptakes. Arm ref, reference-equivalent "
     "mode (no r*, one computation repeated). Arms tan and tanatt, tangent route. Arm small, tangent route on the small "
     "networks and the EnvZ-OmpR ladder. Each run has its JSON (modules, ACR, ACRR, seconds) and /usr/bin/time log, "
     "each arm provenance_*.txt and hardware_*.txt, and slurm_jobs.txt lists the SLURM jobs"),
    ("results/kimonu/uptake_modes", "kinetic_variants/uptake_modes", (PLACED, "kinetic_variants/uptake_modes"),
     "Fig. 3A, B, Results", "both routes with detached and attached uptakes, per network"),
    ("results/kimonu/hierarchy_detached", "kinetic_variants/hierarchy_detached", (PLACED, "kinetic_variants/hierarchy_detached"),
     "Fig. 3A, B", "the counts of Fig. 3A, B, from uptake_modes and member_counts (hierarchy_detached.py)"),
    ("results/kimonu/member_counts", "kinetic_variants/member_counts", (PLACED, "kinetic_variants/member_counts"),
     "Fig. 3C, D, Tables S2, S3, Results", "complexes in kinetic modules, metabolites with ACR or in an ACRR group, group sizes"),
    ("results/kimonu/fb_vs_ma", "kinetic_variants/fb_vs_ma", (PLACED, "kinetic_variants/fb_vs_ma"), "Results",
     "reactions coupled under flux balance against rate classes under mass action"),
    ("results/kimonu/uptake_lists", RERUN + "uptake_lists", (ENV, "UPTAKE_LISTS"), "Data S2, Results",
     "ACR and ACRR lists of the iterated route per uptake mode"),
    ("results/kimonu/gap_roots", RERUN + "gap_roots", (ENV, "GAP_OUT"), "Results, Supplementary Text S6",
     "the open questions the three gap networks come down to"),
    ("results/kimonu/second_order", "kinetic_variants/second_order_results", (PLACED, "kinetic_variants/second_order_results"),
     "Supplementary Text S12", "exact second-order test modulo a prime on the gap networks"),
    ("results/kimonu/gap_experiments/injectivity", "kinetic_variants/injectivity_results", (PLACED, "kinetic_variants/injectivity_results"),
     "Supplementary Text S6, S7", "sign condition, injectivity and tangent dimensions on the gap networks"),
    ("results/kimonu/gap_experiments/descartes", "kinetic_variants/descartes_*.json", (PLACED, "kinetic_variants"),
     "Supplementary Text S6", "Descartes step on the gap networks"),
    ("results/kimonu/gap_experiments/autonomous_blocks", "kinetic_variants/autonomous_blocks_33.jsonl", (PLACED, "kinetic_variants"),
     "Supplementary Text S6", "relations proven by autonomous blocks"),
    ("results/kimonu/uniform_flux", RERUN + "sz_flux_33", (ENV, "SZ_DIR"), "Supplementary Text S8, Methods",
     "both routes at a uniformly drawn r* (64 bits) against the LP point, 33 networks"),
    ("results/kimonu/genericity_bound", RERUN + "genericity_bound", (ENV, "GB_OUT"), "Results, Supplementary Text S8",
     "Schwartz-Zippel bound of the uniform r* per network and route (summary.json)"),
    ("results/kimonu/scale_check", RERUN + "scale_check", (ENV, "SCALE_DIR"), "Supplementary Text S8",
     "rounding of the LP point (scales 100, 10^4, 10^5, 10^6, unrounded) on iNJ661 and iHN637, and the same timed on "
     "an exclusive node (exclusive/, three replicates per network)"),
    ("results/kimonu/float_comparison/parts", SCRATCH + "kimonu_float_four/parts", (ENV, "KF_OUT"),
     "Fig. 3E, F, figs. S5, S6", "floating-point comparison, per network, nullspace and tolerance"),
    ("results/kimonu/symbolic_baselines", SCRATCH + "kimonu_local_acr", (ENV, "KIMONU_LACR"),
     "Fig. 4B, Supplementary Text S11", "symbolic local-ACR criteria (minors, elimination) on the small networks and the ladder"),
    ("results/kimonu/preprocessing_time", SCRATCH + "kimonu_timing_seeds/fva/preprocess_*.log", (ENV, "FVA_OUT"),
     "Materials and Methods", "preprocessing time of iAF1260b, load to factorization"),
    ("results/kimonu/one_lp", "kinetic_variants/one_lp_results.log", (PLACED, "kinetic_variants"), "Materials and Methods",
     "one LP against FVA, blocked reactions and times"),
    ("results/cocoa_jl/efficient_false", "sweep/cocoa_v110", (PLACED, "sweep/cocoa_v110"), "Fig. 4, figs. S1, S4, Data S1",
     "COCOA.jl v1.1.0, efficient=false, three seeds (cell logs and module files)"),
    ("results/cocoa_jl/efficient_true", "sweep/cocoa_eff", (PLACED, "sweep/cocoa_eff"), "Fig. 4, figs. S1, S4, Data S1",
     "COCOA.jl v1.1.0, efficient=true, three seeds"),
    ("results/cocoa_jl/kinetic_stage_iLJ478", SCRATCH + "cocoa_kin_timing_2026-11-03", None, "Fig. 4F",
     "COCOA.jl's kinetic stage on iLJ478 from saved concordance labels, three replicates (the 4.5 s of Fig. 4F), with "
     "the timing wrapper"),
    ("results/reference/exhaustive", "reference_final", (PLACED, "reference_final"), "Fig. 4, figs. S1, S2, Data S1",
     "the reference (Langary et al.), exhaustive, with concordance results, R-stage logs and results, retested partitions"),
    ("results/reference/sampled", "reference_final_sampled", (PLACED, "reference_final_sampled"), "Fig. 4, fig. S4, Data S1",
     "the reference, sampled, three seeds"),
    ("results/reference/kinetic_stage_iLJ478", SCRATCH + "ref_rstage_timing_2026-11-03", None, "Fig. 4F",
     "the reference's R stages on iLJ478, three replicates (the 3.5 min of Fig. 4F), with the timing wrappers"),
    ("results/hardware", "hardware_record.json", (PLACED, "."), "Materials and Methods",
     "nodes, CPUs, memory reserved and used, threads and software versions of the benchmark runs, from sacct, the job "
     "logs and the result files (figures/hardware_record.py)"),
    ("figures/source_data", "figdata.json", (PLACED, "."), "all figures, Data S1", "the values the figures show (gather_figdata.py)"),
    ("figures/source_data", "figures/pair_counts.json", (PLACED, "figures"), "Fig. 4C, D, Data S1", "relation counts per tool"),
    ("figures/source_data", "figures/pair_sets.json", (PLACED, "figures"), "Fig. 4C, D", "set containment against Kimonu"),
    ("figures/source_data", "figures/seed_spread.json", (PLACED, "figures"), "fig. S4", "spread over seeds"),
]
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
assert os.path.basename(HERE) == "benchmarks", "run from benchmarks/"
git = lambda *a: subprocess.run(["git", *a], capture_output=True, text=True, check=True).stdout.strip()
REV = os.environ.get("KIMONU_REV", "HEAD")
commit, commit_full = git("rev-parse", "--short=7", REV), git("rev-parse", REV)
meta = {k: os.environ.get(k, d) for k, d in dict(ZENODO_DOI="10.5281/zenodo.XXXXXXX", ARTICLE_DOI="10.1126/sciadv.XXXXXXX",
        RELEASE_TAG="vX.Y.Z", ARCHIVE_VERSION="X.Y", RELEASE_DATE="YYYY-MM-DD").items()}
bundle, upload = os.path.join(out, "bundle"), os.path.join(out, "upload")
assert not os.path.exists(bundle), f"{bundle} exists"

# ---- collect: archive path -> source, and source -> archive path for the check
plan, where = {}, {}
for dest, src, look, items, what in PARTS:
    hit = files_of(src)
    assert hit, f"no files for {dest} ({src})"
    for p, rel in hit.items():
        if p.endswith(CODE) and dest not in KEEP_CODE or p.endswith(SKIP):
            continue
        a = f"{dest}/{rel}"
        assert a not in plan, f"two sources for {a}"
        plan[a] = p
        where[os.path.relpath(p, HERE) if not p.startswith(SCRATCH) else p] = a
for f in glob.glob(os.path.join(figdir, "*.svg")) + glob.glob(os.path.join(figdir, "supp", "*")):
    if f.endswith((".svg", ".csv")):
        plan[f"figures/{os.path.relpath(f, figdir)}"] = f

# tracked sources must equal REV, which the snapshot holds
paths = sorted({src for _, src, _, _, _ in PARTS if not src.startswith("/")})
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


def loc(src, look):
    """where the scripts of benchmarks/ read a part"""
    if look is None:
        return "no script (the values are quoted in figures/fig4EF_methods_iLJ478.py)"
    if look[0] == ENV:
        return f"`${look[1]}`"
    return f"`{look[1]}/`" if os.path.isdir(src) else f"`{os.path.normpath(os.path.join(look[1], os.path.basename(src)))}`"


def tally(pred):
    hit = [a for a in plan if pred(a)]
    return f"{len(hit):,} | {sum(os.path.getsize(os.path.join(bundle, a)) for a in hit) / 1e6:.2f}"


parts_table = "| archive path | content | paper | where the scripts of benchmarks/ read it |\n|---|---|---|---|\n" + "\n".join(
    f"| `{dest}/` | {what} | {items} | {loc(src, look)} |" for dest, src, look, items, what in PARTS)
counts = "| folder | files | MB |\n|---|---|---|\n" + "\n".join(
    f"| `{d}/` | {tally(lambda a: a.startswith(d + '/'))} |" for d in dict.fromkeys(p[0] for p in PARTS))
counts += f"\n| `figures/` (published figures, Data S1, S2) | {tally(lambda a: a.startswith('figures/') and not a.startswith('figures/source_data/'))} |"
fill = {"@COMMIT@": commit, "@COMMIT_FULL@": commit_full, "@RUN_COMMIT@": RUN_COMMIT, "@DATE@": datetime.date.today().isoformat(),
        "@PARTS@": parts_table, "@COUNTS@": counts, "@SNAPSHOT@": snapshot, "@YEAR@": meta["RELEASE_DATE"][:4]}
fill.update({f"@{k}@": v for k, v in meta.items()})
PLACEHOLDER = ("XXXXXXX", "X.Y", "YYYY")


def render(src, dst, t=None):
    t = open(src, encoding="utf-8").read() if t is None else t
    for k, v in fill.items():
        t = t.replace(k, v)
    assert not re.search(r"@[A-Z_]+@", t), f"placeholder left in {dst}"
    open(dst, "w", encoding="utf-8").write(t)


render("ZENODO_README.md", os.path.join(bundle, "README.md"))
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
# link_records.sh: links the archive into a checkout, where the scripts of benchmarks/ look for each part
link = ["#!/usr/bin/env bash",
        "# Links this archive into a checkout of Kimonu.py (the snapshot Kimonu.py-*.zip, unpacked, or a clone at the same",
        "# commit), where the scripts of benchmarks/ look for each part, and prints the environment variables for the rest.",
        "#   eval \"$(bash link_records.sh <checkout>)\"      then run the scripts from <checkout>/benchmarks",
        "# Files the checkout holds already (tracked in git) are kept and compared with the archive copy.",
        "set -eu", 'A=$(cd "$(dirname "$0")" && pwd); B=$(cd "$1" && pwd)/benchmarks; test -d "$B"',
        'one() { if [ -e "$2" ]; then cmp -s "$1" "$2" || echo "differs from the archive: $2" >&2; '
        'else mkdir -p "$(dirname "$2")"; ln -s "$1" "$2"; fi; }',
        'put() { if [ -d "$A/$1" ]; then (cd "$A/$1" && find . -type f | sed "s|^\\./||") | while read -r f; do one "$A/$1/$f" "$B/$2/$f"; done; '
        'else one "$A/$1" "$B/$2"; fi; }']
envs = {}
for dest, src, look, items, what in PARTS:
    if look is None:
        continue
    if look[0] == PLACED:
        if os.path.isdir(src):
            link.append(f"put {dest} {look[1]}")
        else:
            for p, rel in sorted(files_of(src).items()):
                link.append(f"put {dest}/{rel} {'' if look[1] == '.' else look[1] + '/'}{rel}")
    else:
        envs[look[1]] = dest[:-len("/parts")] if dest.endswith("/parts") else dest
envs.update(KIMONU_DET_DIR=envs["KIMONU_TIMING"] + "/det", KIMONU_ATT_DIR=envs["KIMONU_TIMING"] + "/att",
            KIMONU_REF_DIR=envs["KIMONU_TIMING"] + "/ref")
link.append("cat <<EOF")
link += [f"export {k}=$A/{v}" for k, v in sorted(envs.items())]
link += ["EOF", 'echo "linked into $B" >&2', ""]
open(os.path.join(bundle, "link_records.sh"), "w").write("\n".join(link))
os.chmod(os.path.join(bundle, "link_records.sh"), 0o755)

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
print("scale_check:", ", ".join(a[len("results/kimonu/scale_check/"):] for a in sorted(plan) if a.startswith("results/kimonu/scale_check/")))
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
          f"of the repository, so that the scripts in benchmarks/ rebuild the figures, tables and numbers. The Kimonu "
          f"results were computed with commit {RUN_COMMIT} of <a href=\"{REPO}\">{REPO}</a>, and later commits up to the "
          f"snapshot change the package only in docstrings. Data, records and figures are licensed under CC BY 4.0, and "
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
    "notes": f"Built by benchmarks/make_zenodo_bundle.py of {REPO} at commit {commit}."}}
json.dump(record, open(os.path.join(out, "zenodo_metadata.json"), "w"), indent=2, ensure_ascii=False)
open(os.path.join(out, "description.html"), "w").write(f"<!doctype html><meta charset=utf-8><title>Zenodo description</title>\n{desc}\n")
fill.update({"@ZIPMB@": f"{os.path.getsize(z) / 1e6:.0f}", "@OUT@": os.path.abspath(out), "@TITLE@": TITLE,
             "@KEYWORDS@": ", ".join(f"`{k}`" for k in KEYWORDS)})
render("ZENODO_UPLOAD.md", os.path.join(out, "UPLOAD.md"))
print("wrote zenodo_metadata.json, description.html, UPLOAD.md")
