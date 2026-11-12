#!/usr/bin/env python3
"""Assemble the Zenodo archive of the paper from the records folder.

    python make_zenodo_bundle.py <out dir> [recorded lists ...]      # from scripts/
    bash zenodo_bundle.sh [out dir]          # all steps in one command, on the login node (figures, build, rebuild test)

The archive is the records folder (KIMONU_RECORDS, records.py) as it is, plus README.md (ZENODO_README.md), LICENSE
(ZENODO_LICENSE), CITATION.cff, the source snapshot (git archive of $KIMONU_REV, default HEAD) and SHA256SUMS. PARTS
names every folder of the records with the paper items it feeds, what it holds and the script or SLURM job script that
produced it. Every file of the records must lie in a part, or in figures/ as one of PUBLISHED (Figure1.svg to
Figure4.svg, supp/FigS1_stage_timing.svg to supp/FigS6_float_errors.svg, supp/DataS1_values.csv,
supp/DataS2_uptake_acr.csv), and every part must hold files, or the script stops. drafts/ and experiments/ of the
records (figure drafts, runs that no paper item reads) are left out. Code is not part of the records, except the
wrappers of KEEP_CODE.

The recorded lists are those of figures/record_inputs.py, one per figure or number script, as zenodo_bundle.sh writes
them. Every file a script read must be a file of the archive ("records/...") or a tracked file of the snapshot
("repo/..."), or the script stops. ZENODO_DOI (default the reserved DOI below), ARTICLE_DOI, ARCHIVE_VERSION (the
version of this record) and RELEASE_DATE fill the placeholders of README.md and CITATION.cff. Every archive path and
every repository script that README.md names in backticks must exist, or the script stops. Writes into <out dir>

    top/      README.md, LICENSE, CITATION.cff, Kimonu.py-<commit>.zip, SHA256SUMS (every file of the archive)
    upload/   Kimonu_records.zip (the records and top/ as one deflate zip, top folder Kimonu_records/) and README.md
    zenodo_metadata.json, description.html, UPLOAD.md (ZENODO_UPLOAD.md), when $MANUSCRIPT names the manuscript
"""
import datetime, fnmatch, glob, hashlib, json, os, re, shutil, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records

CODE = (".py", ".pyc", ".sbatch", ".sh", ".R", ".m", ".jl")
RUN_COMMIT = "451a748"          # Kimonu commit of the final runs (provenance_*.txt); later src gives the same results
# The author, the one place to change it: README.md, CITATION.cff, the Zenodo metadata and UPLOAD.md take it from here,
# as they take the title of the article.
AUTHOR = dict(family="Schneider Chuzel", given="Philipp C.", orcid="0000-0001-8858-428X")
PAPER = "Exact nullspaces certify kinetic modules and other mass-action invariants in metabolic networks"   # the title of the article
AFFILIATION = "Department of Molecular Biology, Massachusetts General Hospital Research Institute, Massachusetts General Hospital, Boston, MA 02114, USA"
LEFT_OUT = (records.DRAFTS, records.EXPERIMENTS)

# archive folder, the paper items it feeds, content, and what produced it (paths under scripts/, SLURM job ids)
TE = "scripts/kinetic_variants/timing_exclusive.sbatch"
PARTS = [
    ("networks/unblocked", "every tool, Fig. 3C, D, Tables S2, S3",
     "the 33 networks with blocked reactions removed, the input of every tool (crn.json, triplets of A, Y and S)",
     "`scripts/write_unblocked_crn.py`, run by `scripts/make_unblocked.sbatch` (FVA with SCIP), from `networks/with_blocked_reactions/`"),
    ("networks/with_blocked_reactions", "Materials and Methods",
     "the same networks before blocked-reaction removal, converted from the elementary-step models of the reference",
     "`scripts/export_split_model.py`, from the elementary-step models (.mat) of Langary et al. (2025)"),
    ("networks/reference_coupling_matrices", "input of the reference",
     "exact coupling matrices that the reference's R stage read", "`scripts/matlab_reference/write_coupling_exact.py`"),
    ("networks/deficiency", "Table S2, Supplementary Text S10",
     "deficiency, linkage classes and ranks of the 33 networks (deficiency_33.json) and of the 12 networks of Küken et al. (2022) (deficiency_kuken12.json)",
     "`scripts/deficiency/deficiency_split.py`, `scripts/deficiency/deficiency_kuken.py`"),
    ("networks/e_coli_core", "Fig. 4B, Supplementary Text S11",
     "e_coli_core of BiGG (.mat), the largest network of the symbolic baselines and of the timing arm small",
     "BiGG Models, read by `scripts/kinetic_variants/local_acr_symbolic.py`"),
    ("results/kimonu/timing", "Fig. 4A to D, F, figs. S1 to S4, Data S1, Materials and Methods",
     "timed runs, each one process on one core of an exclusive node, replicates rep1 to rep3 at seeds 0, 1, 2 of r*. "
     "Arms det and att, iterated route (the default) with detached and attached uptakes. Arm ref, reference-equivalent "
     "mode (no r*, one computation repeated). Arms tan and tanatt, tangent route. Arm small, tangent route on the small "
     "networks and the EnvZ-OmpR ladder. Each run has its JSON (modules, ACR, ACRR, seconds) and /usr/bin/time log, "
     "each arm provenance_*.txt and hardware_*.txt, and slurm_jobs.txt lists the SLURM jobs",
     f"`{TE}`, which runs `scripts/run_kimonu_crn.py` (det, att), `scripts/run_kimonu_refmode.py` (ref) and "
     "`scripts/kinetic_variants/tangent_cell.py` (tan, tanatt, small), SLURM jobs 572311 (det, att, ref, small), 572312 (tan), 572313 (tanatt)"),
    ("results/kimonu/uptake_modes", "Fig. 3A, B (through hierarchy_detached), Data S1, Results",
     "both routes with detached and attached uptakes, per network",
     "`scripts/kinetic_variants/uptake_modes.py`, run by `scripts/kinetic_variants/uptake_modes.sbatch`, SLURM job 572113"),
    ("results/kimonu/hierarchy_detached", "Fig. 3A, B", "per network, the counts of both routes that Fig. 3A, B shows (detached uptakes)",
     "`scripts/kinetic_variants/hierarchy_detached.py`, from `results/kimonu/uptake_modes/` and `results/kimonu/member_counts/`"),
    ("results/kimonu/member_counts", "Fig. 3C, D, Tables S2, S3, Data S1, Results",
     "complexes in kinetic modules, metabolites with ACR or in an ACRR group, group sizes",
     "`scripts/kinetic_variants/member_counts.py`, run by `scripts/kinetic_variants/member_counts.sbatch`, SLURM jobs 572115 (detached), 572116 (attached)"),
    ("results/kimonu/flux_balance_vs_mass_action", "Results",
     "reactions coupled under flux balance against rate classes under mass action",
     "`scripts/kinetic_variants/fb_vs_ma.py`, run by `scripts/kinetic_variants/fb_vs_ma.sbatch` (SLURM job 572117), summary by `scripts/kinetic_variants/fb_vs_ma_summary.py`"),
    ("results/kimonu/acr_lists_by_uptake_mode", "Data S2, Results",
     "ACR and ACRR lists of the iterated route per uptake mode, metabolites as indices into the network's own list",
     "`scripts/kinetic_variants/uptake_lists.py`, run by `scripts/kinetic_variants/uptake_lists.sbatch`, SLURM job 572114"),
    ("results/kimonu/gap_roots", "Results, Supplementary Text S6", "the open questions the three gap networks come down to",
     "`scripts/kinetic_variants/gap_roots.py`, run by `scripts/kinetic_variants/gap_roots.sbatch` with `scripts/kinetic_variants/gap_roots_cases.txt`, SLURM job 572118"),
    ("results/kimonu/second_order", "Supplementary Text S12", "exact second-order test modulo a prime on the gap networks",
     "`scripts/kinetic_variants/second_order_modp.py`, run by `scripts/kinetic_variants/second_order_modp.sbatch`, SLURM job 572119"),
    ("results/kimonu/gap_experiments/injectivity", "Supplementary Text S6, S7",
     "sign condition, injectivity and tangent dimensions on the gap networks",
     "`scripts/kinetic_variants/injectivity_gap.py`, run by `scripts/kinetic_variants/injectivity.sbatch`"),
    ("results/kimonu/gap_experiments/descartes", "Supplementary Text S6", "Descartes step on the gap networks",
     "`scripts/kinetic_variants/descartes_vs_tangent.py`, run by `scripts/kinetic_variants/descartes_array.sbatch`"),
    ("results/kimonu/gap_experiments/autonomous_blocks", "Supplementary Text S6", "relations proven by autonomous blocks",
     "`scripts/kinetic_variants/autonomous_blocks.py`, run by `scripts/kinetic_variants/autonomous_blocks.sbatch` (one line per network)"),
    ("results/kimonu/gap_experiments/two_point", "Supplementary Text S6",
     "two steady states at the same rate constants on the gap networks: the relations and bases of the prepare step "
     "(<network>_relations.json, <network>_V.npz, read by the injectivity, autonomous-block and gap-root runs), one "
     "record per relation of iHN637 and the summary two_point_iHN637.json",
     "`scripts/kinetic_variants/two_point_gap.py`, run by `scripts/kinetic_variants/two_point.sbatch`"),
    ("results/kimonu/gap_experiments/gap_reduction", "Supplementary Text S6",
     "the iHN637 witness: the network reduced around the tangent-only relation of glyc3p, stage by stage",
     "`scripts/kinetic_variants/reduce_gap.py`, `reduce_contract.py` and `reduce_stage3.py` in the same folder, each run by the sbatch file of its name"),
    ("results/kimonu/uniform_flux", "Supplementary Text S8, Materials and Methods",
     "both routes at a uniformly drawn r* (64 bits) against the LP point, 33 networks",
     "`scripts/kinetic_variants/sz_flux.py`, run by `scripts/kinetic_variants/sz_flux33.sbatch`, SLURM job 572120"),
    ("results/kimonu/genericity_bound", "Results, Supplementary Text S8",
     "Schwartz-Zippel bound of the uniform r* per network and route (summary.json)",
     "`scripts/kinetic_variants/genericity_bound.py`, from `results/kimonu/uniform_flux/`, SLURM job 572316"),
    ("results/kimonu/lp_point_rounding", "Supplementary Text S8",
     "rounding of the LP point (scales 100, 10^4, 10^5, 10^6, unrounded) on iNJ661 and iHN637, and the same timed on "
     "an exclusive node (exclusive/, three replicates per network)",
     "`scripts/kinetic_variants/scale_check.py` (SLURM jobs 572317, iNJ661, and 572318, iHN637), `scripts/kinetic_variants/scale_check_1e4.py` "
     "(iHN637_1e4.json, job 572319) and `scripts/kinetic_variants/scale_check_exclusive.sbatch` (exclusive/, job 572320)"),
    ("results/kimonu/float_comparison/parts", "Fig. 3E, F, figs. S5, S6, Results",
     "floating-point comparison, per network, nullspace and tolerance",
     "`scripts/kinetic_variants/float_exactcmp.py` (the parts *_exactcmp*), run by `scripts/kinetic_variants/float_exactcmp.sbatch`, and "
     "`scripts/kinetic_variants/float_four.py` (the other parts), run by `scripts/kinetic_variants/float_four.sbatch` with `scripts/kinetic_variants/float_four_parts.txt`"),
    ("results/kimonu/symbolic_local_acr", "Fig. 4B, Supplementary Text S11",
     "symbolic local-ACR criteria (minors, elimination) on the small networks and the ladder",
     "`scripts/kinetic_variants/local_acr_symbolic.py`, run by `scripts/kinetic_variants/local_acr_ladder.sbatch`"),
    ("results/kimonu/preprocessing_time", "Materials and Methods", "preprocessing time of iAF1260b, load to factorization",
     "`scripts/kinetic_variants/preprocess_time.py`, run by `scripts/kinetic_variants/preprocess_time.sbatch`"),
    ("results/kimonu/one_lp_vs_fva", "Materials and Methods", "one LP against FVA, blocked reactions and times",
     "`scripts/kinetic_variants/one_lp.py`, run by `scripts/kinetic_variants/one_lp.sbatch`"),
    ("results/kimonu/upstream_set_route_2026-09", "none (the arm kimonu_up of figdata.json, drawn by no figure)",
     "the upstream-set route of September 2026 (superseded): runs per network (kimonu_runs/), the largest module under "
     "overlap-only merging (overlap_merge_logs/) and ACR/ACRR by the pairwise rule of the reference (shinar_feinberg_rule/)",
     "the upstream-set route of `scripts/compat/upstream.py`, by scripts removed in November 2026 (`scripts/HISTORY.md`); read by `scripts/figures/gather_figdata.py`"),
    ("results/cocoa_jl/efficient_false", "Fig. 4A, C, D, figs. S1, S4, Data S1",
     "COCOA.jl v1.1.0, efficient=false, three seeds (cell logs and module files)",
     "`scripts/run_cocoa_acr.jl` with the Julia project `scripts/julia_v110/`, run by `scripts/cocoa_v110.sbatch`, `scripts/cocoa_v110_seeds.sbatch`, "
     "`scripts/cocoa_v110_backfill.sbatch`, `scripts/cocoa_gap_false.sbatch` and `scripts/cocoa_iaf_false.sbatch` (cell lists `scripts/*_cells.txt`)"),
    ("results/cocoa_jl/efficient_true", "Fig. 4A, C, D, F, figs. S1, S4, Data S1", "COCOA.jl v1.1.0, efficient=true, three seeds",
     "`scripts/run_cocoa_acr.jl` with `scripts/julia_v110/`, run by `scripts/cocoa_eff.sbatch`, `scripts/cocoa_eff_big.sbatch`, `scripts/cocoa_eff_all.sbatch` "
     "and `scripts/cocoa_gap_true.sbatch`"),
    ("results/cocoa_jl/kinetic_stage_iLJ478", "Fig. 4F",
     "COCOA.jl's kinetic stage on iLJ478 from saved concordance labels, three replicates (the 4.5 s of Fig. 4F)",
     "`run.sh` and `kin_only.jl` in this folder (SLURM job 568621), not in the repository"),
    ("results/cocoa_jl/kinetic_stage_on_kimonu_partition", "Supplementary Text",
     "COCOA.jl's kinetic stage run on Kimonu's exact concordance partition, seven networks (comparison_7models.json), and "
     "the isolation of the iCN900 difference (iCN900_isolation/)",
     "`scripts/cocoa_loop/kimonu_side.py`, `cocoa_side.jl` and `compare.py`; `scripts/matlab_reference/cocoa_isolation/isolate.jl`"),
    ("results/reference/exhaustive", "Fig. 4A, C, D, F, figs. S1, S2, Data S1",
     "the reference (Langary et al.), exhaustive, with concordance results, R-stage logs and results, retested partitions",
     "`scripts/matlab_reference/run_reference.m` under Octave, run by `scripts/ref_final.sbatch` (`scripts/final_models.txt`); R stage "
     "`scripts/reference/run_chain.sh` (`scripts/reference/rerun_chain.sbatch`); retests by `scripts/ref_finish.sbatch` with "
     "`scripts/matlab_reference/retest_pairs.py` and `scripts/matlab_reference/retest_rejected.m`"),
    ("results/reference/sampled", "Fig. 4A, C, D, fig. S4, Data S1", "the reference, sampled, three seeds",
     "`scripts/matlab_reference/run_reference.m`, run by `scripts/ref_final_sampled.sbatch` (`scripts/final_sampled_cells.txt`); R stage `scripts/reference/run_chain.sh`"),
    ("results/reference/kinetic_stage_iLJ478", "Fig. 4F", "the reference's R stages on iLJ478, three replicates (the 3.5 min of Fig. 4F)",
     "`run.sh` in this folder (SLURM job 568620), not in the repository"),
    ("results/reference/fork_gate_runs", "Supplementary Text S3",
     "the gate runs of the adaptations of the reference (fixed/, fixed2/), exhaustive, on iIS312 Amastigote and Trypomastigote",
     "`scripts/matlab_reference/run_reference.m`, run by `scripts/ref_fixed_validate.sbatch` and `scripts/ref_fixed2_validate.sbatch`"),
    ("results/hardware", "Materials and Methods",
     "nodes, CPUs, memory reserved and used, threads and software versions of the benchmark runs",
     "`scripts/figures/hardware_record.py`, from sacct, the job logs and the result files"),
    ("figures/source_data", "Fig. 4A to F, figs. S1 to S4, Tables S2, S3, Data S1",
     "figdata.json (the values the figures show, per network, tool and mode), pair_counts.json (relation counts per tool), "
     "pair_sets.json (set containment against Kimonu), seed_spread.json (relations that change with the seed), "
     "reference_rstage_end.json (the end of each R stage of the exhaustive reference, the modification time of its ACR "
     "log, which a zip keeps to 2 s only)",
     "`scripts/figures/gather_figdata.py` (figdata.json), `scripts/figures/pair_sets.py` (pair_counts.json, pair_sets.json), "
     "`scripts/figures/seed_spread.py`, `scripts/figures/reference_rstage_end.py`"),
]
# the published figures, Data S1 and S2 as figures/ of the records holds them, each with the script that draws or writes it
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


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


out, listings = sys.argv[1], sys.argv[2:]
HERE = os.getcwd()
assert os.path.basename(HERE) == "scripts", "run from scripts/"
REC = records.RECORDS
git = lambda *a: subprocess.run(["git", *a], capture_output=True, text=True, check=True).stdout.strip()
REV = os.environ.get("KIMONU_REV", "HEAD")
commit, commit_full = git("rev-parse", "--short=7", REV), git("rev-parse", REV)
meta = {k: os.environ.get(k, d) for k, d in dict(ZENODO_DOI="10.5281/zenodo.24145925", ARTICLE_DOI="",
        ARCHIVE_VERSION="1.0", RELEASE_DATE="2026-11-12").items()}
# Until the article is published, ARTICLE_DOI stays empty and the archive names it as submitted to Science Advances.
ARTICLE_LINE = (f"Article DOI {meta['ARTICLE_DOI']}, archive DOI {meta['ZENODO_DOI']}." if meta["ARTICLE_DOI"]
                else f"Archive DOI {meta['ZENODO_DOI']}. The article is submitted to Science Advances.")
top, upload = os.path.join(out, "top"), os.path.join(out, "upload")
for d in (top, upload):
    assert not os.path.exists(d), f"{d} exists"

# ---- the records: every file in a part or published, the left-out folders counted
every = sorted(os.path.relpath(os.path.join(d, f), REC) for d, _, fs in os.walk(REC, followlinks=True) for f in fs)
left = {d: sum(1 for a in every if a.startswith(d + "/")) for d in LEFT_OUT}
files = [a for a in every if not a.startswith(tuple(d + "/" for d in LEFT_OUT))]
part_of = lambda a: next((d for d, *_ in PARTS if a.startswith(d + "/")), None)
stray = [a for a in files if part_of(a) is None and a not in {f"figures/{p}" for p in PUBLISHED}]
assert not stray, "files of the records in no part of the archive (add a part or move them):\n" + "\n".join(stray[:50])
empty = [d for d, *_ in PARTS if not any(part_of(a) == d for a in files)]
assert not empty, f"parts without files: {empty}"
missing_fig = [p for p in PUBLISHED if f"figures/{p}" not in files]
assert not missing_fig, f"published files missing in {REC}/figures: {missing_fig}"
code = [a for a in files if a.endswith(CODE) and not a.startswith(KEEP_CODE)]
assert not code, "code in the records (it belongs to the repository):\n" + "\n".join(code[:50])
links = [a for a in files if os.path.islink(os.path.join(REC, a))]
assert not links, f"symbolic links in the records: {links[:10]}"

# ---- the recorded lists: every input in the archive or in the snapshot
tracked = set(git("ls-tree", "-r", "--name-only", REV, "--full-tree").split("\n"))
recorded = {l.strip() for f in listings for l in open(f) if l.strip()}
outside = sorted(p for p in recorded if not (p.startswith("records/") and p[8:] in files or p.startswith("repo/") and p[5:] in tracked))
assert not outside, "read by a script but neither in the archive nor in the snapshot:\n" + "\n".join(outside[:50])
read_rec = {p[8:] for p in recorded if p.startswith("records/")}
snap = sorted(p[5:] for p in recorded if p.startswith("repo/"))
unread = {}                                                    # archived files no recorded script read, per part
for a in files:
    if a not in read_rec and part_of(a):
        unread[part_of(a)] = unread.get(part_of(a), 0) + 1

# ---- top/: snapshot, README, LICENSE, CITATION
os.makedirs(top)
snapshot = f"Kimonu.py-{commit}.zip"
subprocess.run(["git", "archive", "--format=zip", f"--prefix=Kimonu.py-{commit}/", "-o", os.path.abspath(os.path.join(top, snapshot)),
                REV], cwd=git("rev-parse", "--show-toplevel"), check=True)     # from the top, the whole repository


def tally(d):
    hit = [a for a in files if part_of(a) == d]
    return f"{len(hit):,} ({sum(os.path.getsize(os.path.join(REC, a)) for a in hit) / 1e6:.2f} MB)"


layout = "| archive folder | content | feeds | produced by | files |\n|---|---|---|---|---|\n"
layout += "".join(f"| `{d}/` | {what} | {items} | {made} | {tally(d)} |\n" for d, items, what, made in PARTS)
layout += "".join(f"| `figures/{f}` | the published file | {item} | {made} | 1 |\n" for f, (item, made) in PUBLISHED.items())
NAME = f"{AUTHOR['given']} {AUTHOR['family']}"                 # as a byline; "family, given" for Zenodo
fill = {"@PAPER@": PAPER, "@AUTHOR@": NAME, "@AUTHOR_FAMILY@": AUTHOR["family"], "@AUTHOR_GIVEN@": AUTHOR["given"], "@ORCID@": AUTHOR["orcid"],
        "@AFFILIATION@": AFFILIATION, "@COMMIT@": commit, "@COMMIT_FULL@": commit_full, "@RUN_COMMIT@": RUN_COMMIT, "@DATE@": datetime.date.today().isoformat(),
        "@LAYOUT@": layout.rstrip(), "@SNAPSHOT@": snapshot, "@YEAR@": meta["RELEASE_DATE"][:4]}
fill.update({f"@{k}@": v for k, v in meta.items()})
fill["@ARTICLE_LINE@"] = ARTICLE_LINE
fill["@ARTICLE_CFF@"] = (f'  journal: "Science Advances"\n  year: {meta["RELEASE_DATE"][:4]}\n  doi: "{meta["ARTICLE_DOI"]}"\n'
                       if meta["ARTICLE_DOI"] else '  journal: "Science Advances"\n  status: submitted\n')
PLACEHOLDER = ("XXXXXXX", "X.Y", "YYYY")


def render(src, dst, t=None):
    t = open(src, encoding="utf-8").read() if t is None else t
    for k, v in fill.items():
        t = t.replace(k, v)
    assert not re.search(r"@[A-Z_]+@", t), f"placeholder left in {dst}"
    open(dst, "w", encoding="utf-8").write(t)


render("ZENODO_README.md", os.path.join(top, "README.md"))


def named_paths_exist(md):
    """every backticked path of the README must exist: in the archive, in the repository at REV, or (a bare name or
    pattern) as the name of a file in either"""
    folders = {os.path.dirname(a) for a in files} | {os.path.dirname(f) for f in tracked}
    names = ({os.path.basename(f) for f in tracked} | {os.path.basename(a) for a in files} | set(os.listdir(top))
             | {os.path.basename(d) for d in folders})
    bad = []
    for tok in sorted(set(re.findall(r"`([^`\s$=()]+)`", md))):
        if "/" not in tok and not re.search(r"\.(py|sbatch|sh|jl|m|R|txt|json|svg|csv|log|zip)$", tok):
            continue                                                       # a word or value, not a path
        parts = re.match(r"(.*)\{([^}]*)\}(.*)", tok)
        for e in ([parts.group(1) + x + parts.group(3) for x in parts.group(2).split(",")] if parts else [tok]):
            e = e.rstrip("/")
            ok = (e in files or any(a.startswith(e + "/") or fnmatch.fnmatch(a, e) for a in files)
                  or any(fnmatch.fnmatch(f, e) or f.startswith(e + "/") for f in tracked) or e in folders
                  or "/" not in e and any(fnmatch.fnmatch(n, e) for n in names))
            if not ok:
                bad.append(e)
    assert not bad, "README.md names paths that exist neither in the archive nor in the repository: " + ", ".join(bad)


render("ZENODO_LICENSE", os.path.join(top, "LICENSE"))
render(None, os.path.join(top, "CITATION.cff"), f"""cff-version: 1.2.0
message: "If you use these records, please cite the article and this archive."
type: dataset
title: 'Run records and figure data for "@PAPER@"'
authors:
  - family-names: "@AUTHOR_FAMILY@"
    given-names: "@AUTHOR_GIVEN@"
    affiliation: "@AFFILIATION@"
    orcid: "https://orcid.org/@ORCID@"
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
  title: "@PAPER@"
  authors:
    - family-names: "@AUTHOR_FAMILY@"
      given-names: "@AUTHOR_GIVEN@"
      orcid: "https://orcid.org/@ORCID@"
@ARTICLE_CFF@""")
named_paths_exist(open(os.path.join(top, "README.md"), encoding="utf-8").read())

# ---- SHA256SUMS over the records and top/, then one zip, reopened and verified member by member
members = [(a, os.path.join(REC, a)) for a in files] + [(f, os.path.join(top, f)) for f in sorted(os.listdir(top))]
with open(os.path.join(top, "SHA256SUMS"), "w") as f:
    for a, p in sorted(members):
        f.write(f"{sha256(p)}  {a}\n")
os.makedirs(upload)
z = os.path.join(upload, "Kimonu_records.zip")
with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
    for a, p in sorted(members) + [("SHA256SUMS", os.path.join(top, "SHA256SUMS"))]:
        zf.write(p, f"Kimonu_records/{a}")
shutil.copy2(os.path.join(top, "README.md"), os.path.join(upload, "README.md"))
sums = dict(reversed(l.rstrip("\n").split("  ", 1)) for l in open(os.path.join(top, "SHA256SUMS")))
with zipfile.ZipFile(z) as zf:
    assert zf.testzip() is None
    names = [n[len("Kimonu_records/"):] for n in zf.namelist()]
    assert sorted(names) == sorted([a for a, _ in members] + ["SHA256SUMS"])
    for a, _ in members:
        assert hashlib.sha256(zf.read(f"Kimonu_records/{a}")).hexdigest() == sums[a], a
size = sum(os.path.getsize(p) for _, p in members)
print(f"commit {commit}; records {REC}: {len(files)} files; archive {len(members) + 1} files, {size / 1e6:.0f} MB; "
      f"{z} {os.path.getsize(z) / 1e6:.0f} MB, verified")
print("left out of the records:", ", ".join(f"{d}/ {n} files" for d, n in left.items()))
print(f"recorded inputs {len(recorded)}: in the archive {len(read_rec)}, in the snapshot only {len(snap)} ({', '.join(snap[:5])}{', ...' if len(snap) > 5 else ''})")
print("archived files that no recorded script read, per part:", ", ".join(f"{d} {n}" for d, n in unread.items()))
held = [f"{k}={v}" for k, v in meta.items() if any(s in v for s in PLACEHOLDER)]
print("placeholders left in README.md and CITATION.cff:", ", ".join(held) or "none")

# ---- Zenodo metadata (API format), the description and the upload steps, next to upload/
MS = os.environ.get("MANUSCRIPT")                              # the manuscript (Markdown), for title and abstract
if not MS:
    print("MANUSCRIPT not set, so no zenodo_metadata.json, description.html, UPLOAD.md")
    sys.exit(0)
ms = open(MS, encoding="utf-8").read()
abstract = ms[ms.index("# Abstract"):].split("\n")[2].strip()
paper = re.search(r"\*\*Title:\*\* (.+)", ms).group(1).strip()
assert paper == PAPER, f"the manuscript's title differs from PAPER of make_zenodo_bundle.py: {paper}"
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
                  "networks before removal, the exact coupling matrices the reference read, their deficiencies and "
                  "e_coli_core for the symbolic baselines"),
    ("results/kimonu/", "the Kimonu runs of 2026-11-11 on the final code. These are the timed runs on exclusive nodes "
                        "(iterated route with detached and attached uptakes, reference-equivalent mode, tangent route, "
                        "small networks), both routes in both uptake modes, coverage counts, ACR and ACRR lists, the gap "
                        "networks, the uniformly drawn flux with its Schwartz-Zippel bound, the rounding check, the "
                        "floating-point comparison and the symbolic baselines"),
    ("results/cocoa_jl/", "COCOA.jl v1.1.0 with efficient=false and efficient=true, three seeds each (frozen records)"),
    ("results/reference/", "the reference implementation of Langary et al. (2025), exhaustive and sampled (frozen records)"),
    ("results/hardware/", "nodes, CPUs, memory, threads and software versions of every benchmark run"),
    ("figures/", "the published figures (SVG), Data S1 and S2 (CSV) and the source data the figures are drawn from"),
    (snapshot, f"a snapshot of the source code and the scripts at commit {commit}")]
esc = lambda s: s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
desc = (f"<p>Run records and figure data of the article \"{esc(paper)}\" by {esc(NAME)} (Science Advances, "
        + (f"DOI {meta['ARTICLE_DOI']}" if meta["ARTICLE_DOI"] else "submitted") + "). The archive holds every file that the figures, tables, Data S1 and S2 and the "
        f"numbers of the article are computed from, and its README.md maps each file to the paper item it supports and "
        f"to the script that reads it.</p><p><strong>Abstract of the article.</strong> {esc(abstract)}</p>"
        f"<p><strong>Contents of Kimonu_records.zip.</strong></p><ul>"
        + "".join(f"<li><code>{esc(a)}</code> {esc(b)}</li>" for a, b in LAYOUT)
        + f"</ul><p>SHA256SUMS lists every file with its checksum. The scripts in scripts/ of the repository read the "
          f"unpacked archive as their records folder (environment variable KIMONU_RECORDS) and rebuild the figures, "
          f"tables and numbers. The Kimonu results were computed with commit {RUN_COMMIT} of <a href=\"{REPO}\">{REPO}</a>, "
          f"and later commits up to the snapshot restructure the package without changing its results "
          f"(scripts/kinetic_variants/refactor_equivalence.py). Data, records and figures are licensed under CC BY 4.0, and "
          f"the code in the snapshot under the Apache License 2.0.</p>")
record = {"metadata": {
    "upload_type": "dataset", "title": TITLE,
    "creators": [{"name": f"{AUTHOR['family']}, {AUTHOR['given']}", "orcid": AUTHOR["orcid"], "affiliation": AFFILIATION}],
    "description": desc, "access_right": "open", "license": "cc-by-4.0", "keywords": KEYWORDS,
    "related_identifiers": [
        {"identifier": REPO, "relation": "isSupplementTo", "resource_type": "software", "scheme": "url"},
        {"identifier": f"{REPO}/tree/{commit_full}", "relation": "isCompiledBy", "resource_type": "software", "scheme": "url"},
    ] + ([{"identifier": meta["ARTICLE_DOI"], "relation": "isSupplementTo", "resource_type": "publication-article",
           "scheme": "doi"}] if meta["ARTICLE_DOI"] else []),
    "references": REFERENCES, "version": meta["ARCHIVE_VERSION"], "publication_date": meta["RELEASE_DATE"],
    "language": "eng",
    "notes": f"Built by scripts/make_zenodo_bundle.py of {REPO} at commit {commit}."}}
json.dump(record, open(os.path.join(out, "zenodo_metadata.json"), "w"), indent=2, ensure_ascii=False)
open(os.path.join(out, "description.html"), "w").write(f"<!doctype html><meta charset=utf-8><title>Zenodo description</title>\n{desc}\n")
fill.update({"@ZIPMB@": f"{os.path.getsize(z) / 1e6:.0f}", "@OUT@": os.path.abspath(out), "@TITLE@": TITLE,
             "@KEYWORDS@": ", ".join(f"`{k}`" for k in KEYWORDS), "@RECORDS@": REC})
render("ZENODO_UPLOAD.md", os.path.join(out, "UPLOAD.md"))
print("wrote zenodo_metadata.json, description.html, UPLOAD.md")
