# Removed and renamed files

On 2026-11-12 the folder `benchmarks/` was renamed `scripts/` (git mv, history kept), because it holds every script
of the paper, not only the benchmarks. Paths in this file are relative to that folder. In commits before the rename it
is `benchmarks/`, so a file of an older commit is `git show <commit>:benchmarks/<path>`, and `git log --follow
scripts/<path>` shows the history of a file under both names.

On 2026-11-11 the files listed below were removed from the working tree to keep `benchmarks/` readable. They
record earlier stages of the work. These are the August pipeline, pilot runs, reference runs before the
adaptations of the fork (pull request 1 of the reference repository), runs on shared nodes, the exploration that
led to the iterated and tangent routes, earlier figure drafts and September notes. Later scripts replaced each of
them, and no result in the paper depends on them. Every path is present at commit `3066b9b` (and up to `7fafd4a`).
To look at or restore one, run `git show 3066b9b:benchmarks/<path>` or `git checkout 3066b9b -- benchmarks/<path>` (then `git mv` it into
`scripts/`).
The result folders these scripts wrote are in the records folder or, superseded, at the commit named in "Results out
of git" below.

Four superseded modules stay because live scripts import them. `compare.py` is imported by
`matlab_reference/retest_pairs.py`, `kinetic_variants/tangent_fast.py` and `tangent_one_shot.py` by the
second-order test, the symbolic baselines and the Supplement experiments, and `kinetic_variants/branch_certificate.py`
by `kinetic_variants/test_tangent_one_shot.py`.

## August pipeline (before the elementary-step models)
- `run_kimonu.py`, `run_cocoa.jl`, `prepare_models.jl`, `export_crn.jl`, `staged_sweep.py`, `three_way.py`
- `supplement_figure.jl`, `archive_results.sh`
- `julia/Project.toml`, `julia_revision/Manifest.toml`, `julia_revision/Project.toml`,
  `julia_solvers/Manifest.toml`, `julia_solvers/Project.toml` (COCOA v1.0.0 and development projects; the paper
  used `julia_v110/`)

## COCOA v1.0.0 and pilot runs (never reported)
- `cocoa_determinism.sbatch`, `cocoa_fullnode.sbatch`, `cocoa_one.sbatch`, `cocoa_pilot.sbatch`,
  `cocoa_sweep.sbatch`, `cocoa_unblocked.sbatch`, `inf517_conc.sbatch`, `inf517_modules.sbatch`
- `cocoa_loop/cocoa_loop.sbatch`
- `v110_ext_cells.txt`, `v110_one_1.txt`, `v110_one_2.txt` (cell lists no script reads)

## Reference runs before the adaptations of the fork
- `reference_sampled.sbatch`, `reference_sweep.sbatch`, `ref_full.sbatch`, `ref_full_serial.sbatch`,
  `ref_sampled_seed.sbatch`, `refseed_cells.txt`
- `reference_acr/km_fixed.R`, `reference_acr/run_chain.sh` (R stage on six networks, replaced by
  `reference_final/acr/`)
- `matlab_reference/tentative_upstream.py` (four-phase variant of the upstream sets)
- `matlab_reference/PR_DRAFT.md` (draft of pull request 1, which now exists in the reference repository)

## Kimonu runs on shared nodes (replaced by `kinetic_variants/timing_exclusive.sbatch`)
- `kimonu_triplicate.sbatch`, `kimonu_refmode.sbatch`
- `kinetic_variants/one_core.sbatch`, `kinetic_variants/one_core_cells.txt`

## Exploration that led to the iterated and tangent routes (2026-09-21 to 2026-10-05)
- `kinetic_variants/psi_closure.py`, `psi_extras.py`, `psi_nullspace.py`, `psi_validate.py`, `psi_jobs.sbatch`
- `kinetic_variants/giant_merge.py`, `measure_phases.py`
- `kinetic_variants/tangent_vs_iterated.py`, `tangent_array.sbatch`, `fast_vs_old.py`, `fast_array.sbatch`
- `kinetic_variants/branch_global.py`, `branch_search.py`, `certify_model.py`, `certify_model.sbatch` (sign
  certificates, replaced by `two_point.py`)
- `kinetic_variants/hierarchy_data.py`, `hierarchy.sbatch` (attached uptakes only, replaced by `uptake_modes.py`
  and `hierarchy_detached.py`)

## Floating-point and second-order predecessors
- `kinetic_variants/float_tradeoff.py`, `float_tradeoff.sbatch`, `tangent_float.py`, `tangent_float.sbatch`
  (replaced by `float_four.py` and `float_exactcmp.py`)
- `kinetic_variants/second_order.py`, `second_order.sbatch` (float64, replaced by `second_order_modp.py` and
  `second_order_modp.sbatch`)

## Figure drafts and earlier edit rounds
- `figures/make_main.py`, `make_prelim.py`, `make_figure2.py`, `make_tangent_benchmark.py`,
  `make_deficiency_panel.py`, `pair_counts.py` (draft figures, the combined Figure 2 of October, and the pair
  counter that `pair_sets.py` replaced)
- `docs/examples/network/complex_network.py`, `even_offsets.py`, `figure1_equations.py`, `figure1_harpoons.py`,
  `figure1_latex.py`, `figure1_legend2.py`, `figure1_modules.py`, `figure1_round5.py`, `figure1_round6.py`,
  `figure2_update.py` (earlier edit rounds of the Illustrator masters of Figures 1 and 2)

## September pages and builders
- `docs/figures_page.html`, `docs/pipeline_ledger.html`, `docs/theory_methods.html`
- `docs/examples/figure1_gen.py`, `langary_example.html`, `net_fork.py`, `net_reduced.py`,
  `nullspace_walkthrough.html`, `walkthrough_gen.py` (walkthrough of the Langary et al. example, replaced by
  `figures/textS1_walkthrough.py`, Supplementary Text S1)
- `docs/paper/md_to_science_tex.py` (replaced by `md_to_sciadv.py`)

## Notes and plans (removed in a later commit the same day, also present at `3066b9b`)
- `../ROADMAP.md` (the plan of August and September 2026)
- `MANIFEST.md`, `BENCHMARK_DESIGN.md` (the benchmark map of 2026-09-17 and the design of August, replaced by
  `README.md`)
- `FINDINGS.md`, `matlab_reference/FINDINGS.md` (dated lab records to 2026-09-24; the adaptations of the reference
  they led to are in pull request 1 of the fork and in Supplementary Text S3)

## Removed with the refactor of the package (2026-11-11, present at `f534bf3`)
- `kinetic_variants/s210_criterion.py`. Its function `s210`, the pairwise Shinar–Feinberg rule of the reference, is
  now `kimonu.reference.shinar_feinberg`. Its main ran the upstream-set route of September (`compat/upstream.py` with
  `Y`) and wrote `kinetic_variants/s210_*.json`, which `figures/gather_figdata.py` still reads for the `kimonu_up` arm.
- `compat/test_upstream.py` moved to `../tests/test_reference.py`, as tests of `kimonu.reference_route`. Its tests of
  the September merges (`Y` given) were dropped.

## Renamed on 2026-11-12 (git mv, history kept)
Every drawing, table and data script is named by the paper item it produces, and every supplementary figure file
carries its number. `git log --follow <new path>` shows the history under the old name.

| old name | new name |
|---|---|
| `figures/make_hierarchy.py` | `figures/fig3AB_hierarchy.py` |
| `figures/figure3_ratios.py` | `figures/fig3B_update_ratios.py` |
| `figures/figure3_legend.py` | `figures/fig3AB_update_legend.py` |
| `figures/extend_figure3.py` | `figures/fig3CD_coverage.py` |
| `figures/make_figure2_bottom.py` | `figures/fig3EF_float.py` |
| `figures/split_figure2.py` | `figures/fig2_fig3_split.py` |
| `figures/make_figure3.py` | `figures/fig4_figS1_S2_benchmark.py` |
| `figures/make_fig3_methods.py` | `figures/fig4EF_methods_iLJ478.py` |
| `figures/figure4_update.py` | `figures/fig4ABF_update.py` |
| `figures/figure34_terms.py` | `figures/fig3_fig4_update_terms.py` |
| `figures/make_supp_figures.py` | `figures/figS3_S4_uptakes_seeds.py` |
| `figures/make_float_four.py` | `figures/figS5_S6_float.py` |
| `figures/make_table_s2.py` | `figures/tableS2_S3.py` |
| `figures/make_data_s1_s2.py` | `figures/dataS1_S2.py` |
| `figures/supp_walkthrough.py` | `figures/textS1_walkthrough.py` |
| `docs/examples/network/figure1_lower.py` | `docs/examples/network/fig1CDE_matrices.py` |
| `docs/examples/network/figure1_layers.py` | `docs/examples/network/fig1_update_layers.py` |
| `docs/examples/network/figure1_matrix_labels.py` | `docs/examples/network/fig1_update_matrix_labels.py` |
| `docs/examples/network/figure1_panelG_shared_W.py` | `docs/examples/network/fig1G_update_shared_W.py` |
| `docs/examples/network/figure12_panels.py` | `docs/examples/network/fig1_fig2_panels.py` |
| `docs/examples/network/figure2.py` | `docs/examples/network/fig2_kinase.py` |
| `docs/examples/network/figure2_conventions.py` | `docs/examples/network/fig2_update_conventions.py` |
| `docs/examples/network/figure2_loglog.py` | `docs/examples/network/fig2C_update_loglog.py` |
| `docs/examples/network/figure2_panelE.py` | `docs/examples/network/fig2E_update_table.py` |
| `docs/examples/network/figure3_phil_2026-10-06.svg` | `docs/examples/network/fig2AB_layout_phil_2026-10-06.svg` |
| `figures/main/Figure3_draft.*` (October draft of Figure 4) | `figures/main/Figure4_draft.*` |
| `figures/supp/FigS_stage_timing.*`, `FigS_reference_agreement.*` (October drafts) | `figures/supp/FigS1_stage_timing.*`, `FigS2_reference_agreement.*` |

The figure files the scripts write were renamed with them: `FigS_stage_timing`, `FigS_reference_agreement`,
`FigS_uptakes`, `FigS_seeds`, `FigS_certificate` and `FigS_float_errors` are `FigS1_stage_timing` to
`FigS6_float_errors`, `Figure3_draft` (Figure 4) is `Figure4_draft`, `Figure3CD_draft` (Figure 4E, F) is
`Figure4EF_draft`, `Fig_certificate` (Figure 3E, F) is `Figure3EF_draft` and `hierarchy` (Figure 3A, B) is
`Figure3AB_draft`.

## One records folder for every result (2026-11-12)
Every script now reads its inputs from the records folder (`KIMONU_RECORDS`) and writes its results to `KIMONU_OUT`, in
the layout of the Zenodo archive (`records.py`). These variables of earlier scripts are gone, each replaced by a folder
of the records.

| variable | folder of the records |
|---|---|
| `KIMONU_TIMING`, `TIMING_OUT`, `KIMONU_DET_DIR`, `KIMONU_ATT_DIR`, `KIMONU_REF_DIR` | `results/kimonu/timing/` (arms det, att, ref, ...) |
| `UM_OUT`, `UM_DIR` | `results/kimonu/uptake_modes/` |
| `MC_OUT`, `MC_DIR` | `results/kimonu/member_counts/` |
| `HD_OUT`, `HIERARCHY_DIR` | `results/kimonu/hierarchy_detached/` |
| `FB_OUT` | `results/kimonu/flux_balance_vs_mass_action/` |
| `UPTAKE_LISTS`, `UL_OUT` | `results/kimonu/acr_lists_by_uptake_mode/` |
| `GAP_OUT` | `results/kimonu/gap_roots/` |
| `SO_OUT` | `results/kimonu/second_order/` |
| `INJ_OUT` | `results/kimonu/gap_experiments/injectivity/` |
| `TWO_POINT_OUT`, `TWO_POINT_LOCAL_OUT` | `results/kimonu/gap_experiments/two_point/` (and its `local/`) |
| `SZ_OUT`, `SZ_DIR` | `results/kimonu/uniform_flux/` |
| `GB_OUT` | `results/kimonu/genericity_bound/` |
| `SCALE_DIR` | `results/kimonu/lp_point_rounding/` |
| `KF_OUT` | `results/kimonu/float_comparison/parts/` |
| `KIMONU_LACR`, `LACR_OUT` | `results/kimonu/symbolic_local_acr/` |
| `FVA_OUT` | `results/kimonu/preprocessing_time/` (`experiments/fva_time/` for the FVA networks) |
| `CM_OUT`, `BF_OUT`, `RV_OUT`, `SEED_OUT`, `SR_OUT`, `ZM_OUT`, `REFMOD_OUT` | `experiments/<script>/` (runs no paper item reads) |

`link_records.sh` of the archive is gone as well, since the scripts read the unpacked archive directly.

## Results out of git (2026-11-12)
Git keeps code, documentation and the small inputs the scripts need (cell and model lists, the Julia project, the
compat shims, hand-made figure inputs). Every result file it tracked was removed from the working tree on 2026-11-12,
1,361 files and 312 MB. Each is present at commit `509f968` (`git show 509f968:scripts/<path>`). The 891 files that the
paper reads, or that back a Supplement statement, are in the records folder (and the Zenodo archive), byte-identical,
at these places.

| in git until 509f968 (under `scripts/`) | folder of the records |
|---|---|
| `unblocked/*_fixed.crn.json` | `networks/unblocked/` |
| `sweep/*_fixed.crn.json` | `networks/with_blocked_reactions/` |
| `coupling_exact/*.mat` | `networks/reference_coupling_matrices/` |
| `deficiency/deficiency_33.json`, `deficiency_kuken12.json` | `networks/deficiency/` |
| `kinetic_variants/uptake_modes/`, `hierarchy_detached/`, `member_counts/` | `results/kimonu/` (same names) |
| `kinetic_variants/fb_vs_ma/` | `results/kimonu/flux_balance_vs_mass_action/` |
| `kinetic_variants/second_order_results/` | `results/kimonu/second_order/` |
| `kinetic_variants/injectivity_results/`, `descartes_*.json`, `autonomous_blocks_33.jsonl` | `results/kimonu/gap_experiments/` (`injectivity/`, `descartes/`, `autonomous_blocks/`) |
| `kinetic_variants/two_point_*.json`, `witness_*.json` | `results/kimonu/gap_experiments/two_point/`, `gap_reduction/` |
| `kinetic_variants/one_lp_results.log` | `results/kimonu/one_lp_vs_fva/` |
| `sweep_shared/`, `sweep_tentative/`, `kinetic_variants/s210_*.json` | `results/kimonu/upstream_set_route_2026-09/` (`kimonu_runs/`, `overlap_merge_logs/`, `shinar_feinberg_rule/`) |
| `reference_final/`, `reference_final_sampled/` | `results/reference/exhaustive/`, `sampled/` |
| `reference_fixed/`, `reference_fixed2/` | `results/reference/fork_gate_runs/fixed/`, `fixed2/` |
| `cocoa_loop/comparison_7models.json`, `matlab_reference/cocoa_isolation/iCN900_kimonu_partition.json` | `results/cocoa_jl/kinetic_stage_on_kimonu_partition/` |
| `hardware_record.json` | `results/hardware/` |
| `figdata.json`, `figures/pair_counts.json`, `pair_sets.json`, `seed_spread.json`, `reference_rstage_end.json` | `figures/source_data/` |

The other 470 files were superseded, and no script reads them. These are `kinetic_variants/hierarchy/` (attached
uptakes, October), `local_acr_results/`, `tangent_float_results/`, `closure_*`, `fast_*`, `psi_*` and `tangent_*`
JSON files (September and October experiments), `natural_example_hits.json` and `terminal_join_33.jsonl` (one-off
checks), `sweep/kimonu_*_fixed.json` and `sweep/sweep.json` (August), `reference/`, `reference_full/`,
`reference_sampled/`, `reference_sampled_seeds/` and `reference_acr/` (reference runs before the adaptations of the
fork, never reported), the figure drafts `figures/main/`, `figures/prelim/` and `figures/supp/` (the published figures
are in the vault and the archive), and the drafts and pages of `docs/examples/` (`Figure1_draft.*`,
`Figure2_toy_draft.*`, `figure1_top_draft.svg`, `figure1_phil_2026-10-14.svg`, `langary_network_*.svg`,
`net_fork.svg`, `walkthrough.css`). The git-ignored results of the checkout (`sweep/cocoa_v110/`, `sweep/cocoa_eff/`,
the stage logs of `reference_final*/`, the August and September folders `triplicate*/`, `sweep_unblocked/`,
`archive/`, `processed/`, `published/`, `results/`, `reference/`, `seedcheck_psi/` and the logs) were moved to
`/scratch/sadreyev_lab/pschneider/kimonu_checkout_leftovers_2026-11-12/scripts/`, the first three after being copied
into the records.

## The reference's code is fetched, not copied (2026-11-12)
The four copies of the fork's R code, `reference_final/acr/Code/kinetic_modules/` and
`reference_final_sampled/acr/Code/kinetic_modules/` (`code_kineticModule_analysis.R`,
`Get_excel_tables_kinetic_modules.R`), were removed. They were identical to the fork
https://github.com/VonAlphaBisZulu/Upstream_Algorithm at 8a6ed17 (pull request 1 of ankueken/Upstream_Algorithm), and
`reference/fetch_upstream.sh` clones that commit. The R stage runner moved with its two variants merged:
`reference_final/acr/run_chain.sh` and `reference_final_sampled/acr/run_chain.sh` are `reference/run_chain.sh`
(exhaustive or sampled), `reference_final_sampled/acr/run_all.sh` is `reference/run_all_sampled.sh` and
`reference_final/rerun_chain.sbatch` is `reference/rerun_chain.sbatch`. No change to the reference exists only in this
repository: the corrected negative branch that `matlab_reference/test_corrected_branch.m` tries and the indexing that
`test_kinetic_na_indexing.R` patches are fork commits 8cde9eb and 28350a9.
