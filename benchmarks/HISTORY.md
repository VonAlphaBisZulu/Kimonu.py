# Removed files

On 2026-11-11 the files listed below were removed from the working tree to keep `benchmarks/` readable. They
record earlier stages of the work. These are the August pipeline, pilot runs, reference runs before the
adaptations of the fork (pull request 1 of the reference repository), runs on shared nodes, the exploration that
led to the iterated and tangent routes, earlier figure drafts and September notes. Later scripts replaced each of
them, and no result in the paper depends on them. Every path is present at commit `3066b9b` (and up to `7fafd4a`).
To look at or restore one, run `git show 3066b9b:benchmarks/<path>` or `git checkout 3066b9b -- benchmarks/<path>`.
The result folders these scripts wrote stay in place and are listed in `README.md` under "Data folders".

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
  `figures/supp_walkthrough.py`, Supplementary Text S1)
- `docs/paper/md_to_science_tex.py` (replaced by `md_to_sciadv.py`)

## Notes and plans (removed in a later commit the same day, also present at `3066b9b`)
- `../ROADMAP.md` (the plan of August and September 2026)
- `MANIFEST.md`, `BENCHMARK_DESIGN.md` (the benchmark map of 2026-09-17 and the design of August, replaced by
  `README.md`)
- `FINDINGS.md`, `matlab_reference/FINDINGS.md` (dated lab records to 2026-09-24; the adaptations of the reference
  they led to are in pull request 1 of the fork and in Supplementary Text S3)
