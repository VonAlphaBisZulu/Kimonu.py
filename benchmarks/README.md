# Benchmarks and the records behind the paper

This folder holds every script that produced a figure, table or number of the paper "Exact nullspaces find and
certify kinetic modules and other mass-action invariants in metabolic networks", the frozen harness that ran
COCOA.jl and the reference implementation of Langary et al. (2025), and the small result files. The Kimonu package
itself is in `../src/kimonu`. Large results live on the cluster scratch space and in the Zenodo archive
(`ZENODO_README.md`). Files removed on 2026-11-11 are listed in `HISTORY.md`.

## How to run

- Run every script from `benchmarks/`. Most scripts put `benchmarks/` and `../src` on `sys.path` themselves, and the
  sbatch files add `benchmarks/` to `PYTHONPATH`.
- Python 3.12 with StrainDesign 1.19.2 (`pip install "straindesign>=1.19.2"`, or on the cluster
  `PYTHONPATH=$HOME/sd-kernel-nonzeros`, the tag v1.19.2), cobra, SciPy and, for `cert_modp.py`, python-flint.
  The cluster environment is `~/miniconda3/envs/straindesign/bin/python`.
- Every `.sbatch` file carries absolute paths to this checkout and writes to a directory under
  `/scratch/sadreyev_lab/pschneider/` named in its header. Most take an environment variable for the output
  directory (`TIMING_OUT`, `SO_OUT`, `SEED_OUT`, `GAP_OUT`, ...) and `SD_TREE` for the StrainDesign tree. Several
  recorded runs used the tree `~/sd-scip-pricing` (1.19.1 plus a branch on SCIP's LP pricing). Their counts agree with
  those of the timing runs on 1.19.2 (modules, largest module, ACR and ACRR of `uptake_modes/` against `figdata.json`,
  both uptake modes, all 33 networks). The rerun uses 1.19.2 throughout.
- Timed runs use `kinetic_variants/timing_exclusive.sbatch`, one process pinned to one core on an exclusive node.

To rerun the Kimonu results (COCOA.jl and the reference are not rerun), in this order:
1. `sbatch kinetic_variants/timing_exclusive.sbatch det att ref`, then the arms `tan`, `tanatt` and `small`, each with
   a new `TIMING_OUT`.
2. `sbatch kinetic_variants/uptake_modes.sbatch`, then `python kinetic_variants/hierarchy_detached.py`.
3. `member_counts.sbatch`, `uptake_lists.sbatch`, `fb_vs_ma.sbatch` (then `fb_vs_ma_summary.py`), `gap_roots.sbatch`,
   `sz_flux33.sbatch` (then `genericity_bound.py`) (all in `kinetic_variants/`). The modules of the reference-equivalent
   mode and the seed check come from the timing arms ref and det.
4. Only if the reference flux or the readout changes, `float_exactcmp.sbatch`, `float_four.sbatch` and
   `second_order_modp.sbatch`.
5. Rebuild as in `ZENODO_README.md`, with the `KIMONU_*` variables pointing at the new directories, and compare
   the output of `figures/paper_numbers.py` with the text.

## Paper item, script, input and output

"K" marks a Kimonu run, "F" a frozen COCOA.jl or reference result. Scripts under `figures/` draw, scripts under
`kinetic_variants/` compute. The vault is Phil's manuscript folder, which holds the final figures.
Each drawing, table and data script is named by the paper item it makes: `fig<N><panels>_<what>.py` draws (for
example `fig3AB_hierarchy.py`, `fig4_figS1_S2_benchmark.py`), `fig<N><panels>_update_<what>.py` edits the published
Illustrator file in place, and `figS<n>_...`, `tableS2_S3.py`, `dataS1_S2.py` and `textS1_walkthrough.py` make the
Supplement items. Scripts that gather data for several items are named by what they write (`gather_figdata.py` writes
`figdata.json`, `pair_sets.py` `pair_sets.json` and `pair_counts.json`, `seed_spread.py` `seed_spread.json`,
`hardware_record.py` `hardware_record.json`, `paper_numbers.py` prints the numbers). `HISTORY.md` lists the old names.
In the Zenodo archive the same map has the archive folder of every input and the script or SLURM job script that
produced it (`ZENODO_README.md`, filled in by `make_zenodo_bundle.py`).

| paper item | script | input | output |
|---|---|---|---|
| Fig. 1 | Illustrator master, edited by `docs/examples/network/fig1_update_layers.py`, `fig1_update_matrix_labels.py`, `fig1G_update_shared_W.py`, `fig1_fig2_panels.py`. Matrices from `fig1CDE_matrices.py` (SymPy) | none | vault `figures/Figure1.svg` |
| Fig. 2 | `docs/examples/network/fig2_kinase.py` (networks, steady states, values), then `fig2_update_conventions.py`, `fig2C_update_loglog.py`, `fig2E_update_table.py`, `fig1_fig2_panels.py` | closed form | vault `Figure2.svg` |
| Fig. 3A, B (K) | `figures/fig3AB_hierarchy.py`; panel B in k - 1 ratios on the vault SVG by `figures/fig3B_update_ratios.py`, the shared legend without colons by `figures/fig3AB_update_legend.py` | `kinetic_variants/hierarchy_detached/` from `hierarchy_detached.py` from `uptake_modes/` (`uptake_modes.py`, `.sbatch`) and `member_counts/` | vault `Figure3.svg` (through `figures/fig2_fig3_split.py`) |
| Fig. 3C, D (K) | `figures/fig3CD_coverage.py`, in place on the vault SVG (redraws C, D when they exist) | `kinetic_variants/member_counts/` (`member_counts.py`, `.sbatch`), `unblocked/` | vault `Figure3.svg` |
| Fig. 3E, F (K) | `figures/fig3EF_float.py` (uses `figS5_S6_float.draw_certificate`) | scratch `kimonu_float_four/parts` from `float_exactcmp.py` (E) and `float_four.py` (F), their `.sbatch` and part lists | vault `Figure3.svg` |
| Fig. 4A (K, F) | `figures/fig4_figS1_S2_benchmark.py` | `figdata.json` from `figures/gather_figdata.py`, which reads scratch `kimonu_rerun_1111/timing/{det,att,ref}` (`timing_exclusive.sbatch` running `run_kimonu_crn.py`, `run_kimonu_refmode.py`) and `sweep/cocoa_v110`, `sweep/cocoa_eff`, `reference_final`, `reference_final_sampled` | vault `Figure4.svg` |
| Fig. 4B (K) | `figures/fig4_figS1_S2_benchmark.py` | scratch `kimonu_rerun_1111/timing/{tan,small}` (`tangent_cell.py`), `kimonu_local_acr` (`local_acr_symbolic.py`, `local_acr_ladder.sbatch`) | same |
| Fig. 4C, D (K, F) | `figures/fig4_figS1_S2_benchmark.py` | `figures/pair_counts.json`, `pair_sets.json` from `figures/pair_sets.py`, which runs Kimonu's default and reads the modules of the reference-equivalent mode in scratch `kimonu_rerun_1111/timing/ref` and the module files of COCOA.jl and the reference | same |
| Fig. 4E, F (F) | `figures/fig4EF_methods_iLJ478.py`, called by `fig4_figS1_S2_benchmark.py` | `figdata.json`, `reference_final/reference_iLJ478.json`, `sweep/cocoa_eff/cell_iLJ478_seed1234.log` | same |
| Fig. 4 on new results (K) | `figures/fig4ABF_update.py`, in place on the vault SVG, which carries Phil's Illustrator edits: moves the markers, fit lines and slopes of A, the exact-nullspace markers of B and Kimonu's times and the size of iLJ478 in F | `figdata.json` (new) and that of commit 451a748 (old), the timing arms small and tan of the rerun and of `kimonu_timing_seeds` | vault `Figure4.svg` |
| Figs. 3, 4 in the words of the text | `figures/fig3_fig4_update_terms.py`, in place on both vault SVGs: nullspace for kernel, networks for models, the iIS312 stages as in `figures/names.py` | none | vault `Figure3.svg`, `Figure4.svg` |
| fig. S1 (K, F) | `figures/fig4_figS1_S2_benchmark.py` | `KIMONU_REF_DIR`, `reference_final/acr/*.log`, COCOA.jl logs | `figures/supp/FigS1_stage_timing.svg` |
| fig. S2 (K, F) | `figures/fig4_figS1_S2_benchmark.py` (asserts that the reference-equivalent mode equals the exhaustive reference) | `figdata.json` | `figures/supp/FigS2_reference_agreement.svg` |
| fig. S3 (K) | `figures/figS3_S4_uptakes_seeds.py` | `figdata.json` (arms `kimonu`, `kimonu_det`) | `figures/supp/FigS3_uptakes.svg` |
| fig. S4 (K, F) | `figures/figS3_S4_uptakes_seeds.py` | `figures/seed_spread.json` from `figures/seed_spread.py`, which reads the COCOA.jl and sampled-reference seeds and the three seeds of Kimonu's timing runs, scratch `kimonu_rerun_1111/timing/det` | `figures/supp/FigS4_seeds.svg` |
| figs. S5, S6 (K) | `figures/figS5_S6_float.py` | scratch `kimonu_float_four/parts` | `figures/supp/FigS5_certificate.svg`, `FigS6_float_errors.svg` |
| Text S1 (K) | `figures/textS1_walkthrough.py` | Kimonu at r* = K·1 | Markdown for the Supplement |
| Tables S2, S3 | `figures/tableS2_S3.py` | `unblocked/`, `deficiency/deficiency_33.json` (`deficiency_split.py`), `figdata.json`, `kinetic_variants/member_counts/` | Markdown tables |
| Data S1 (K, F) | `figures/dataS1_S2.py` | `figdata.json`, `figures/pair_counts.json`, Kimonu's pairs and tangent-route rows from `kinetic_variants/uptake_modes/`, `member_counts/` and scratch `kimonu_rerun_1111/timing/{tan,ref}` | `DataS1_values.csv` |
| Data S2 (K) | `figures/dataS1_S2.py` | scratch `kimonu_rerun_1111/uptake_lists` (`uptake_lists.py`, `.sbatch`) | `DataS2_uptake_acr.csv` |
| numbers of Results, Methods and Supplementary Text | `figures/paper_numbers.py` | all of the above | one line per number with its source |

Runs behind single sentences, all in `kinetic_variants/` and read by `paper_numbers.py`: `fb_vs_ma.py` with
`fb_vs_ma_summary.py` (flux coupling under mass action), `gap_roots.py` with `gap_roots.sbatch` and
`gap_roots_cases.txt` (the open questions), `genericity_bound.py` (the Schwartz–Zippel bound of the uniform r*), `scale_check.py` and
`scale_check_1e4.py` (the rounding of the LP point, S8), `sz_flux.py` with
`sz_flux.sbatch` and `sz_flux33.sbatch` (the uniformly drawn flux), `one_lp.py` (one LP against FVA,
`one_lp_results.log`), `preprocess_time.py` with `crn_to_sbml.py` (preprocessing time), `fva_time.sbatch` with
`fva_one_process.py` (FVA time), `second_order_modp.py` (second-order test, `second_order_results/`).

## Paper scripts (Kimonu runs)

- Runners. `run_kimonu_crn.py` runs concordance and the iterated route (`kimonu.iterated_route`) on a crn.json
  (uptakes attached unless `--detach-uptakes`) and holds the crn.json loader (`triplets_to_sparse`, `to_exact`) that
  most scripts import, and `crn_of`, the `kimonu.CRN` of a script's own A, Y and exact N. The scripts build the
  `kimonu.CRN` of a crn.json with `kimonu.crn_from_matrices(A, Y)`, whose N = Y·A equals the S of the file on all 33
  networks. `run_kimonu_refmode.py` times the stages of the reference-equivalent mode
  (`kimonu.reference_route`, upstream sets merged by overlap only and the Shinar–Feinberg rule).
  `kinetic_variants/tangent_cell.py` times the tangent route.
- `kinetic_variants/`. `timing_exclusive.sbatch`, `uptake_modes.py`, `hierarchy_detached.py`, `member_counts.py`,
  `uptake_lists.py`, `fb_vs_ma.py`, `fb_vs_ma_summary.py`, `gap_roots.py`, `cert_modp.py`, `float_exactcmp.py`,
  `float_four.py`, `sz_flux.py`, `one_lp.py`, `preprocess_time.py`, `second_order_modp.py`, `seed_invariance.py`, each
  with its `.sbatch`. `models33.txt` and `models_gap.txt` are the array indices.
- `figures/`. `gather_figdata.py`, `pair_sets.py`, `seed_spread.py`, `fig4_figS1_S2_benchmark.py`, `fig4EF_methods_iLJ478.py`,
  `fig3AB_hierarchy.py`, `fig3CD_coverage.py`, `fig3EF_float.py`, `fig2_fig3_split.py` (the one-off split of the
  combined figure into Figures 2 and 3), `figS5_S6_float.py`, `figS3_S4_uptakes_seeds.py`, `textS1_walkthrough.py`,
  `dataS1_S2.py`, `tableS2_S3.py`, `paper_numbers.py`, `fig3_fig4_update_terms.py`. `names.py` holds the one display name
  per network that every figure, table and data file uses.
- `deficiency/deficiency_split.py` (deficiency, linkage classes and ranks of the 33 networks).
- `docs/examples/network/` (the edit scripts of Figures 1 and 2).

## Frozen comparison harness (COCOA.jl, reference)

These ran in August and September and are kept as the record. Their results cannot be rerun in the time of the
paper, so do not change them or their inputs.
- Inputs. `export_split_model.py` converted the reference's elementary-step models (.mat) into `sweep/*_fixed.crn.json`.
  `write_unblocked_crn.py` (`make_unblocked.sbatch`, FVA with SCIP) removed the blocked reactions and wrote
  `unblocked/`, the input of every tool. `matlab_reference/write_coupling_exact.py` wrote `coupling_exact/`, the exact
  coupling matrices that the reference's R stage reads.
- COCOA.jl v1.1.0. `run_cocoa_acr.jl` with the Julia project `julia_v110/`, run by `cocoa_v110*.sbatch`
  (efficient=false), `cocoa_eff*.sbatch` (efficient=true), `cocoa_gap_*.sbatch` and `cocoa_iaf_false.sbatch`, with
  their cell lists (`*_cells.txt`). Results in `sweep/cocoa_v110/` and `sweep/cocoa_eff/` (git-ignored).
- Reference. `matlab_reference/run_reference.m` under Octave with the shims in `matlab_reference/compat/` and the
  CPLEX bridge `matlab_reference/compat/lp_server.py`, run by `ref_final.sbatch` (exhaustive, `final_models.txt`) and
  `ref_final_sampled.sbatch` (`final_sampled_cells.txt`). The R stage is `reference_final/acr/run_chain.sh` with the
  fork's R code in `reference_final*/acr/Code/` (`reference_final/rerun_chain.sbatch` reran it where a coupling matrix
  was missing). `ref_finish.sbatch` with `matlab_reference/retest_pairs.py` and `retest_rejected.m` retested the
  rejected pairs. `ref_fixed_validate.sbatch`, `ref_fixed2_validate.sbatch`, `matlab_reference/test_*.m` and
  `test_kinetic_na_indexing.R` test the adaptations. `patch_mat_names.m` adds metabolite names for the R stage.
- The settings of the reference that needed adjusting at genome scale are in pull request 1 of the fork
  (https://github.com/ankueken/Upstream_Algorithm/pull/1) and in Supplementary Text S3. `compare.py` is an August
  module kept because `matlab_reference/retest_pairs.py` imports it.

## Experiments behind Supplement statements

One-off runs whose result the Supplement states in words. `kinetic_variants/README.md` records them by date.
- Two steady states at the same rate constants. `two_point.py`, `two_point_gap.py`, `two_point_local.py`.
- Sign condition. `injectivity.py`, `injectivity_gap.py`.
- Descartes step. `descartes.py`, `descartes_vs_tangent.py`.
- Shrinking the iHN637 gap. `reduce_gap.py`, `reduce_contract.py`, `reduce_stage3.py`.
- Basis of the tolerance comparison. `ortho_vs_rref.py`.
- Deficiency one. `cocoa_loop/deficiency_one_check.py`. COCOA.jl's kinetic stage on Kimonu's partition, `cocoa_loop/`
  and `matlab_reference/cocoa_isolation/isolate.jl`. Containment of the upstream sets in the iterated route,
  `cocoa_loop/proof_upstream_in_psi.md`.
- One-off checks with no number in the paper. `autonomous_blocks.py`, `binomial_formulas.py`, `relation_values.py`,
  `terminal_join.py`, `zero_module.py`, `single_round.py`, `lumped_toric.py`, `toric_check.py`,
  `natural_example_search.py`, `deficiency/deficiency_kuken.py`, the August diagnosis scripts in `matlab_reference/`,
  `compare_to_published.py`, `validate_sweep.py`, `validate_acr_published.py`, `cocoa_s5000.sbatch`.
- Imported modules of earlier stages. `kinetic_variants/tangent_one_shot.py` and `tangent_fast.py` (the first exact
  readouts of the tangent route, used by the scripts above and by `second_order_modp.py` and `local_acr_symbolic.py`),
  `kinetic_variants/branch_certificate.py`, and `compat/upstream.py` (the extra merges of the upstream-set route of
  September on top of `kimonu.reference`, used by `validate_acr_published.py`,
  `matlab_reference/compare_kinetic_partitions.py`, `kinetic_variants/autonomous_blocks.py` and `terminal_join.py`).

## Tools and loaders

- `figures/svg_finish.py`, `svg_panels.py` (fonts, size, panel groups of every SVG), `svg_same.py` (two drawings of
  the same figure compared, text and pixels, for the rebuild check), `fonts_single_family.py`,
  `docs/examples/network/latex_text.py` (formulas as outlines).
- `figures/record_inputs.py` (records the files a script reads) and `make_zenodo_bundle.py` (the Zenodo archive, with
  `ZENODO_README.md`, `ZENODO_LICENSE` and the upload steps `ZENODO_UPLOAD.md`). The archive checks that it holds every
  file the recorded scripts read, and `link_records.sh` in it links the records into a checkout. `zenodo_bundle.sh`
  runs it all in one command (`bash benchmarks/zenodo_bundle.sh`, on the login node), from copying the vault figures
  to rebuilding the figures, tables and numbers from the finished zip.
- `docs/paper/` (`make_bib.py`, `kimonu.bib`, `md_to_sciadv.py`, `md_to_kimonu_tex.py`, the manuscript in LaTeX).

## Tests

`../tests/` is the package test suite (pytest, run by CI), including the tests of the reference-equivalent mode
(`test_reference.py`). `kinetic_variants/test_tangent_one_shot.py` tests the code of the Supplement experiments
(it needs CPLEX and is not run by CI).

## Data folders

| folder | content | status |
|---|---|---|
| `unblocked/`, `sweep/*_fixed.crn.json`, `coupling_exact/` | the 33 networks and the coupling matrices the reference read | input, frozen |
| `figdata.json`, `figures/pair_counts.json`, `pair_sets.json`, `seed_spread.json` | the values the figures show | paper (K, F) |
| `kinetic_variants/uptake_modes/`, `hierarchy_detached/`, `member_counts/`, `fb_vs_ma/`, `second_order_results/`, `one_lp_results.log` | Kimonu results small enough for git | paper (K) |
| `deficiency/deficiency_33.json` | deficiency of the 33 networks | paper (Table S2) |
| `reference_final/`, `reference_final_sampled/` | the reference, exhaustive and sampled (stage logs git-ignored) | paper (F) |
| `sweep/cocoa_v110/`, `sweep/cocoa_eff/` | COCOA.jl v1.1.0 (git-ignored) | paper (F) |
| `reference_fixed/`, `reference_fixed2/` | gate runs of the fork | record behind Supplementary Text S3 |
| `kinetic_variants/injectivity_results/`, `two_point_*.json`, `witness_*.json`, `descartes_*.json`, `deficiency/deficiency_kuken12.json` | Supplement experiments | Supplement |
| `kinetic_variants/hierarchy/`, `local_acr_results/`, `tangent_float_results/`, `logs/`, `closure_*`, `fast_*`, `psi_*`, `s210_*`, `tangent_*` | September and October experiments (`s210_*` still read by `gather_figdata.py` for an arm no figure draws) | superseded |
| `sweep_shared/`, `sweep_tentative/` | upstream-set route of September (read by `gather_figdata.py` for the same arm) | superseded |
| `triplicate*/`, `sweep_unblocked/`, `seedcheck_psi/`, `processed/`, `published/`, `results/`, `archive/`, `sweep/cocoa/`, `sweep/cocoa_s5000/`, `sweep/modules/`, `sweep/kimonu_*` (mostly git-ignored) | August and September runs | superseded |
| `reference/`, `reference_full/`, `reference_sampled/`, `reference_sampled_seeds/`, `reference_acr/`, `reference_glpk_broken/` | reference runs before the adaptations of the fork | superseded, never reported |
| `figures/main/`, `figures/prelim/`, `figures/supp/` | figure drafts of September and October (the default output folders of the drawing scripts; the published figures are in the vault and the Zenodo archive) | superseded |

Kimonu results on scratch (`/scratch/sadreyev_lab/pschneider/`). The paper's runs on the final code are in `kimonu_rerun_1111/` (timing
arms det, att, ref, tan, tanatt, small; uptake_modes, member_counts, fb_vs_ma, uptake_lists, gap_roots, second_order,
sz_flux_33, genericity_bound, scale_check; STATUS.md). Earlier: `kimonu_timing_seeds/` (det, att, tan, tanatt, small,
fva), `kimonu_timing_final/ref`, `kimonu_uptake_modes/lists`, `kimonu_float_four/parts`, `kimonu_cert_modp`,
`kimonu_gaproots`, `kimonu_seeds`, `kimonu_trip_2026-11-02/refmode_modules`, `kimonu_second_order`, `kimonu_local_acr`,
`kimonu_sz_flux`, `kimonu_sz_flux_33`, `kimonu_fb_vs_ma`, `kimonu_member_counts`. The Zenodo bundle of 2026-11-06 is in
`kimonu_zenodo`, that of the rerun in `kimonu_zenodo_2026-11`.
