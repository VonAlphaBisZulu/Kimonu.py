# Scripts behind the paper

This folder holds every script that produced a figure, table or number of the paper "Exact nullspaces find and
certify kinetic modules and other mass-action invariants in metabolic networks", and the frozen harness that ran
COCOA.jl and the reference implementation of Langary et al. (2025). The Kimonu package itself is in `../src/kimonu`.
No results are kept in git: they are in the records folder, which is uploaded to Zenodo as it is (`ZENODO_README.md`).
Removed and renamed files are listed in `HISTORY.md`.

## The records folder

Every result lives in one folder, the records, in the layout of the Zenodo archive: `networks/` (the 33 networks,
the input of every tool), `results/kimonu/`, `results/cocoa_jl/`, `results/reference/`, `results/hardware/`,
`figures/source_data/` and `figures/` (the published figures, Data S1 and S2). `records.py` names every folder of
it, and every script reads and writes through it. Two environment variables place it.
- `KIMONU_RECORDS`, the records folder. Every script reads its inputs from it, the networks as well as earlier
  results. Default `records/` at the top of the checkout (ignored by git). The unpacked Zenodo archive
  (`Kimonu_records/`) is such a folder. On the cluster `records` is a link to
  `/scratch/sadreyev_lab/pschneider/kimonu_records/`, filled on 2026-11-12 from the frozen results (copies, the
  sources left in place) and identical to the archive of build 572517.
- `KIMONU_OUT`, where the scripts write, in the same layout. Default `KIMONU_RECORDS`. Set it to a new folder to rerun
  a script without overwriting the records. A script reads from `KIMONU_RECORDS` only, so to feed new results into
  later scripts, point `KIMONU_RECORDS` at a folder that holds them.
- Figure drafts go to `drafts/` of `KIMONU_OUT` and runs that no paper item reads to `experiments/`. Neither is part
  of the archive, and `make_zenodo_bundle.py` leaves them out. Any other file of the records outside the folders of
  the archive stops the build.
- In an sbatch file, `python records.py read TIMING det` and `python records.py write UPTAKE_MODES` print a folder
  (and make it, for write). The `#SBATCH --output` lines of the frozen harness name the records folder on scratch,
  because SLURM reads no variables there; with another `KIMONU_OUT`, pass `-o` on the command line.

## How to run

- Run every script from `scripts/`. Most scripts put `scripts/` and `../src` on `sys.path` themselves, and the
  sbatch files add `scripts/` to `PYTHONPATH`.
- Python 3.12 with StrainDesign 1.19.2 (`pip install "straindesign>=1.19.2"`, or on the cluster
  `PYTHONPATH=$HOME/sd-kernel-nonzeros`, the tag v1.19.2), cobra, SciPy and, for `cert_modp.py`, python-flint.
  The cluster environment is `~/miniconda3/envs/straindesign/bin/python`.
- Every `.sbatch` file carries absolute paths to this checkout and its SLURM log path, and writes its results to
  `KIMONU_OUT`. `SD_TREE` selects the StrainDesign tree. Several recorded runs used the tree `~/sd-scip-pricing`
  (1.19.1 plus a branch on SCIP's LP pricing). Their counts agree with those of the timing runs on 1.19.2 (modules,
  largest module, ACR and ACRR of `uptake_modes/` against `figdata.json`, both uptake modes, all 33 networks). The
  rerun uses 1.19.2 throughout.
- Timed runs use `kinetic_variants/timing_exclusive.sbatch`, one process pinned to one core on an exclusive node.

To rerun the Kimonu results (COCOA.jl and the reference are not rerun), with `KIMONU_OUT` set to a new folder, in this
order:
1. `sbatch kinetic_variants/timing_exclusive.sbatch det att ref`, then the arms `tan`, `tanatt` and `small`.
2. `sbatch kinetic_variants/uptake_modes.sbatch`, then `python kinetic_variants/hierarchy_detached.py`.
3. `member_counts.sbatch`, `uptake_lists.sbatch`, `fb_vs_ma.sbatch` (then `fb_vs_ma_summary.py`), `gap_roots.sbatch`,
   `sz_flux33.sbatch` (then `genericity_bound.py`) (all in `kinetic_variants/`). The modules of the reference-equivalent
   mode and the seed check come from the timing arms ref and det.
4. Only if the reference flux or the readout changes, `float_exactcmp.sbatch`, `float_four.sbatch` and
   `second_order_modp.sbatch`.
5. Copy the new folders over a copy of the records, point `KIMONU_RECORDS` at it, rebuild as in `ZENODO_README.md`
   and compare the output of `figures/paper_numbers.py` with the text.

## Paper item, script, input and output

"K" marks a Kimonu run, "F" a frozen COCOA.jl or reference result. Scripts under `figures/` draw, scripts under
`kinetic_variants/` compute. The vault is Phil's manuscript folder, which holds the final figures.
Each drawing, table and data script is named by the paper item it makes: `fig<N><panels>_<what>.py` draws (for
example `fig3AB_hierarchy.py`, `fig4_figS1_S2_benchmark.py`), `fig<N><panels>_update_<what>.py` edits the published
Illustrator file in place, and `figS<n>_...`, `tableS2_S3.py`, `dataS1_S2.py` and `textS1_walkthrough.py` make the
Supplement items. Scripts that gather data for several items are named by what they write (`gather_figdata.py` writes
`figdata.json`, `pair_sets.py` `pair_sets.json` and `pair_counts.json`, `seed_spread.py` `seed_spread.json`,
`hardware_record.py` `hardware_record.json`, `paper_numbers.py` prints the numbers). `HISTORY.md` lists the old names.
Inputs are folders of the records. In the Zenodo archive the same map has the archive folder of every input and the
script or SLURM job script that produced it (`ZENODO_README.md`, filled in by `make_zenodo_bundle.py`).

| paper item | script | input | output |
|---|---|---|---|
| Fig. 1 | Illustrator master, edited by `docs/examples/network/fig1_update_layers.py`, `fig1_update_matrix_labels.py`, `fig1G_update_shared_W.py`, `fig1_fig2_panels.py`. Matrices from `fig1CDE_matrices.py` (SymPy) | none | vault `figures/Figure1.svg` |
| Fig. 2 | `docs/examples/network/fig2_kinase.py` (networks, steady states, values), then `fig2_update_conventions.py`, `fig2C_update_loglog.py`, `fig2E_update_table.py`, `fig1_fig2_panels.py` | closed form | vault `Figure2.svg` |
| Fig. 3A, B (K) | `figures/fig3AB_hierarchy.py`; panel B in k - 1 ratios on the vault SVG by `figures/fig3B_update_ratios.py`, the shared legend without colons by `figures/fig3AB_update_legend.py` | `results/kimonu/hierarchy_detached/` from `kinetic_variants/hierarchy_detached.py` from `uptake_modes/` (`uptake_modes.py`, `.sbatch`) and `member_counts/` | vault `Figure3.svg` (through `figures/fig2_fig3_split.py`) |
| Fig. 3C, D (K) | `figures/fig3CD_coverage.py`, in place on the vault SVG (redraws C, D when they exist) | `results/kimonu/member_counts/` (`member_counts.py`, `.sbatch`), `networks/unblocked/` | vault `Figure3.svg` |
| Fig. 3E, F (K) | `figures/fig3EF_float.py` (uses `figS5_S6_float.draw_certificate`) | `results/kimonu/float_comparison/parts/` from `float_exactcmp.py` (E) and `float_four.py` (F), their `.sbatch` and part lists | vault `Figure3.svg` |
| Fig. 4A (K, F) | `figures/fig4_figS1_S2_benchmark.py` | `figures/source_data/figdata.json` from `figures/gather_figdata.py`, which reads `results/kimonu/timing/{det,att,ref}` (`timing_exclusive.sbatch` running `run_kimonu_crn.py`, `run_kimonu_refmode.py`), `results/cocoa_jl/` and `results/reference/` | vault `Figure4.svg` |
| Fig. 4B (K) | `figures/fig4_figS1_S2_benchmark.py` | `results/kimonu/timing/{tan,small}` (`tangent_cell.py`), `results/kimonu/symbolic_local_acr/` (`local_acr_symbolic.py`, `local_acr_ladder.sbatch`) | same |
| Fig. 4C, D (K, F) | `figures/fig4_figS1_S2_benchmark.py` | `pair_counts.json`, `pair_sets.json` (in `figures/source_data/`) from `figures/pair_sets.py`, which runs Kimonu's default and reads the modules of the reference-equivalent mode in `results/kimonu/timing/ref` and the module files of COCOA.jl and the reference | same |
| Fig. 4E, F (F) | `figures/fig4EF_methods_iLJ478.py`, called by `fig4_figS1_S2_benchmark.py` | `figdata.json`, `results/reference/exhaustive/reference_iLJ478.json`, `results/cocoa_jl/efficient_true/cell_iLJ478_seed1234.log` | same |
| Fig. 4 on new results (K) | `figures/fig4ABF_update.py`, in place on the vault SVG, which carries Phil's Illustrator edits: moves the markers, fit lines and slopes of A, the exact-nullspace markers of B and Kimonu's times and the size of iLJ478 in F | `figdata.json` (new) and that of commit 451a748 (old), the timing arms small and tan of the rerun and of `kimonu_timing_seeds` on scratch | vault `Figure4.svg` |
| Figs. 3, 4 in the words of the text | `figures/fig3_fig4_update_terms.py`, in place on both vault SVGs: nullspace for kernel, networks for models, the iIS312 stages as in `figures/names.py` | none | vault `Figure3.svg`, `Figure4.svg` |
| fig. S1 (K, F) | `figures/fig4_figS1_S2_benchmark.py` | `results/kimonu/timing/ref`, `results/reference/exhaustive/acr/*.log`, `reference_rstage_end.json` (the modification times of the ACR logs, from `figures/reference_rstage_end.py`), COCOA.jl logs | `FigS1_stage_timing.svg` (`drafts/supp/` of `KIMONU_OUT`) |
| fig. S2 (K, F) | `figures/fig4_figS1_S2_benchmark.py` (asserts that the reference-equivalent mode equals the exhaustive reference) | `figdata.json` | `FigS2_reference_agreement.svg` |
| fig. S3 (K) | `figures/figS3_S4_uptakes_seeds.py` | `figdata.json` (arms `kimonu`, `kimonu_det`) | `FigS3_uptakes.svg` |
| fig. S4 (K, F) | `figures/figS3_S4_uptakes_seeds.py` | `seed_spread.json` from `figures/seed_spread.py`, which reads the COCOA.jl and sampled-reference seeds and the three seeds of Kimonu's timing runs, `results/kimonu/timing/det` | `FigS4_seeds.svg` |
| figs. S5, S6 (K) | `figures/figS5_S6_float.py` | `results/kimonu/float_comparison/parts/` | `FigS5_certificate.svg`, `FigS6_float_errors.svg` |
| Text S1 (K) | `figures/textS1_walkthrough.py` | Kimonu at r* = K·1 | Markdown for the Supplement |
| Tables S2, S3 | `figures/tableS2_S3.py` | `networks/unblocked/`, `networks/deficiency/deficiency_33.json` (`deficiency_split.py`), `figdata.json`, `results/kimonu/member_counts/` | Markdown tables |
| Data S1 (K, F) | `figures/dataS1_S2.py` | `figdata.json`, `pair_counts.json`, Kimonu's pairs and tangent-route rows from `results/kimonu/uptake_modes/`, `member_counts/` and `timing/{tan,ref}` | `DataS1_values.csv` (`figures/supp/` of `KIMONU_OUT`) |
| Data S2 (K) | `figures/dataS1_S2.py` | `results/kimonu/acr_lists_by_uptake_mode/` (`uptake_lists.py`, `.sbatch`) | `DataS2_uptake_acr.csv` |
| numbers of Results, Methods and Supplementary Text | `figures/paper_numbers.py` | all of the above | one line per number with its source |
Runs behind single sentences, all in `kinetic_variants/` and read by `paper_numbers.py`: `fb_vs_ma.py` with
`fb_vs_ma_summary.py` (flux coupling under mass action), `gap_roots.py` with `gap_roots.sbatch` and
`gap_roots_cases.txt` (the open questions), `genericity_bound.py` (the Schwartz–Zippel bound of the uniform r*), `scale_check.py` and
`scale_check_1e4.py` (the rounding of the LP point, S8), `sz_flux.py` with
`sz_flux.sbatch` and `sz_flux33.sbatch` (the uniformly drawn flux), `one_lp.py` (one LP against FVA,
`one_lp_results.log`), `preprocess_time.py` with `crn_to_sbml.py` (preprocessing time), `fva_time.sbatch` with
`fva_one_process.py` (FVA time), `second_order_modp.py` (second-order test, `results/kimonu/second_order/`).

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
- Inputs. `export_split_model.py` converted the reference's elementary-step models (.mat) into
  `networks/with_blocked_reactions/`. `write_unblocked_crn.py` (`make_unblocked.sbatch`, FVA with SCIP) removed the
  blocked reactions and wrote `networks/unblocked/`, the input of every tool. `matlab_reference/write_coupling_exact.py`
  wrote `networks/reference_coupling_matrices/`, the exact coupling matrices that the reference's R stage reads.
- COCOA.jl v1.1.0. `run_cocoa_acr.jl` with the Julia project `julia_v110/`, run by `cocoa_v110*.sbatch`
  (efficient=false), `cocoa_eff*.sbatch` (efficient=true), `cocoa_gap_*.sbatch` and `cocoa_iaf_false.sbatch`, with
  their cell lists (`*_cells.txt`). Results in `results/cocoa_jl/efficient_false/` and `efficient_true/`.
- Reference. Its code is not in this repository. `reference/fetch_upstream.sh` clones the fork
  https://github.com/VonAlphaBisZulu/Upstream_Algorithm at 8a6ed17, the head of pull request 1, into
  `reference/Upstream_Algorithm/` (git-ignored; code only, `--full` adds the models and published results, about 2 GB,
  best on scratch), and every script that runs the reference takes the clone from `UPSTREAM_ALGORITHM`, default that
  folder. `matlab_reference/run_reference.m` runs it under Octave with the shims in `matlab_reference/compat/` and the
  CPLEX bridge `matlab_reference/compat/lp_server.py`, run by `ref_final.sbatch` (exhaustive, `final_models.txt`) and
  `ref_final_sampled.sbatch` (`final_sampled_cells.txt`). The R stage is `reference/run_chain.sh`, which runs the
  fork's two R files unmodified on `acr/` of the reference's records (`reference/run_all_sampled.sh` for every
  sampled stem, `reference/rerun_chain.sbatch` reran it where a coupling matrix was missing). `ref_finish.sbatch` with `matlab_reference/retest_pairs.py` and `retest_rejected.m` retested the
  rejected pairs. `ref_fixed_validate.sbatch`, `ref_fixed2_validate.sbatch`, `matlab_reference/test_*.m` and
  `test_kinetic_na_indexing.R` test the adaptations. `patch_mat_names.m` adds metabolite names for the R stage.
- The settings of the reference that needed adjusting at genome scale are in pull request 1 of the fork
  (https://github.com/ankueken/Upstream_Algorithm/pull/1) and in Supplementary Text S3. `compare.py` is an August
  module kept because `matlab_reference/retest_pairs.py` imports it.

## Experiments behind Supplement statements

One-off runs whose result the Supplement states in words, with their records in `results/kimonu/gap_experiments/`
(and `results/reference/fork_gate_runs/`, `networks/deficiency/deficiency_kuken12.json`). `kinetic_variants/README.md`
records them by date. The one-off checks write to `experiments/` of `KIMONU_OUT`.
- Two steady states at the same rate constants. `two_point.py`, `two_point_gap.py`, `two_point_local.py`.
- Sign condition. `injectivity.py`, `injectivity_gap.py`.
- Descartes step. `descartes.py`, `descartes_vs_tangent.py`.
- Shrinking the iHN637 gap. `reduce_gap.py`, `reduce_contract.py`, `reduce_stage3.py`.
- Basis of the tolerance comparison. `ortho_vs_rref.py`.
- Deficiency one. `cocoa_loop/deficiency_one_check.py`. COCOA.jl's kinetic stage on Kimonu's partition, `cocoa_loop/`
  and `matlab_reference/cocoa_isolation/isolate.jl` (records in `results/cocoa_jl/kinetic_stage_on_kimonu_partition/`). Containment of the upstream sets in the iterated route,
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
- `records.py` (the folders of the records), `figures/record_inputs.py` (records the files a script reads) and
  `make_zenodo_bundle.py` (the Zenodo archive, the records folder with `ZENODO_README.md`, `ZENODO_LICENSE`, a
  generated `CITATION.cff`, `SHA256SUMS` and the source snapshot, and the upload steps `ZENODO_UPLOAD.md`). The build
  checks that every file a recorded script read is in the archive or in the snapshot. `zenodo_bundle.sh` runs it all
  in one command (`bash scripts/zenodo_bundle.sh`, on the login node), from copying the vault figures into the records
  to rebuilding the figures, tables and numbers from the finished zip.
- `docs/paper/` (`make_bib.py`, `kimonu.bib`, `md_to_sciadv.py`, `md_to_kimonu_tex.py`, the manuscript in LaTeX).

## Tests

`../tests/` is the package test suite (pytest, run by CI), including the tests of the reference-equivalent mode
(`test_reference.py`). `kinetic_variants/test_tangent_one_shot.py` tests the code of the Supplement experiments
(it needs CPLEX and is not run by CI).

## Results on the cluster

The records folder (`/scratch/sadreyev_lab/pschneider/kimonu_records/`, layout above) holds every result the paper
reads. It was filled from these frozen folders on scratch (`/scratch/sadreyev_lab/pschneider/`), which stay where they
are: the paper's runs on the final code in `kimonu_rerun_1111/` (timing arms det, att, ref, tan, tanatt, small;
uptake_lists, gap_roots, sz_flux_33, genericity_bound, scale_check; STATUS.md), `kimonu_float_four/parts`,
`kimonu_local_acr`, `kimonu_timing_seeds/fva`, `kimonu_twopoint`, `cocoa_kin_timing_2026-11-03` and
`ref_rstage_timing_2026-11-03`, and from the result files git tracked until 2026-11-12 (`HISTORY.md`). Earlier runs
that no paper item reads stay on scratch: `kimonu_timing_seeds/` (det, att, tan, tanatt, small), `kimonu_timing_final/ref`,
`kimonu_uptake_modes/lists`, `kimonu_cert_modp`, `kimonu_gaproots`, `kimonu_seeds`,
`kimonu_trip_2026-11-02/refmode_modules`, `kimonu_second_order`, `kimonu_sz_flux`, `kimonu_fb_vs_ma`,
`kimonu_member_counts`. The Zenodo bundles of 2026-11-06 and of the rerun are in `kimonu_zenodo` and
`kimonu_zenodo_2026-11`.
