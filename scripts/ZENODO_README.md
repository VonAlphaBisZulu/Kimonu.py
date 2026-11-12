# Run records and figure data for "Exact nullspaces find and certify kinetic modules and other mass-action invariants in metabolic networks"

Philipp C. Schneider, Department of Molecular Biology, Massachusetts General Hospital, Boston, MA, USA
(ORCID 0000-0001-8858-428X).
Article DOI @ARTICLE_DOI@, archive DOI @ZENODO_DOI@.
Built on @DATE@ by `benchmarks/zenodo_bundle.sh` and `benchmarks/make_zenodo_bundle.py` of Kimonu.py at commit @COMMIT@.

## What this archive holds

This archive holds every file that the figures, tables, Data S1 and S2 and the numbers of the paper are computed from.
These are the 33 input networks, the Kimonu runs of 2026-11-11 on the final code, the runs of COCOA.jl v1.1.0 and of
the reference implementation of Langary et al. (2025), which are frozen records of August to November 2026, the figure
source data, the published figures with Data S1 and S2, and a snapshot of the source code (`@SNAPSHOT@`). The scripts
are in the repository https://github.com/VonAlphaBisZulu/Kimonu.py and in the snapshot. `SHA256SUMS` lists every file
with its SHA-256 checksum.

Paths in this file that start with `benchmarks/` or `src/` are paths of the repository. Paths that start with
`networks/`, `results/` or `figures/` are paths of this archive. Each script is named by the paper item it makes
(`benchmarks/figures/fig3AB_hierarchy.py` draws Fig. 3A, B, `figS3_S4_uptakes_seeds.py` figs. S3 and S4,
`tableS2_S3.py` Tables S2 and S3, `dataS1_S2.py` Data S1 and S2). Scripts named `*_update_*` edit the published
figure in place, because Figures 1 to 4 carry hand edits made in Adobe Illustrator.

Superseded runs are left out. Two kinds of files that the scripts read are in the source snapshot only, where git
tracks them. These are the list of the 33 networks (`benchmarks/kinetic_variants/models33.txt`) and the upstream-set
arm of September, which `benchmarks/figures/gather_figdata.py` still reads into `figdata.json` but no figure, table or
number shows (`benchmarks/sweep_shared/`, `benchmarks/sweep_tentative/`, `benchmarks/kinetic_variants/s210_*.json`).
The records of the Supplement experiments that the text states in words, without a number, are also in the snapshot
only (`benchmarks/kinetic_variants/two_point_*.json` and similar, listed in `benchmarks/README.md`).

## Paper item by item

Each paper item, the script that draws or computes it, the files of this archive it reads, and the script or SLURM
job script that produced each of them. The table "Layout" below gives the same for every folder.

| paper item | drawn or computed by | archive files it reads | each produced by |
|---|---|---|---|
| Fig. 1 | Illustrator master (`figures/Figure1.svg`). Matrices by `benchmarks/docs/examples/network/fig1CDE_matrices.py`, edits by `fig1_update_layers.py`, `fig1_update_matrix_labels.py`, `fig1G_update_shared_W.py` and `fig1_fig2_panels.py` in the same folder | none (the example network, in closed form) | none |
| Fig. 2 | Illustrator master (`figures/Figure2.svg`). Panels by `benchmarks/docs/examples/network/fig2_kinase.py`, edits by `fig2_update_conventions.py`, `fig2C_update_loglog.py`, `fig2E_update_table.py` and `fig1_fig2_panels.py` in the same folder, and split from the combined figure of early November by `benchmarks/figures/fig2_fig3_split.py` | none (closed form) | none |
| Fig. 3A, B | `benchmarks/figures/fig3AB_hierarchy.py`, then on the published file `fig3B_update_ratios.py` and `fig3AB_update_legend.py` | `results/kimonu/hierarchy_detached/` | `benchmarks/kinetic_variants/hierarchy_detached.py`, from `results/kimonu/uptake_modes/` (`benchmarks/kinetic_variants/uptake_modes.sbatch`, job 572113) and `results/kimonu/member_counts/` (`benchmarks/kinetic_variants/member_counts.sbatch`, jobs 572115, 572116) |
| Fig. 3C, D | `benchmarks/figures/fig3CD_coverage.py`, on the published file | `results/kimonu/member_counts/`, `networks/unblocked/` | `benchmarks/kinetic_variants/member_counts.sbatch`; `benchmarks/write_unblocked_crn.py` |
| Fig. 3E, F | `benchmarks/figures/fig3EF_float.py`, with `draw_certificate` of `figS5_S6_float.py` | `results/kimonu/float_comparison/parts/` | `benchmarks/kinetic_variants/float_exactcmp.sbatch` and `benchmarks/kinetic_variants/float_four.sbatch` |
| Fig. 4A | `benchmarks/figures/fig4_figS1_S2_benchmark.py`, then on the published file `fig4ABF_update.py` and `fig3_fig4_update_terms.py` | `figures/source_data/figdata.json` | `benchmarks/figures/gather_figdata.py`, from `results/kimonu/timing/{det,att,ref}` (`benchmarks/kinetic_variants/timing_exclusive.sbatch`, job 572311), `results/cocoa_jl/efficient_false/`, `results/cocoa_jl/efficient_true/` (`benchmarks/cocoa_v110*.sbatch`, `benchmarks/cocoa_eff*.sbatch`), `results/reference/exhaustive/` and `results/reference/sampled/` (`benchmarks/ref_final.sbatch`, `benchmarks/ref_final_sampled.sbatch`) |
| Fig. 4B | as Fig. 4A | `results/kimonu/timing/{tan,small}`, `results/kimonu/symbolic_local_acr/` | `benchmarks/kinetic_variants/timing_exclusive.sbatch` (jobs 572312, 572311); `benchmarks/kinetic_variants/local_acr_ladder.sbatch` |
| Fig. 4C, D | as Fig. 4A | `figures/source_data/pair_counts.json`, `figures/source_data/pair_sets.json` | `benchmarks/figures/pair_sets.py`, from `results/kimonu/timing/ref` and the module files in `results/cocoa_jl/` and `results/reference/` |
| Fig. 4E, F | `benchmarks/figures/fig4EF_methods_iLJ478.py`, drawn into Figure 4 by `fig4_figS1_S2_benchmark.py` | `figures/source_data/figdata.json`, `results/reference/exhaustive/reference_iLJ478.json`, `results/cocoa_jl/efficient_true/cell_iLJ478_seed1234.log`, `results/cocoa_jl/kinetic_stage_iLJ478/`, `results/reference/kinetic_stage_iLJ478/` | `benchmarks/figures/gather_figdata.py`; `benchmarks/ref_final.sbatch`; `benchmarks/cocoa_eff*.sbatch`; the `run.sh` in each `kinetic_stage_iLJ478/` folder, whose times the script holds as constants |
| fig. S1 | `benchmarks/figures/fig4_figS1_S2_benchmark.py` (`figures/supp/FigS1_stage_timing.svg`) | `results/kimonu/timing/ref`, `results/reference/exhaustive/acr/*.log`, `figures/source_data/reference_rstage_end.json`, the COCOA.jl logs of `results/cocoa_jl/efficient_false/` and `results/cocoa_jl/efficient_true/`, `figures/source_data/figdata.json` | `benchmarks/kinetic_variants/timing_exclusive.sbatch` (job 572311); `benchmarks/reference_final/acr/run_chain.sh`; `benchmarks/figures/reference_rstage_end.py` (the end times of the R stages, which a zip cannot carry as file times); `benchmarks/cocoa_v110*.sbatch`, `benchmarks/cocoa_eff*.sbatch`; `benchmarks/figures/gather_figdata.py` |
| fig. S2 | `benchmarks/figures/fig4_figS1_S2_benchmark.py` (`figures/supp/FigS2_reference_agreement.svg`) | `figures/source_data/figdata.json` | `benchmarks/figures/gather_figdata.py` |
| fig. S3 | `benchmarks/figures/figS3_S4_uptakes_seeds.py` (`figures/supp/FigS3_uptakes.svg`) | `figures/source_data/figdata.json`, arms kimonu and kimonu_det | `benchmarks/figures/gather_figdata.py`, from `results/kimonu/timing/{att,det}` |
| fig. S4 | `benchmarks/figures/figS3_S4_uptakes_seeds.py` (`figures/supp/FigS4_seeds.svg`) | `figures/source_data/seed_spread.json` | `benchmarks/figures/seed_spread.py`, from `results/kimonu/timing/det`, `results/cocoa_jl/` and `results/reference/sampled/` |
| figs. S5, S6 | `benchmarks/figures/figS5_S6_float.py` (`figures/supp/FigS5_certificate.svg`, `figures/supp/FigS6_float_errors.svg`) | `results/kimonu/float_comparison/parts/` | `benchmarks/kinetic_variants/float_exactcmp.sbatch`, `benchmarks/kinetic_variants/float_four.sbatch` |
| Tables S2, S3 | `benchmarks/figures/tableS2_S3.py` | `networks/unblocked/`, `networks/deficiency/`, `figures/source_data/figdata.json`, `results/kimonu/member_counts/` | `benchmarks/write_unblocked_crn.py`; `benchmarks/deficiency/deficiency_split.py`; `benchmarks/figures/gather_figdata.py`; `benchmarks/kinetic_variants/member_counts.sbatch` |
| Data S1 | `benchmarks/figures/dataS1_S2.py` (`figures/supp/DataS1_values.csv`) | `figures/source_data/figdata.json`, `figures/source_data/pair_counts.json`, `results/kimonu/uptake_modes/`, `results/kimonu/member_counts/`, `results/kimonu/timing/{tan,ref}`, `networks/unblocked/` | `benchmarks/figures/gather_figdata.py`; `benchmarks/figures/pair_sets.py`; `benchmarks/kinetic_variants/uptake_modes.sbatch`; `benchmarks/kinetic_variants/member_counts.sbatch`; `benchmarks/kinetic_variants/timing_exclusive.sbatch` (jobs 572312, 572311) |
| Data S2 | `benchmarks/figures/dataS1_S2.py` (`figures/supp/DataS2_uptake_acr.csv`) | `results/kimonu/acr_lists_by_uptake_mode/` | `benchmarks/kinetic_variants/uptake_lists.sbatch` (job 572114) |
| Text S1 | `benchmarks/figures/textS1_walkthrough.py` | none (Kimonu run on the example network) | none |
| numbers of the Results, Materials and Methods and Supplementary Text | `benchmarks/figures/paper_numbers.py`, one line per number with its source | all of the above, and `results/kimonu/flux_balance_vs_mass_action/`, `results/kimonu/gap_roots/`, `results/kimonu/second_order/`, `results/kimonu/gap_experiments/`, `results/kimonu/uniform_flux/`, `results/kimonu/genericity_bound/`, `results/kimonu/lp_point_rounding/`, `results/kimonu/preprocessing_time/`, `results/kimonu/one_lp_vs_fva/`, `results/hardware/` | see "Layout" |

`benchmarks/figures/fig4_figS1_S2_benchmark.py` draws the values of Figure 4 without the hand edits, and figs. S1 to S6
are published as their scripts draw them.

## Layout

Top level: `README.md` (this file), `LICENSE`, `CITATION.cff`, `SHA256SUMS`, `link_records.sh` (links the archive
into a checkout, see below) and `@SNAPSHOT@`. Every folder, what it holds, the paper items it feeds and the script or
SLURM job script that produced it (SLURM job ids as in the `provenance_*.txt` of the folder).

@LAYOUT@

Each Kimonu result folder of the reruns holds `provenance_*.txt` (date, node, SLURM job, Kimonu commit, StrainDesign
tree and version, Python), and each timing arm `hardware_*.txt` (hostname, lscpu, free -g).

## How the scripts find the records

`link_records.sh` links the records into a checkout (`eval "$(bash link_records.sh <checkout>)"`). It links each
folder where the scripts read it and keeps the files the checkout holds already (the networks, `figdata.json` and the
small Kimonu results are tracked in git, identical to the copies here). For the larger Kimonu runs the scripts read an
environment variable instead, and `link_records.sh` prints an `export` line for each.

@READING@

The variables, the archive folder `link_records.sh` sets each to, the run of the paper the scripts read without it,
and the scripts that read it.

@VARIABLES@

## Software and hardware

- Kimonu. The runs used commit 451a748 (`provenance_*.txt`). Later commits up to the snapshot @COMMIT@ restructure
  the package (`src/kimonu`), with new module and function names, without changing its results. Up to commit ed9240e
  they change it only in docstrings. After that, `benchmarks/kinetic_variants/refactor_equivalence.py` (run by
  `refactor_equivalence.sbatch`) ran concordance, the iterated and the tangent route with detached and attached
  uptakes, certify and kinetic flux coupling on all 33 networks at commit ed9240e and after the restructuring, and
  every result was identical (modules, ACR, ACRR, counts, K, M, B, the flux r* and the coupling groups). The latest
  stable version of Kimonu is under the releases of the repository,
  https://github.com/VonAlphaBisZulu/Kimonu.py/releases.
- StrainDesign 1.19.2 (git tag v1.19.2, commit 7143904), whose sparse RREF in rational arithmetic computes every
  nullspace.
- Python 3.12.12 with NumPy 2.4.3, SciPy 1.17.1 (HiGHS), SymPy 1.14.0, cobra 0.30.0, optlang 1.8.3, python-flint 0.9.0
  (FLINT 3.6.0), PySCIPOpt 6.2.1 with SCIP 10.0, CPLEX 22.2 and pandas 2.3.3. The figures were drawn with
  matplotlib 3.8.0 under Python 3.10, because matplotlib 3.10 changes the SVG font attributes that
  `benchmarks/figures/svg_finish.py` reads.
- COCOA.jl v1.1.0 (git tag v1.1.0, commit 26e9359) under Julia 1.12.6 with HiGHS and 64 workers
  (`benchmarks/julia_v110/` holds the Julia project).
- Reference implementation of Langary et al. (2025), https://github.com/ankueken/Upstream_Algorithm, run from the
  branch fix/negative-branch-charnes-cooper of the fork https://github.com/VonAlphaBisZulu/Upstream_Algorithm (pull
  request 1, head 8a6ed17) with the settings that needed adjusting at genome scale, described in Supplementary Text S3.
  Concordance stage under GNU Octave 10.3.0 with CPLEX 22.2 (bridge `benchmarks/matlab_reference/compat/lp_server.py`),
  kinetic stage in R 4.4.3.
- Hardware. Cluster nodes with two 16-core Intel Xeon Silver 4314 at 2.40 GHz (64 logical cores, 251 GB by free -g).
  Each timed Kimonu run was one process pinned to one core of an exclusive node (node01, node03, node04), with 1 GB
  reserved (3 GB for the tangent route). COCOA.jl ran 54 of its 183 runs on nodes with Intel Xeon Silver 4514Y and the
  others on the 4314 nodes, with 64 workers, and the reference ran on 2 cores of the 4314 nodes.
  `results/hardware/hardware_record.json` holds the node, CPU, memory and thread of every run and the versions above.

## Reproducing the figures, tables and numbers

    unzip Kimonu_records.zip && cd Kimonu_records && sha256sum -c SHA256SUMS
    unzip @SNAPSHOT@                          # or: git clone https://github.com/VonAlphaBisZulu/Kimonu.py && git -C Kimonu.py checkout @COMMIT@
    pip install "straindesign==1.19.2" cobra scipy sympy python-flint pandas "matplotlib==3.8.*"
    eval "$(bash link_records.sh Kimonu.py-@COMMIT@)"     # links the records into the checkout, sets the variables
    cd Kimonu.py-@COMMIT@/benchmarks
    python figures/gather_figdata.py; python figures/pair_sets.py; python figures/seed_spread.py   # figures/source_data
    python figures/fig4_figS1_S2_benchmark.py out/main   # Fig. 4 (out/main), figs. S1, S2 (out/supp)
    python figures/figS3_S4_uptakes_seeds.py out         # figs. S3, S4
    python figures/figS5_S6_float.py out                 # figs. S5, S6 (and the draft of Fig. 3E, F)
    python figures/fig3AB_hierarchy.py out               # Fig. 3A, B
    python figures/dataS1_S2.py out                      # Data S1, S2
    python figures/tableS2_S3.py out/TableS2.md out/TableS3.md
    python figures/paper_numbers.py                      # every number of the text, with its source

`benchmarks/zenodo_bundle.sh`, which built this archive, rebuilds the source data, figs. S1 to S6, Tables S2 and S3,
Data S1 and S2 and every number in this way from the zip alone and compares them with the files here. To rerun the
Kimonu computations themselves, follow "How to run" in `benchmarks/README.md`.

## License

Data, records and figures CC BY 4.0, code Apache 2.0 (see `LICENSE`).
