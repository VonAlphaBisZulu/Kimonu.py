# Run records and figure data for "Exact nullspaces find and certify kinetic modules and other mass-action invariants in metabolic networks"

Philipp C. Schneider, Department of Molecular Biology, Massachusetts General Hospital, Boston, MA, USA
(ORCID 0000-0001-8858-428X).
Article DOI @ARTICLE_DOI@, archive DOI @ZENODO_DOI@, version @ARCHIVE_VERSION@.
Built on @DATE@ by `benchmarks/zenodo_bundle.sh` and `make_zenodo_bundle.py` of Kimonu.py at commit @COMMIT@.

## What this archive holds

This archive holds every file that the figures, tables, Data S1 and S2 and the numbers of the paper are computed from.
These are the 33
input networks, the Kimonu runs of 2026-11-11 on the final code, the runs of COCOA.jl v1.1.0 and of the reference
implementation of Langary et al. (2025), which are frozen records of August to November 2026, the figure source data,
the published figures with Data S1 and S2, and a snapshot of the source code (`@SNAPSHOT@`). The scripts that made each
item are in the repository https://github.com/VonAlphaBisZulu/Kimonu.py, and `benchmarks/README.md` there maps every
paper item to its script, inputs and outputs. `SHA256SUMS` lists every file with its SHA-256 checksum.

Superseded runs are left out. Two kinds of files that the scripts read are in the source snapshot only, where git
tracks them. These are the list of the 33 networks (`benchmarks/kinetic_variants/models33.txt`) and the upstream-set
arm of September, which `gather_figdata.py` still reads into `figdata.json` but no figure, table or number shows
(`sweep_shared/`, `sweep_tentative/`, `kinetic_variants/s210_*.json`). The records of the Supplement experiments that
the text states in words, without a number, are also in the snapshot only (`kinetic_variants/two_point_*.json` and
similar, listed in `benchmarks/README.md`).

## Layout

Top level: `README.md` (this file), `LICENSE`, `CITATION.cff`, `SHA256SUMS`, `link_records.sh` (links the archive
into a checkout, see below) and `@SNAPSHOT@`.

@PARTS@

`figures/` holds the published figures (`Figure1.svg` to `Figure4.svg`, `supp/FigS_*.svg`) and Data S1 and S2
(`supp/DataS1_values.csv`, `supp/DataS2_uptake_acr.csv`).

@COUNTS@

Each Kimonu result folder holds `provenance_<job>.txt` (date, node, SLURM job, Kimonu commit, StrainDesign tree and
version, Python), and each timing arm `hardware_<job>.txt` (hostname, lscpu, free -g).

## Paper item by item

| paper item | files of this archive | script (in `benchmarks/`) |
|---|---|---|
| Fig. 1, 2 | `figures/Figure1.svg`, `Figure2.svg` (closed-form examples, no data) | edit scripts in `docs/examples/network/` |
| Fig. 3A, B | `results/kimonu/hierarchy_detached/` (from `uptake_modes/` and `member_counts/`) | `figures/make_hierarchy.py`, `figure3_ratios.py` |
| Fig. 3C, D | `results/kimonu/member_counts/`, `networks/unblocked/` | `figures/extend_figure3.py` |
| Fig. 3E, F | `results/kimonu/float_comparison/parts/` | `figures/make_figure2_bottom.py` |
| Fig. 4A | `figures/source_data/figdata.json` (from `results/kimonu/timing/{det,att,ref}`, `results/cocoa_jl/`, `results/reference/`) | `figures/gather_figdata.py`, `make_figure3.py` |
| Fig. 4B | `results/kimonu/timing/{tan,small}`, `results/kimonu/symbolic_baselines/` | `figures/make_figure3.py` |
| Fig. 4C, D | `figures/source_data/pair_counts.json`, `pair_sets.json` (from `results/kimonu/timing/ref` and the module files of COCOA.jl and the reference) | `figures/pair_sets.py`, `make_figure3.py` |
| Fig. 4E, F | `figdata.json`, `results/reference/exhaustive/reference_iLJ478.json`, `results/cocoa_jl/efficient_true/cell_iLJ478_seed1234.log`, `results/*/kinetic_stage_iLJ478/` | `figures/make_fig3_methods.py` |
| fig. S1 | `results/kimonu/timing/ref`, `results/reference/exhaustive/acr/*.log`, COCOA.jl logs | `figures/make_figure3.py` |
| fig. S2 | `figdata.json` | `figures/make_figure3.py` |
| fig. S3 | `figdata.json` | `figures/make_supp_figures.py` |
| fig. S4 | `figures/source_data/seed_spread.json` (from `results/kimonu/timing/det`, `results/cocoa_jl/`, `results/reference/sampled/`) | `figures/seed_spread.py`, `make_supp_figures.py` |
| figs. S5, S6 | `results/kimonu/float_comparison/parts/` | `figures/make_float_four.py` |
| Tables S2, S3 | `networks/unblocked/`, `networks/deficiency/`, `figdata.json`, `results/kimonu/member_counts/` | `figures/make_table_s2.py` |
| Data S1 | `figdata.json`, `pair_counts.json`, result `figures/supp/DataS1_values.csv` | `figures/make_data_s1_s2.py` |
| Data S2 | `results/kimonu/uptake_lists/`, result `figures/supp/DataS2_uptake_acr.csv` | `figures/make_data_s1_s2.py` |
| numbers of the Results, Materials and Methods and Supplementary Text | all of the above, also `fb_vs_ma/`, `gap_roots/`, `second_order/`, `gap_experiments/`, `uniform_flux/`, `genericity_bound/`, `scale_check/`, `preprocessing_time/`, `one_lp/` | `figures/paper_numbers.py` (one line per number with its source) |

The published Figures 1 to 4 carry hand edits made in Adobe Illustrator. Figure 3B to F and Figure 4 were updated in
place on those files by `figures/figure3_ratios.py`, `extend_figure3.py`, `make_figure2_bottom.py`,
`figure3_legend.py` and `figure4_update.py`. `figures/make_figure3.py` draws the same values of Figure 4 without the
hand edits.

## Software and hardware

- Kimonu. The runs used commit 451a748 (`provenance_*.txt`). Later commits up to the snapshot @COMMIT@ change the
  package (`src/kimonu`) only in docstrings. Release tag @RELEASE_TAG@.
- StrainDesign 1.19.2 (git tag v1.19.2, commit 7143904), whose sparse RREF in rational arithmetic computes every
  nullspace.
- Python 3.12.12 with NumPy 2.4.3, SciPy 1.17.1 (HiGHS), SymPy 1.14.0, cobra 0.30.0, optlang 1.8.3, python-flint 0.9.0
  (FLINT 3.6.0), PySCIPOpt 6.2.1 with SCIP 10.0, CPLEX 22.2 and pandas 2.3.3. The figures were drawn with
  matplotlib 3.8.0 under Python 3.10, because matplotlib 3.10 changes the SVG font attributes that
  `figures/svg_finish.py` reads.
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
    python figures/make_figure3.py out        # Fig. 4, figs. S1, S2
    python figures/make_supp_figures.py out   # figs. S3, S4
    python figures/make_float_four.py out     # figs. S5, S6
    python figures/make_hierarchy.py out      # Fig. 3A, B
    python figures/make_data_s1_s2.py out     # Data S1, S2
    python figures/make_table_s2.py out/TableS2.md out/TableS3.md
    python figures/paper_numbers.py           # every number of the text, with its source

`link_records.sh` links the COCOA.jl and reference records into `benchmarks/sweep/` and `benchmarks/reference_final*/`,
where the scripts read them, and keeps the files the checkout holds already (the networks, `figdata.json` and the
small Kimonu results are tracked in git, identical to the copies here). It prints `export` lines for the variables
that point the scripts at the rest (`KIMONU_TIMING`, `KIMONU_DET_DIR`, `KIMONU_ATT_DIR`, `KIMONU_REF_DIR`,
`KIMONU_LACR`, `KF_OUT`, `GB_OUT`, `GAP_OUT`, `SZ_DIR`, `SCALE_DIR`, `FVA_OUT`, `UPTAKE_LISTS`).
`benchmarks/zenodo_bundle.sh`, which built this archive, rebuilds the source data, tables S2 and S3, Data S1 and S2
and every number in this way from the zip alone and compares them with the files here. To rerun the Kimonu
computations themselves, follow "How to run" in `benchmarks/README.md`.

## License

Data, records and figures CC BY 4.0, code Apache 2.0 (see `LICENSE`).
