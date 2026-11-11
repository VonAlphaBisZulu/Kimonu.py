# Records behind "Exact nullspaces find and certify kinetic modules and other mass-action invariants in metabolic networks"

The code is at https://github.com/VonAlphaBisZulu/Kimonu.py, and `benchmarks/README.md` there maps every paper item to
the script that makes it. This archive holds the records. `SHA256SUMS` lists every file. Paths under `results/` are
the run directories of the cluster scratch space, all other paths keep their place under `benchmarks/`.

| Path | Content | Paper |
|---|---|---|
| `sweep/*_fixed.crn.json` | the 33 networks before blocked-reaction removal | Methods |
| `unblocked/` | the same networks with blocked reactions removed, the input of every tool | Table S2 |
| `results/kimonu_timing_seeds/det`, `att` | Kimonu's default (iterated route), detached and attached uptakes, seeds 0, 1, 2, with modules, ACR and ACRR | Fig. 4A, figs. S3, S4, Data S1 |
| `results/kimonu_timing_seeds/tan`, `tanatt`, `small` | tangent route on the 33 networks (detached, attached) and on the small networks and the size ladder | Fig. 4B |
| `results/kimonu_timing_final/ref` | Kimonu's reference-equivalent mode, three runs | Fig. 4A, figs. S1, S2 |
| `results/kimonu_trip_2026-11-02/refmode_modules` | the modules of the reference-equivalent mode | Fig. 4C, D |
| `results/kimonu_timing_seeds/fva` | one LP against FVA, preprocessing and FVA times | Methods |
| `sweep/cocoa_v110`, `sweep/cocoa_eff` | COCOA.jl v1.1.0, efficient=false and efficient=true | Fig. 4, figs. S1, S4, Data S1 |
| `reference_final`, `reference_final_sampled` | the reference implementation, exhaustive and sampled | Fig. 4, figs. S1, S2, S4, Data S1 |
| `kinetic_variants/uptake_modes`, `kinetic_variants/hierarchy_detached` | both routes and both uptake modes per network | Fig. 3A, B, Results |
| `kinetic_variants/member_counts` | complexes in kinetic modules and metabolites with ACR or ACRR | Fig. 3C, D |
| `kinetic_variants/fb_vs_ma` | flux coupling under flux balance and under mass action | Results |
| `results/kimonu_uptake_modes/lists` | ACR and ACRR lists of the iterated route per uptake mode | Data S2, Supplementary Text |
| `results/kimonu_gaproots` | the open questions the gap networks come down to | Results |
| `results/kimonu_float_four/parts` | floating-point comparison, per network and tolerance | Fig. 3E, F, figs. S5, S6 |
| `results/kimonu_local_acr` | symbolic local-ACR criteria | Fig. 4B |
| `results/kimonu_seeds` | Kimonu under three seeds of the reference flux | fig. S4 |
| `results/kimonu_cert_modp`, `results/kimonu_sz_flux`, `results/kimonu_sz_flux_33` | second evaluation modulo a prime, and the uniformly drawn flux | Methods |
| `results/kimonu_second_order`, `results/kimonu_twopoint`, `kinetic_variants/injectivity_results` | second-order test, two-point criterion, sign condition on the gap networks | Methods, Supplementary Text |
| `figdata.json`, `figures/pair_counts.json`, `figures/pair_sets.json`, `figures/seed_spread.json` | the values the figures show | all |

To rebuild the figures and tables, place the `benchmarks/` scripts of the repository next to these files and run from
this directory:

    R=$PWD/results
    export KIMONU_DET_DIR=$R/kimonu_timing_seeds/det KIMONU_ATT_DIR=$R/kimonu_timing_seeds/att KIMONU_REF_DIR=$R/kimonu_timing_final/ref
    export KIMONU_TIMING=$R/kimonu_timing_seeds KIMONU_LACR=$R/kimonu_local_acr KF_OUT=$R/kimonu_float_four
    export KIMONU_REFMOD_DIR=$R/kimonu_trip_2026-11-02/refmode_modules SEED_OUT=$R/kimonu_seeds UPTAKE_LISTS=$R/kimonu_uptake_modes/lists
    python figures/gather_figdata.py; python figures/pair_sets.py; python figures/seed_spread.py
    python figures/make_figure3.py out        # Figure 4 and figs. S1, S2
    python figures/make_supp_figures.py out   # figs. S3, S4
    python figures/make_float_four.py out     # figs. S5, S6
    python figures/make_hierarchy.py out      # Figure 3A, B
    python figures/make_data_s1_s2.py out; python figures/make_table_s2.py out/TableS2.md
    python figures/paper_numbers.py           # every number of the Results and Methods, with its source

Figure 3C to F were drawn into the Illustrator figure by `figures/extend_figure3.py` and `figures/make_figure2_bottom.py`
(after `figures/split_figure2.py`). Figures 1 and 2 are drawn by hand, and the scripts in `docs/examples/network/`
record each edit.
