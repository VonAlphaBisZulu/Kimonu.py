# Records behind "Exact nullspaces find and certify kinetic modules and other mass-action invariants in metabolic networks"

The code is at https://github.com/VonAlphaBisZulu/Kimonu.py. This archive holds the records. `SHA256SUMS` lists every file.

| Path | Content |
|---|---|
| `sweep/*_fixed.crn.json` | the 33 networks before blocked-reaction removal |
| `unblocked/` | the same networks with blocked reactions removed, the input of every tool |
| `results/kimonu_timing_seeds/det`, `att` | Kimonu's default mode, detached and attached uptakes, seeds 0, 1, 2, with modules, ACR and ACRR |
| `results/kimonu_timing_seeds/tan`, `tanatt`, `small` | tangent route, on the 33 networks and on the small networks and size ladder |
| `results/kimonu_timing_final/ref` | Kimonu's reference-equivalent mode, three runs |
| `results/kimonu_timing_seeds/fva` | one LP against FVA, and preprocessing times |
| `sweep/cocoa_v110`, `sweep/cocoa_eff` | COCOA.jl v1.1.0, efficient=false and efficient=true |
| `reference_final`, `reference_final_sampled` | the reference implementation, exhaustive and sampled |
| `results/kimonu_float_four/parts` | floating-point comparison, per network and tolerance (Figure 2H, I) |
| `results/kimonu_local_acr` | symbolic local-ACR criteria (Figure 3B) |
| `results/kimonu_cert_modp`, `results/kimonu_sz_flux` | genericity certificate and uniform random flux |
| `results/kimonu_second_order`, `results/kimonu_twopoint`, `kinetic_variants/injectivity_results` | second-order test, two-point criterion, sign condition on the gap networks |
| `kinetic_variants/hierarchy` | both routes per network (Figure 2F, G) |
| `figdata.json` | the values the figures show |

To rebuild the figures, place the `benchmarks/figures` scripts next to these files and run from this directory:

    R=$PWD/results
    export KIMONU_DET_DIR=$R/kimonu_timing_seeds/det KIMONU_REF_DIR=$R/kimonu_timing_final/ref KIMONU_ATT_DIR=$R/kimonu_timing_seeds/att
    export KIMONU_TIMING=$R/kimonu_timing_seeds KIMONU_LACR=$R/kimonu_local_acr KF_OUT=$R/kimonu_float_four
    python figures/gather_figdata.py
    python figures/make_figure3.py out; python figures/make_fig3_methods.py out; python figures/make_float_four.py out
    python figures/make_supp_figures.py out; python figures/make_hierarchy.py out
