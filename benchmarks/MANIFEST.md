# Benchmark manifest — where everything is, and its state (2026-09-17)

Start here. `FINDINGS.md` is the running record with its retractions; this file is the map.
Snapshot of every directory below, checksummed: `/scratch/sadreyev_lab/pschneider/kimonu_benchmark_snapshots/benchmarks_2026-09-17.tar.gz` (sha256 `88acd39b…b55a47b`, 19 MB).
Live pages: ledger https://claude.ai/artifact/RGt8asxhRiyDcJCpmVKrrF · figures https://claude.ai/artifact/JvoAa9HHVaCUoJmM2VyjNv (sources in `docs/`).

## Inputs

| path | what | tracked |
|---|---|---|
| `unblocked/<m>_fixed.crn.json` | the 33 networks, blocked reactions removed. **Every arm reads these.** | yes |
| `sweep/<m>_fixed.crn.json` | same networks before blocked-reaction removal; iJB785 and iSynCJ816 still carry 363 / 280. Do not analyse from here without the FVA. | yes |
| `coupling_exact/` | full flux-coupling matrices from the exact nullspace for the 6 sound models (F2C2-shaped, code 1 only) | yes |

## Kimonu

| path | what | state |
|---|---|---|
| `sweep_shared/kimonu_<m>.json` | **the figure source**: Kimonu on `unblocked/`, `--assume-unblocked`, 3 timing repeats | 33/33, identical results to `sweep/` |
| `sweep/kimonu_<m>_fixed.json` | Kimonu on `sweep/` with its own FVA (older) | 33/33, superseded by `sweep_shared/` |
| `sweep_tentative/<m>.log` | the §S.2.2 four-phase algorithm run inside Kimonu (`matlab_reference/tentative_upstream.py`) | 33/33; giant agrees with S2-1 on 25 |

## COCOA v1.1.0

| path | what | state |
|---|---|---|
| `sweep/cocoa_v110/` | `efficient=false`, seeds 1234/2026/31337 | 30/33 at 3 seeds; iAF1260b, iSDY_1059, iS_1188 queued | 
| `sweep/cocoa_eff/` | `efficient=true` | 14/33; array 547571 running (20 left, 3 concurrent) |
| `sweep/cocoa/` | v1.0.0 — **stale, never report** | 31/33, kept only as history |
| `matlab_reference/cocoa_isolation/` | COCOA's kinetic stage fed Kimonu's exact partition → 250/19/591/8 = Kimonu | done, iCN900 |

COCOA JSON `complex_ids` write coefficients as `2_h[c]` / `4.5_h[c]`; `concordance_module == 0` is balanced.

## Reference (Upstream_Algorithm under Octave + CPLEX bridge)

All runs before commit `3ea2bf7` used an order-dependent bridge and are suspect at marginal pairs. Balanced complexes are unaffected (no pair LP) and match Kimonu 9/9.

| path | what | algorithm / bridge | state |
|---|---|---|---|
| `reference_full/` | exhaustive, run to completion, 5 mid-size models | unfixed / old | ground truth for balanced; concordance sound on iCN900, adjudicated on iJB785 |
| `reference/` | exhaustive, 24 h cap; 4 iIS312 models complete (tracked), 29 truncated (ignored) | unfixed / old | the 4 complete are ground truth for balanced |
| `reference_sampled/` (seed 1234) + `reference_sampled_seeds/` (2026, 31337) | sampled mode, 9 models | unfixed / old | 8 at 3 seeds, iLJ478 at 1 |
| `reference_fixed/` | 3 concordance fixes, old bridge — Trypomastigote, Amastigote | fixed / old | both reproduce Kimonu exactly |
| `reference_fixed2/` | 3 fixes + fixed bridge — **the gate** | fixed / fixed | both models exact vs Kimonu, 0 degraded LPs — passed |
| `reference_acr/` | corrected R kinetic + ACR stage on the 6 sound models, with README on reading the tables | corrected R | done 6/6; every ACRR pair ⊆ Kimonu's |
| `reference_final/`, `reference_final_sampled/` | **the final reruns** — jobs 549287–549295 + 549302/549303 (11 exhaustive; 14-day cap, 21 for iNF517/iSB619) and arrays 549296 + 549304 (33 sampled cells, 7-day cap; 13 cells OOM at 4 GB and resubmitted per model as 549880-549884, see FINDINGS) | fixed / fixed | launched 2026-09-18 after both gate models matched Kimonu exactly |

Code: `matlab_reference/run_reference.m` (harness; now saves metabolite names), `compat/` (CPLEX bridge, two-arg `round`), `test_determinism.m`, `test_fixed_verdict.m`, `patch_mat_names.m`, `write_coupling_exact.py`, `compare_kinetic_partitions.py`, `PR_DRAFT.md`.
Upstream: fork `VonAlphaBisZulu/Upstream_Algorithm`, branch `fix/negative-branch-charnes-cooper`, 4 commits, pushed, no PR. `~/pipelines/2026-08_Kimonu/Upstream_Algorithm_benchmark` is pinned at 088d827 for reproducing old runs only.

## Figures

`figdata.json` (all arms, built by the validated parsers) → `figures/make_prelim.py` → `figures/prelim/F1–F5.{png,svg}`. Median over seeds with seed points. Reference columns are pre-rerun and marked as such.

## In flight

- 549287–549295, 549302/549303 + sampled arrays 549296, 549880–549884 final reference reruns on node01 · 547571 COCOA `efficient=true` · 548401/548402 COCOA gap cells · 547572 iAF1260b `efficient=false` · iNF517/iSB619 now in the final rerun (old-bridge jobs cancelled).
- Node reservation: node05–07 and node12–16 excluded throughout.

## The handful to bring home

1. ~~Gate~~ passed; reruns launched 2026-09-18.
2. Reruns land (~6 days, by ~24 Sep) → regenerate `figdata.json`, re-render, reference columns become final.
3. COCOA `efficient=true` completes on its own.
4. Methods text: Kimonu implements Remark S2-1 (two phases, full-graph terminality); the reference implements the §S.2.2 tentative version; Remark S2-2's equivalence claim has a counterexample.
5. PR when told.
