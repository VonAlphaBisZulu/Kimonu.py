# What running the reference implementation turned up

The reference (ankueken/Upstream_Algorithm) runs under Octave, licence-free, both on the shared
benchmark networks and on its own published elementary-step-split models.

## The Octave port is faithful

`compat/linprog.m` replaces MATLAB's Optimization Toolbox `linprog` with Octave's bundled GLPK.
Rerunning the reference's own `find_balanced_complexes` on its own published models reproduces the
committed balanced sets exactly:

| model | complexes | balanced, rerun | balanced, published | difference |
|---|---|---|---|---|
| iIS312_Trypomastigote | 1046 | 228 | 228 | none |
| iIS312_Amastigote | 1446 | 315 | 315 | none |
| iIS312_Epimastigote | 1738 | 379 | 379 | none |

Sets, not just counts, are identical in all three cases. Octave's `load` fails on some of the
committed result files with a zlib decompression error where scipy reads them without complaint, so
the published side is read in Python (`compare_balanced_to_published.py`) rather than in Octave.

## Result on the shared models

Same preprocessed network and same complex decomposition for all three implementations, scenario
(i), every candidate pair tested (no sampling prefilter).

| model | complexes | balanced (K/C/R) | modules (K/C/R) | concordant pairs (K/C/R) |
|---|---|---|---|---|
| EnvZ-OmpR | 16 | 5 / 5 / 5 | 3 / 3 / 3 | 10 / 10 / 10 |
| e_coli_core | 144 | 8 / 8 / 8 | 43 / 43 / 37 | 153 / 153 / 131 |
| e_coli_core anaerobic | 140 | 7 / 7 / 7 | 41 / 41 / 35 | 158 / 158 / 146 |

K = Kimonu, C = COCOA.jl, R = reference. Balanced complexes agree everywhere, and Kimonu and COCOA
agree exactly. The reference's concordant pairs are a **strict subset** on both e_coli models — it
misses 22 and 12 pairs and invents none.

## Why: the negative branch of the concordance test is infeasible by construction

`find_concordant_complexes.m` decides concordance with a Charnes–Cooper transform: set
`t = 1/(denominator activity)` and `w = t·v`. It runs two branches, `t ∈ [0, 999]` and
`t ∈ [-999, 0]`, both under the same constraint matrix

```matlab
a = sparse([eye(n) -model.ub; -eye(n) model.lb]);   %  a * [w; t] <= 0,  i.e.  lb*t <= w <= ub*t
```

Multiplying `lb <= v <= ub` by `t` only preserves the inequality directions while `t > 0`. For
`t < 0` the correct constraint is `ub*t <= w <= lb*t`, but the code keeps the `t > 0` form. The
resulting interval is empty for every reaction with `lb < ub`, and `t = 0` forces `w = 0`, which
fails the normalisation `A(denominator,:)·w = 1`. **The negative branch is therefore infeasible for
any model whose reactions have a non-degenerate flux interval** — no LP solver enters into it.

Verified on the reference's own published input: `iIS312_Trypomastigote_pre_balanced_fixed.mat`
(1046 complexes, 1339 reactions) has `lb = 0`, `ub = 1000` for all reactions and `lb == ub` for
none of them; the negative branch is infeasible for 8 of 8 tested denominators.

### Confirmed by correcting it

Swapping `lb` and `ub` in the constraint matrix for the negative branch only — the one-line change
the derivation implies — recovers every missed pair:

| model | pairs missed | recovered with the corrected branch | still missed |
|---|---|---|---|
| e_coli_core | 22 | 22 | 0 |
| e_coli_core anaerobic | 12 | 12 | 0 |

The negative branch, infeasible in 0 of 22 and 0 of 12 cases before, is feasible in all of them
after. (`test_corrected_branch.m`.)

### How it reaches the answer

The defect produces false negatives by two routes, and an audit of all 22 missed pairs on
e_coli_core accounts for every one of them:

1. **Denominator with strictly negative activity (`group == 'N'`)** — only the negative branch is
   run, so the pair is always rejected.
2. **Sign-variable denominator (`group == 'U'`)** — both branches run and the update rule requires
   them to agree:
   ```matlab
   CC(i,j) = (round(Max_p,2) == round(Min_p,2) & round(Max_n,2) == round(Min_n,2) & ...
              round(Max_p,2) == round(Min_n,2)) * 2;
   ```
   The positive branch proves a constant ratio, but the infeasible negative branch leaves
   `Max_n = Inf`, `Min_n = -Inf`, so the conjunction fails.

Not one missed pair was rejected because the ratio bounds genuinely differed.

## Which pairs are affected is decided by complex numbering

Candidate pairs come from `nchoosek(1:n, 2)`, so only `i < j` is tested and the higher-numbered
complex is always the denominator. Concordance is symmetric; this test is not. For the pair
(34, 112) in e_coli_core the reference returns *not concordant* as tested, and *concordant with
ratio -1* when the same pair is presented as (112, 34). Renumbering the complexes changes the
answer.

## Exposure at genome scale

Any complex classified `N` or `U` can only be tested as a denominator through the branch that is
infeasible. On the reference's own models that is the majority of complexes:

| model | complexes | P | N | U | balanced | exposed (N or U) |
|---|---|---|---|---|---|---|
| iIS312_Trypomastigote | 1046 | 244 | 227 | 347 | 228 | 574 (54.9%) |
| iIS312_Amastigote | 1446 | 253 | 240 | 638 | 315 | 878 (60.7%) |
| iIS312_Epimastigote | 1738 | 328 | 313 | 718 | 379 | 1031 (59.3%) |

## Against the published numbers

Kimonu analysed `iIS312_Trypomastigote` (ordered binding) directly from the reference's own
pre-split model and compared against the committed result:

| | complexes | balanced | modules | concordant pairs | seconds |
|---|---|---|---|---|---|
| Kimonu | 1046 | 228 | 111 | 5405 | 0.12 |
| published reference | 1046 | 228 | 307 | 429 | — |

Balanced complexes match exactly. The published concordant pairs are a strict subset: 4976 of 5405
are missing. Auditing every one of those 4976 with the reference's own LP test separates the two
mechanisms:

| cause | pairs | share |
|---|---|---|
| `N` denominator — negative branch infeasible | 1924 | 38.7% |
| `U` denominator — positive branch proves it, negative branch infeasible | 699 | 14.0% |
| `P` denominator — the LP accepts the pair, so only the sampling prefilter can explain the miss | 2353 | 47.3% |
| ratio bounds genuinely differ | 0 | 0% |

So roughly half of the gap is the negative-branch defect and roughly half is the stochastic
prefilter discarding candidates before any LP runs. Neither is a numerical-tolerance effect: no
pair was rejected because its activity ratio actually varies.

## Caveats

- Run under Octave with GLPK rather than MATLAB with `linprog`. The mechanism is structural — an
  empty constraint set — so it is not solver-specific, and the correction test reproduces the
  expected pairs exactly. It has not been confirmed against MATLAB itself.
- Only the concordance stage was exercised. The kinetic-module and ACR/ACRR stages are R code and
  were not run; nothing here speaks to them, and the published giant-module and ACR figures are
  downstream of the concordance partition rather than equal to it.
- The published-number comparison is one model and one binding variant.
- COCOA.jl does **not** reproduce this defect: it agrees with Kimonu exactly on all four shared
  models. On this point the comparison runs in COCOA's favour.
