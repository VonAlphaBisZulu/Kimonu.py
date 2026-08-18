# What running the reference implementation turned up

Status: the reference (ankueken/Upstream_Algorithm) runs under Octave, licence-free, on both the
shared benchmark networks and its own published pre-split models.

## Result on the shared models

Same preprocessed network and same complex decomposition for all three implementations, scenario
(i), every candidate pair tested (no sampling prefilter).

| model | complexes | balanced (K/C/R) | modules (K/C/R) | concordant pairs (K/C/R) |
|---|---|---|---|---|
| EnvZ-OmpR | 16 | 5 / 5 / 5 | 3 / 3 / 3 | 10 / 10 / 10 |
| e_coli_core | 144 | 8 / 8 / 8 | 43 / 43 / 37 | 153 / 153 / 131 |
| e_coli_core anaerobic | 140 | 7 / 7 / 7 | 41 / 41 / 35 | 158 / 158 / 146 |

K = Kimonu, C = COCOA.jl, R = reference. Balanced complexes agree everywhere. Kimonu and COCOA
agree exactly everywhere. The reference's concordant pairs are a **strict subset** on both e_coli
models — it misses 22 and 12 pairs respectively, and invents none.

## Why: the negative branch of the concordance test is infeasible by construction

`find_concordant_complexes.m` decides concordance with a Charnes–Cooper transform. It runs two
branches, one with the homogenising variable `t ∈ [0, 999]` and one with `t ∈ [-999, 0]`, under the
constraint matrix

```matlab
a = sparse([eye(n) -model.ub; -eye(n) model.lb]);   %  a * [v; t] <= 0
```

which is `lb*t <= v <= ub*t`. For `t < 0` that interval is empty for every reaction with
`lb < ub`: the constraint requires `(lb - ub) * t <= 0`, and with `lb <= ub` and `t < 0` the
product is non-negative. `t = 0` forces `v = 0`, which then fails the normalisation
`A(denominator,:) * v = 1`. **The negative branch is therefore infeasible for any model whose
reactions have a non-degenerate flux interval** — no LP solver enters into it.

Verified on the reference's own published input, not only on ours:
`iIS312_Trypomastigote_pre_balanced_fixed.mat` (1046 complexes, 1339 reactions) has `lb = 0`,
`ub = 1000` for all 1339 reactions, `lb == ub` for none of them, and the negative branch is
infeasible for 8 of 8 tested denominators.

That single defect produces the missed pairs by two routes, and an audit of all 22 missed pairs on
e_coli_core accounts for every one of them:

1. **Denominator with strictly negative activity (`group == 'N'`)** — 7 pairs. Only the negative
   branch is run for these, so the pair is always rejected.
2. **Sign-variable denominator (`group == 'U'`)** — 15 pairs. Both branches run, and the update
   rule requires the positive *and* negative branch to agree:
   ```matlab
   CC(i,j) = (round(Max_p,2) == round(Min_p,2) & round(Max_n,2) == round(Min_n,2) & ...
              round(Max_p,2) == round(Min_n,2)) * 2;
   ```
   The positive branch proves a constant ratio (max = min = ±1 exactly in all 15 cases), but the
   infeasible negative branch leaves `Max_n = Inf`, `Min_n = -Inf`, so the conjunction fails.

Not one of the 22 missed pairs was rejected because the ratio bounds genuinely differed.

## Which pairs are affected is decided by complex numbering

Candidate pairs come from `nchoosek(1:n, 2)`, so only `i < j` is ever tested and the
higher-numbered complex is always the denominator. Concordance is symmetric, but this test is not:
for the pair (34, 112) in e_coli_core the reference returns *not concordant* as tested, and
*concordant with ratio -1* when the same pair is presented as (112, 34). Renumbering the complexes
changes the answer.

## Caveats

- Run under Octave with GLPK rather than MATLAB with `linprog`. The mechanism above is structural
  — an empty constraint set — so it is not solver-specific, and the failures were 100% consistent
  across every negative-branch LP attempted. It has not been confirmed against MATLAB itself.
- Only the concordance stage was exercised. The kinetic-module and ACR/ACRR stages are R code and
  were not run; nothing here speaks to them.
- Three small models. Whether the same defect shifts the published genome-scale numbers is a
  separate question — it depends on how many pairs have `N` or `U` denominators in those networks.
- COCOA.jl does **not** reproduce this defect: it agrees with Kimonu exactly on all four shared
  models. Whatever disagreements exist between COCOA and the reference in the COCOA paper's
  supplement, this one runs in COCOA's favour.
