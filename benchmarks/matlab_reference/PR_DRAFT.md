# Title

Fix three defects in the pairwise concordance test

# Body

We have been using `Upstream_Algorithm` as the reference implementation while validating an
independent exact method for balanced complexes and concordance modules, and running it on the 33
networks from the paper turned up three problems in `find_concordant_complexes`. Each is a small
change; together they alter the concordance partition substantially on some models. All three are
in the same function and interact, so they are in one PR.

Happy to split them, add tests elsewhere, or adjust anything to your preferences.

## 1. The negative Charnes-Cooper branch uses the positive branch's flux bounds

With `w = t*v`, `lb <= v <= ub` becomes `lb*t <= w <= ub*t` only while `t > 0`; a negative `t`
reverses both inequalities. Both branches used

```matlab
a = sparse([eye(size(model.S,2)) -model.ub; -eye(size(model.S,2)) model.lb]);
```

so in the `t < 0` branch the interval is empty for every reaction with `lb < ub`, and `t = 0`
forces `w = 0`, which violates the normalisation. That branch is therefore infeasible for any
model with non-degenerate flux intervals, independently of the LP solver.

Two kinds of pair are affected. A denominator with strictly negative activity (`group` `'N'`) is
only ever tested through that branch, so it is always reported non-concordant. A sign-variable
denominator is tested through both and the update requires them to agree, so the infeasible branch
overrules the positive one. Because candidate pairs come from `nchoosek`, the higher-numbered
complex is always the denominator, and the verdict depends on the numbering: in `e_coli_core` the
pair `(34, 112)` is non-concordant while `(112, 34)` is concordant with ratio -1.

## 2. The scaling variable is in the objective

Under the same substitution, maximising `(A_i.v + alpha)/(A_j.v)` subject to `A_j.w = 1` becomes
maximising `A_i.w + alpha*t`. The numerator here is `A_i.v` with no constant term, so `alpha = 0`.
The objective vector is `[-model.A(At(i,1),:) 1]`, so the solver trades the ratio off against `t`
and stops at a point that is not the extremum of the ratio that is then computed from the solution
and compared.

The reported range collapses towards zero. On `iLJ478` with blocked reactions removed:

| pair | objective `[-A_i 1]` | objective `[-A_i 0]` |
|---|---|---|
| (1842, 2368) | max 0, min 0 | max 0, min -7.299e-06 |
| (2244, 1625) | max 4.912e-07, min 0 | max 6.090e-01, min 0 |
| (2403, 3295) | max 0, min -1.287e-04 | max 0, min -1.840e+02 |

The first pair is reported as having a constant ratio and accepted as concordant. Its ratio is not
constant. Since the classes are the connected components of this relation, each such acceptance
merges two classes outright.

## 3. The verdict compares the bounds on an absolute grid, and lets an unsolved branch veto

The bounds were compared as `round(Maximum_c_p,2) == round(Minimum_c_p,2)`. That is an absolute
tolerance of 5e-3 on the ratio, so it is loosest around zero: a pair whose ratio ranges over
(-0.005, 0.005) without being constant still lands in one bin. Combined with defect 2, which puts
the reported bounds near zero, this is what lets the false positives through. Replaced with a
relative comparison.

Separately, the sign-variable case requires both branches to agree. An infeasible branch leaves
`Inf` and `-Inf` and fails the conjunction however cleanly the other branch behaved. But the two
branches cover the parts of the flux space where the denominator is positive and negative, so an
infeasible branch means the denominator never takes that sign and carries no information. The
verdict now uses the branches that were solved: constant wherever defined, and equal across both
when both exist.

The per-pair results are also reset at the top of the loop, since only the branch matching the
denominator's sign group is evaluated.

## Effect

On `iLJ478` (3439 complexes, blocked reactions removed, every candidate pair tested), the
published code gives 91 concordance classes with a largest class of 2194 complexes. Checking that
class against the definition directly -- sampling the flux polytope and measuring activity ratios,
independently of any implementation -- its members span 19 distinct proportional directions, where
a concordance class must span one. 188,877 of the 482,007 accepted pairs join complexes that an
exact method separates.

In the other direction, one 12-complex class is returned as six pairs. Those complexes have a
ratio of exactly +-1 to each other, constant to a spread of 6e-11 over 64 flux points.

Balanced complexes are unaffected by all three changes; they agree exactly with two independent
implementations on all nine models we have complete runs for.

## Validation

- `e_coli_core`: the partition agrees exactly with two independent implementations (8 balanced,
  43 modules, 153 concordant pairs, up from 37 and 131), and likewise on an anaerobic variant
  (7, 41, 158, up from 35 and 146) and on the EnvZ-OmpR model, where all three already agreed.
- Eight pairs on `iLJ478` whose status was established from the definition with flux vectors, four
  concordant and four not, are classified correctly only with all three changes applied.

Run under Octave, with `linprog` backed by CPLEX and a two-argument `round` supplied, since Octave
has neither. The mechanisms are structural rather than numerical, so they are not specific to that
setup, but the changes have not been exercised under MATLAB itself.
