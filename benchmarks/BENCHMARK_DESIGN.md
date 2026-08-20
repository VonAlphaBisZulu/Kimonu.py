# Three-way benchmark: design decisions

## Input: one shared, pre-unblocked network per model

Every tool is handed the identical CRN from `benchmarks/unblocked/`, produced once by
`write_unblocked_crn.py`. Blocked reactions are removed there, not inside any arm.

Removing them is not optional. `span{v : Nv = 0, v >= 0}` equals `ker(N)` only when no reaction is
blocked; otherwise Kimonu's kernel test asks about a strictly larger space than the feasible set
and under-reports. Measured on iSynCJ816: 2991 balanced raw against 3770 unblocked, and 3770 is
exactly the published number.

Costing it once rather than inside one arm matters because it is the same FVA for everyone --
charging it to Kimonu while COCOA and the reference receive it free would misattribute it. It is
1-159 s per model with SCIP.

Only **2 of 33 models have any blocked reactions**: iSynCJ816 (280) and iJB785 (363). For the
other 31 the unblocked CRN is a copy of the input.

## What is deliberately NOT applied

**COCOA's wider preprocessing pipeline.** Their README recommends
`normalize_bounds -> remove_blocked_reactions -> remove_orphans -> split_into_elementary ->
split_into_irreversible`. We apply only the blocked-reaction step, because the reference's models
arrive already elementary-split and irreversible, and `normalize_bounds` is a **no-op** on them:
all 33 inputs already have `lb = 0, ub = 1000` on every reaction, so it would change zero entries.

TODO -- worth doing at some point, not now: run the whole chain from the *original* models
(`Models_original/*.xml`, with real uptake and maintenance bounds) rather than from the
pre-split ones. That would exercise `normalize_bounds` (whose `lb > 0 -> 0` rule deletes
maintenance fluxes and turns a polytope back into a cone) and `split_into_elementary`, and would
test whether the decomposition itself is a source of disagreement. It is a larger experiment: it
means reproducing their splitting step, and the complex indices would no longer line up with the
published results.

**An objective bound.** `relative_tolerance_bound(0.999)` leaves a near-optimal face that is still
hugely degenerate at genome scale, so it changes the feasible set enough to break comparability
while removing almost none of the ambiguity -- and it is undefined on the pre-split models, which
carry no objective. Cheap to add later if a reviewer asks, since it is one flag on the same
machinery.

## Hardware, per arm as designed

| arm | solver | CPUs | note |
|---|---|---|---|
| Kimonu | SCIP | 1 | deterministic; verified bit-identical over 3 replicates on all 33 models |
| COCOA | HiGHS | 64 | batches candidate pairs over Distributed workers; seed-dependent |
| reference | GLPK (Octave) + R | 1 | serial by construction -- no parfor/parpool anywhere in its sources |

Core counts are recorded per cell rather than equalised: parallelism is part of COCOA's design,
and the reference has no parallel path at all.

## Reference caveats

- Its ACR/ACRR come from the **full stoichiometric coupling** (F2C2), a different route than the
  span-membership test Kimonu and COCOA share. A reference-vs-others gap is a method difference,
  not an implementation one.
- F2C2 does not run under Octave (needs `null(A,'r')`, shimmed, and then fails on a NaN from its
  own GLPK path), so the coupling matrices published with the reference are used. That omits
  F2C2's runtime from the reference's total, which **understates** it.
- Concordance is run exhaustively over all candidate pairs, not through the reference's stochastic
  prefilter, so it is deterministic. A run that hits the wall budget is marked `*`: its concordance
  is a subset, so modules come out finer and pairs undercounted, while balanced is unaffected.
