# Benchmark findings

Working notes on what the three-way comparison actually shows. Numbers here are from the
v1.1.0 COCOA arm and the current Kimonu sweep; v1.0.0 results are stale and are not used.

## The "giant module" is a poor headline statistic

`giant_module` is the size of the largest block of the kinetic partition. It is an order
statistic, not a property of the network, and it is the least stable number in the table.

On iAF1260b, Kimonu reports a giant of 96 and COCOA reports 40. Read as a ratio that looks like
a 2.4x disagreement. The partitions say otherwise:

| | blocks | non-singleton | complexes covered | largest |
|---|---|---|---|---|
| Kimonu | 13135 | 1091 | 2632 | 91 |
| COCOA (`kinetic_efficient=true`) | 13141 | 1096 | 2631 | 40 |

Six more blocks, five more non-singletons, one fewer complex covered. Exactly one Kimonu block is
split by COCOA, into 8 pieces; the other 13134 are identical. That one block happens to be the
largest, so the max-statistic swings and every other column stays flat. (91 rather than 96 because
name matching drops ~300 complexes, see below.)

The statistic is not even directionally stable. COCOA's giant exceeds Kimonu's on
ArabidopsisCoreModel (28 vs 18) and on iJN678 at `kinetic_efficient=false` (1132 vs 537), and falls
far below it on iAF1260b (40 vs 91).

Report it for comparability with prior work, which publishes it, but carry the claims on the
number of modules, the number of non-singleton modules, and the complexes covered by them. Those
are monotone under fragmentation.

## iAF1260b: the 96-complex module is verified, independently of Kimonu

Checked against the definition using only the CRN and an LP solver -- no exact kernel, no grouping
code. The module is 50 concordant complexes plus 46 balanced ones.

- The 46: LP max and min of `A_c.v` are both exactly 0. Genuinely balanced.
- The 50: all active, activity range up to 51.9. Over 200 flux points, including 100 chosen to
  maximise and minimise their activity, they form a single proportional group. Ratios come out at
  exactly 0.5 and 2, spread 1.7e-16.

Control on the same sampler: none of Kimonu's 4004 balanced complexes show activity across 60
random vertices, while 2669 non-balanced ones do.

COCOA splits that verified class into 7 pieces (14/9/8/6/5/2/1), identically at
`kinetic_efficient` true and false. Both runs share seed 1234, so this is not yet evidence of
seed-stability -- seeds 2026 and 31337 are queued for both modes and will settle it. Either
outcome is informative: stable means a structural blind spot, variable means sampling variance
that needs error bars.

Kimonu reproduces COCOA's own 40-complex module exactly (its module 109, 20 concordant +
20 balanced, all 40 present). The relation on this model is containment, not conflict: COCOA's 40
is "giant" only because the 96 was missed.

Why sampling misses it: the module's complexes are invisible to generic flux sampling. Random
vertices never activate them -- targeted LPs were needed to make them nonzero at all. A pairwise
ratio test has no signal on those rows and defaults to splitting.

## Neither partition refines the other

The clean claim "COCOA's partition refines Kimonu's" is false. Across 27 model-mode pairs there
are 178 COCOA blocks that straddle two Kimonu blocks, i.e. COCOA merges complexes Kimonu
separates. The two error directions are asymmetric in character:

- Splits are rare but large: 0-3 blocks split per model, but one block can fragment into 8
  (iAF1260b), 22 (iSB619, iYS854) or 31 (iJN678 at `true`).
- Merges are more numerous but small, except iJN678 at `false`, where COCOA's largest module is
  1132 against Kimonu's 537 and the block count drops from 5042 to 4718.

`kinetic_efficient=true` has consistently fewer merge violations than `false` on the same model
(iLJ478 16 -> 0, iNF517 10 -> 0, iSB619 5 -> 0, iIT341 9 -> 2); iJN678 is the exception (1 -> 3).

## Two gotchas that cost time

`sweep/iJB785_fixed.crn.json` and `sweep/iSynCJ816_fixed.crn.json` still contain blocked
reactions (363 and 280). The Kimonu runner strips them with its own FVA, and COCOA is handed the
already-stripped copies under `unblocked/`, so the benchmark itself compares like with like -- but
any analysis that loads from `sweep/` without the FVA is working on the wrong network and inflates
Kimonu's modules (iJB785 to 1091, iSynCJ816 to 1267). The published reference giant for iJB785,
1089, is consistent with that same omission.

COCOA's `complex_ids` drop stoichiometric coefficients, so distinct complexes can share a name
(~300 of 14973 on iAF1260b). Cross-tool comparisons by name must either skip the collisions or
accept them only when every candidate on both sides sits in one block.

## The reference: exact on balanced, cascading on concordance

Nine models now have a complete exhaustive reference run -- the four iIS312 variants, plus
Saha2011iRS1597, iCN900, iIT341, iJB785 and iLJ478 from `reference_full/`, each 82-138 h.

**Balanced complexes agree exactly on all nine.** 379, 315, 379, 228, 817, 744, 696, 4285, 778 --
Kimonu reproduces every one, in 0.3-10 s against 3.7-138 h. That stage of the reference uses no
pair LP, and it is solid.

Concordance diverges on three models, and the divergence has a single shape:

| model | ref classes | K classes | ref largest | K largest | LP non-optimal |
|---|---|---|---|---|---|
| iCN900 | 353 | 353 | 71 | 71 | 0.2% |
| iIS312 | 199 | 199 | 99 | 99 | 3.1% |
| iJB785 | 217 | 228 | 149 | 90 | 21.7% |
| iLJ478 | 91 | 268 | 2194 | 777 | 8.9% |
| Saha2011iRS1597 | 169 | 362 | 1916 | 603 | 16.1% |
| iIT341 | 29 | 174 | 2315 | 913 | 9.8% |

Note that the LP failure rate does not predict the divergence: iJB785 has the worst rate and
nearly agrees, iCN900 and iIS312 agree exactly, iLJ478 diverges most at 8.9%.

### Mechanism

`run_reference.m` builds the concordance relation pairwise by LP, then forms classes as connected
components of that graph (`E = (CC == 2) | (CC == -1)`, symmetrised, BFS). Concordance is an
equivalence relation, so imposing transitivity this way is sound only if the pairwise test has no
false positives. It does have them, and a single spurious edge fuses two entire classes
irreversibly.

That predicts the observed shape, and the shape is what we see:

- The damage is concentrated, not diffuse. On iLJ478 only one reference class merges Kimonu
  classes -- but it absorbs 184 of them, Kimonu's true 777-complex class plus 183 others.
- Blowup happens only where a large true class already exists. A class of n complexes carries
  O(n^2) pairs, so it is the likeliest place for a spurious edge to attach, and once attached the
  component grows. Largest true class 603/913/777 becomes 1916/2315/2194; largest true class
  71/99/65 stays 71/99/65.

Checked against the definition rather than against Kimonu: of the 2194 complexes in the iLJ478
blob, 169 carry nonzero activity over 90 random flux vertices, and those 169 span **19 distinct
proportional directions**. A concordance class must span one. The merge is provably wrong, from
the definition alone.

Kimonu is structurally immune. It groups identical normalised rows of M = A.K, computing the
equivalence classes directly, so there is no pairwise verdict that could be wrong and nothing for
transitivity to propagate.

### Root cause: the scaling variable sits in the objective

The rounding below is real but secondary. The primary defect is the objective vector of the pair
LP. Charnes-Cooper turns max (A_i.v + alpha)/(A_j.v) with A_j.w = 1 into max A_i.w + alpha*t; the
numerator has no constant term, so alpha = 0 and the scaling variable must not appear. The code
uses `[-A_i 1]`, so the solver trades the ratio off against t and stops somewhere that is not the
ratio's extremum. The reported range collapses towards zero, always in the same direction:

```
pair (1842,2368)   [-A_i 1]: max 0          min 0            -> "constant", accepted
                   [-A_i 0]: max 0          min -7.299e-06   -> varies
pair (2244,1625)   [-A_i 1]: max 4.912e-07  min 0
                   [-A_i 0]: max 6.090e-01  min 0
pair (2403,3295)   [-A_i 1]: max 0          min -1.287e-04
                   [-A_i 0]: max 0          min -1.840e+02
```

Six orders of magnitude on the last two. This is why the false positives cluster near zero, and
why `round(.,2)` then waves them through: the objective squashes the range into the bin before the
tolerance ever looks at it. Correcting the tolerance alone does not fix the pair (1842,2368),
whose collapsed bounds are exactly equal; correcting the objective does.

### The pairwise test, and why it fails in both directions

`find_concordant_complexes.m` decides a pair by maximising and minimising the activity ratio
under a Charnes-Cooper transformation and comparing the two results **rounded to two decimals**:

```matlab
CC(i,j) = (round(Maximum_c_p,2) == round(Minimum_c_p,2))*2;
```

That single line is an absolute tolerance of 5e-3 on the ratio, and it is the false-positive
source. It is maximally permissive near zero, where a pair whose ratio is not constant at all
still lands in the same bin. Reproducing the reference's own LP on pairs that bridge different
Kimonu classes inside the iLJ478 blob:

```
    2244,1625   max 0.00000049   min  0.00000000   round2 equal -> concordant
    2403,3295   max 0.00000000   min -0.00012865   round2 equal -> concordant
```

That is not a rare event. Of the 482,007 pairs the LP accepts on iLJ478, **188,877 bridge two
different Kimonu classes** -- 39% of every edge it draws. The same count is 285,606 of 659,005 on
iIT341 and 219,009 of 404,760 on Saha2011iRS1597, against 297 of 15,130 on iJB785 and **zero** on
iCN900, the two models that agree.

Sampling 24 of those accepted bridging edges on iLJ478 and asking the definition: 16 have a
demonstrably non-constant ratio, 3 are genuinely concordant, 5 are inactive on the sampled fluxes
and stay undecided. Every failing ratio sits near zero with a spread inside one rounding bin --
2.1e-6, 2.6e-7, 7.3e-6, 1.1e-4, 3.4e-3 -- which is exactly where `round(.,2)` cannot tell a
varying ratio from a constant one.

So the blob is not a cascade from a few bad edges. The pairwise test accepts wrong pairs by the
hundred thousand, and connected components merely assembles them.

The false negatives come from a different line. When a complex is in neither the all-positive nor
the all-negative group, the verdict ANDs the two sign branches:

```matlab
CC(i,j) = (round(Maximum_c_p,2) == round(Minimum_c_p,2) & ...
           round(Maximum_c_n,2) == round(Minimum_c_n,2) & ...
           round(Maximum_c_p,2) == round(Minimum_c_n,2))*2;
```

An infeasible branch leaves `Maximum = Inf, Minimum = -Inf`, so the first clause fails and the
pair is rejected however well the other branch behaved. Reproducing the LP on Kimonu class 204 of
iLJ478, the twelve complexes the reference cuts into six pairs:

```
  kept  (58,58)   t>0: ratio -1.0 constant     t<0: infeasible
  split (58,59)   t>0: infeasible              t<0: ratio 1.0 constant
  split (58,60)   t>0: infeasible              t<0: ratio 1.0 constant
  split (59,60)   t>0: infeasible              t<0: ratio 1.0 constant
```

The pairs it keeps are normalisable with t > 0; the ones it drops are normalisable only with
t < 0. Nothing about them is numerically marginal -- the surviving branch returns exactly 1.0.

### Kimonu is right in both directions

Checked against the definition with flux vectors, not against either tool:

- Every Kimonu class the reference splits is genuinely **one** proportional direction. Over 48-64
  flux points including targeted extremes, the worst ratio spread across all six cases is 6.3e-11,
  and four of the six are exactly 0.
- The reference's iLJ478 blob is genuinely **19** directions among the 169 of its 2194 members
  that carry activity. A concordance class must be one.

So the reference over-merges by rounding and under-merges by branch ANDing, and Kimonu matches the
definition on both counts. Its grouping of identical normalised rows of M = A.K has no pairwise
verdict to get wrong and no tolerance to tune.

## The reference's kinetic stage: one line of NA indexing

The R stage was never run during the earlier work, so its numbers were taken from the published
`.RData` files. Running it is cheap -- about five minutes a model -- and `run_reference.m` already
writes the `.mat` it consumes.

On iCN900 all three implementations agree **exactly** on concordance: 353 classes, largest 71. The
kinetic stage then diverges by a factor of 33.

| | modules | giant |
|---|---|---|
| Kimonu | 250 | 19 |
| COCOA | 250 | 19 |
| reference | 1 | 626 |

The inputs are not the problem: the R graph and its strong components match Kimonu exactly -- 1601
SCCs, largest 8, 1287 terminal, 744 balanced, 353 classes. The cause is a single line in
`autonom`, phase IV:

```r
terminal <- rep(1, length(clust$membership))   % one entry per complex in the whole graph
...
phaseIV <- x[terminal == 0]                    %%which complexes are non-terminal
```

`x` holds complex ids; `terminal` is a mask over every complex in the graph. Indexing `x` by a
logical of a different length selects members of `x` by *position* and pads with `NA` past its
end, so phase IV returns a positional sample of the set plus one NA. The comment states the
intent, which is `terminal[x]`.

Correcting that one subscript, with nothing else touched:

```
as published       NAs 352 | sets 352 total 958 union 625 overlap 333 | blocks 185 largest  93
terminal[x] == 0   NAs   0 | sets 253 total 586 union 586 overlap   0 | blocks 253 largest  19
Kimonu             NAs   - | sets 257 total 594 union 594 overlap   0 | blocks 250 largest  19
```

The NAs disappear, the upstream sets become pairwise disjoint as Kimonu's already are, and the
largest kinetic module matches Kimonu and COCOA exactly. Module counts remain 3 apart, 253 against
250, which is not yet explained.

Worth being clear about why this is not circular: the correction follows from the code's own
comment and from R's indexing semantics, not from comparing against Kimonu. Agreement with Kimonu
is the result, not the criterion, so a corrected R stage still counts as an independent third
measurement.

The NAs also propagate. In the published script they are carried into the pairwise-intersection
matrix and the connected-components merge that follows, which is how 352 sets of median size 3
end up as a single module of 626.

### Kinetic stage, second finding: phase III's output bypasses the terminal filter

With the NA subscript corrected, the reference's kinetic partition still differs from Kimonu's on
iCN900: 246 non-singleton modules against 250, covering 583 complexes against 591, with 96 of the
modules differing by a single complex at their boundary.

Two things had to be matched first. `class_with_balanced` stores only the non-singleton
concordance classes -- `module_of` is -1 for singletons -- so the R stage seeds 352 upstream sets
where Kimonu seeds 372. Adding the 20 singleton classes as seeds brings the aggregate figures into
line: 257 sets, union 594, zero overlap, on both sides.

What remains is directional. Of the complexes the reference keeps and Kimonu does not, **105 of
105 are globally terminal**; of those Kimonu keeps and the reference does not, 7 of 113 are. The
return value explains it:

```r
phaseIII <- setdiff(x1, x)        # what phase III removed, returned unfiltered
phaseIV  <- x[terminal[x] == 0]   # survivors, terminal filter applied here only
return(union(phaseIII, phaseIV))
```

Terminal complexes removed by phase III re-enter through `phaseIII` and never meet the filter that
phase IV exists to apply. Kimonu and COCOA both report 250 modules with a largest of 19 on this
model; the reference is the outlier.

Hypotheses tested and refuted on the way, recorded so they are not retried: Kimonu's coupling and
deficiency merges (identical with Y=None), the iteration count (identical at max_iterations=1),
index misalignment (zero composition mismatches), phase II's out-degree-zero removal (none of the
differing complexes have out-degree zero), the phase structure itself (a two-phase reimplementation
in R gives results identical to the four-phase one), and computing strong components on the induced
subgraph rather than the full graph (far worse: 305 sets spanning 1480 complexes).

Neither `union(phaseIII, phaseIV)` nor `phaseIV` alone reproduces Kimonu -- the latter is too
aggressive, giving 238 sets with a largest of 12 -- so the correct single-line form is not yet
settled and no fix is proposed here.

### Resolved against the specification: the R stage runs the wrong algorithm

Langary et al. (2025), Supplementary Material S.2.2-S.2.3, gives two algorithms. Section S.2.2
presents a *tentative* four-phase version finalising as `C_up = C_-4 u C_-3`; Remark S2-1 then
gives the one the paper names the **Upstream Algorithm**:

```
1. Phase I : exclusion of entry complexes, to a fixed point
2. Phase IV: exclusion of terminal complexes
             identify all terminal strong linkage classes in G_{C_-4}, remove them all
3. Finalization: C_up = C_-4
```

Kimonu implements this. `code_kineticModule_analysis.R` implements the tentative four-phase
version instead, so its `union(phaseIII, phaseIV)` is faithful to the algorithm it implements --
that line is not a bug, and an earlier note here suspecting it was wrong.

Running the tentative algorithm inside Kimonu's own framework reproduces the reference exactly --
246 non-singleton modules, largest 19, covering 583 complexes, against the R stage's 246 / 19 / 583.
Both implementations are therefore correct for the algorithm each one runs; the only thing
separating them is which of the paper's two algorithms they chose.

With that established, the outputs differ by exactly 8 complexes, and the tentative result is a
strict subset of the Upstream Algorithm's:

```
tentative (S.2.2)      246 non-singleton   covered 583   globally terminal inside 0
Upstream Alg (S2-1)    250 non-singleton   covered 591   globally terminal inside 8
```

Lemma S2-4 says any complex in an upstream set is nonterminal, so those 8 should not be there. They
survive because Kimonu's phase IV tests `succ[v] - c` against the whole network: a strong component
whose members have successors outside the candidate set is judged nonterminal and kept, where the
specification asks for terminal strong linkage classes of `G_{C_-4}`, the reduction of the network
to the candidate set. A component with no successor inside the set but successors outside it is
terminal in `G_{C_-4}` and nonterminal under Kimonu's test.

This is small -- 8 complexes of 3329 on iCN900, and the giant is 19 either way -- but it is a
deviation in Kimonu, not in the reference.
