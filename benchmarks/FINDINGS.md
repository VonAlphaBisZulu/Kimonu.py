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

### A naming mismatch invalidated several comparisons here

COCOA writes a stoichiometric coefficient as a prefix on the species, `2_glu__L[c]`, where the
Kimonu export carries it as a separate field. Comparisons that matched complexes by name with the
coefficient dropped therefore failed on every multi-coefficient complex -- 54 of 3328 on iCN900,
and precisely the boundary complexes the partitions disagree about. Several per-complex figures
reported in earlier revisions of this file were produced that way and were wrong.

Rendering both sides with the same convention leaves 2 complexes of 3328 unmatched and no name
collisions. On that bijection:

| | non-singleton | largest | covered | globally terminal inside |
|---|---|---|---|---|
| tentative (S.2.2) | 246 | 19 | 583 | 0 |
| Kimonu (S2-1) | 250 | 19 | 591 | 8 |
| COCOA | 249 | 19 | 589 | 113 |

The giant is 19 in all three. The terminal content is not explained: all three implement a phase IV
that removes terminal complexes, and COCOA's does so against the full network unconditionally
(`setdiff!(current_indices, network.global_terminal)`), yet its reported modules retain 113. The
`kinetic_module` labels in the output files are a partition of every complex rather than the raw
upstream sets, so the difference may sit in what happens after phase IV rather than in phase IV.
Settling it needs COCOA instrumented at that point, not more inference from its output.

`compare_kinetic_partitions.py` holds the corrected bijection and `tentative_upstream.py` the
S.2.2 algorithm implemented against Kimonu's topology, which reproduces the reference exactly
(246 / 19 / 583).

### The CPLEX bridge was order-dependent

`lp_server.py` accepted only CPLEX statuses 1/101/102. Status 5, OPTIMAL_INFEAS -- an optimal
basis carrying small infeasibilities after unscaling, which the pair LP's 1e9 bounds provoke --
was returned as a failed solve, and whether a given LP landed on 1 or 5 depended on the warm basis
the slot carried in. So a pair's verdict depended on which pairs came before it.

Measured on iCN900 with all three concordance fixes in: the first pair of a cold process,
(233,243), came back non-concordant, then concordant on the immediately following identical run.
A cold, independent solve gives t<0 max = min = 1.0 exactly and Kimonu places both complexes in one
10-complex module, so the cold-start verdict was the wrong one. Over 300 random accepted pairs on a
warm process there were no flips, so the exposure is cold starts and slot evictions rather than
steady state; but on Saha2011iRS1597 and iJB785 14.8% of pairs hit status 5.

This also explains why the sampled mode had accepted edges the exhaustive mode rejected (312 on
Saha2011iRS1597, 4864 on iIT341): those are the same LP solved from different slot states. The
prefilter itself only removes pairs from `At` and never sets `CC`, as documented.

Fix: on any status outside 1/101/102 the worker re-solves from scratch (`advance = 0`), which does
not depend on order, and the shim accepts 5 after that and counts it separately (`linprog('degraded')`).
After the patch: (233,243) = 2 in both runs, 15/15 and 300/300 with zero disagreements.
`test_determinism.m` is the check. Every reference concordance run before this patch is suspect
at the level of individual marginal pairs; the two full-model reruns in `reference_fixed2/` are
the gate before any of the nine is redone.

### Kimonu's two deviations from Remark S2-1, measured across the 33

The outer iteration loop in `kinetic_modules` changes the partition on **0 of 33** models. It can
be documented as inert on this benchmark or removed.

Judging phase IV terminality on the reduced graph `G_{C-4}` rather than the full network changes
the partition on **32 of 33** models and the giant on most of them (iJB785 22 to 20, iAF1260b 96
to 93, the iAM_P* family 100 to 98). The choice is not cosmetic. It stays as it is: the definition
of an upstream set requires every complex to have a path leading to 𝒞 \ C_up, and a component that
still feeds complexes outside the candidate set satisfies that regardless of what the reduced graph
shows. COCOA reads the remark the same way and says so in its source.

### The tentative and Upstream algorithms across the 33

Running the §S.2.2 four-phase algorithm inside Kimonu on every model: the giant agrees with the
two-phase Upstream Algorithm on 25 of 33 and differs by 1-4 complexes on the other 8 (the iAM_P*
family, iJR904, iSB619, iYO844). Module counts differ far more on iJB785 (135 vs 1905) and
iSynCJ816 (296 vs 1597) with the giant unchanged; that gap is not yet explained.

### The 8 terminal complexes in Kimonu's modules are by design; COCOA's are not the same thing

Traced on iCN900: all 8 are balanced, sit in two-complex terminal strong linkage classes, and
appear in **none** of the raw `upstream_set` outputs. They enter at `kinetic.py:303`, which adds
every weak linkage class composed entirely of balanced complexes as its own module -- the step the
reference's R code performs at its lines 141-155 and COCOA in `add_singleton_balanced`. So the
earlier note calling this a Kimonu deviation is withdrawn: Phase I reaches a fixed point, no
terminal class straddles the candidate set, and `upstream_set` itself never returns a terminal
complex. Full-graph terminality stays.

COCOA's retained terminal complexes are a different population: of the 113 on iCN900, 664 on
Chlorella_variabilis_iAJ526, 550 on iAB_RBC_283 and 1027 on iAF1260b, only 0-6% are balanced. Its
`compute_global_terminal_complexes` and `add_singleton_balanced` both read correctly, so they enter
through the merge or the per-complex label mapping. Still needs COCOA instrumented.

### The reference ACR stage needs metabolite names in the saved model

`Get_excel_tables_kinetic_modules.R` reads `MODEL_r{1}.metNames`; `run_reference.m` saved only
`A Y S lb ub b c complexes`. Names are now written (`mets`, `metNames`, `rxns`) and the six
existing files were patched with `patch_mat_names.m`.

### Reference ACR stage: first numbers, and how they relate to Kimonu's

With metabolite names in the saved model the reference's ACR script runs. On the five models with
sound concordance it reports ACR = 0 everywhere and ACRR = 1, 1, 1, 1, 3 (iCN900, iIS312,
Epimastigote, Amastigote, Trypomastigote). Two output quirks matter for reading it: `OverviewTable[15]`
is set to 1 when `unique(diff_two_met[,1:2])` collapses to a vector, which happens for exactly one
pair, not zero -- the empty case never enters that branch -- and `write.table` writes a single pair
as a two-row column rather than one row.

The definitions differ. The reference calls a species ACR when two complexes in one kinetic module
differ in exactly that species, and a pair ACRR when two complexes differ in exactly those two,
one each -- the combinatorial sufficient condition of Corollary S2-10. Kimonu tests
e_i in im(Y Delta) on the kernel, which is the necessary-and-sufficient criterion for a coupled
module. So the reference's sets should be contained in Kimonu's, and they are, by name, on all
five models (table in `acr_contain.txt` at the time of writing): every reference ACRR pair is in
Kimonu's `acrr`, and the reference's ACR set is empty where Kimonu's is empty. Kimonu additionally
finds 8 / 5 / 4 / 4 / 1 free ACRR pairs the combinatorial test cannot see.

Kimonu's `acr`/`acrr` count enzyme pseudo-species; the reference's tables also do, so those are
the comparable columns. The `_free` counts are what both prior papers report.

### The module-count gap between the two algorithms is the balanced linkage classes

The 4-phase run inside Kimonu gave 135 modules on iJB785 against Kimonu's 1905, and 296 against
1597 on iSynCJ816, with identical giants. The difference is exactly the number of weak linkage
classes composed entirely of balanced complexes -- 1770 and 1301 -- which Kimonu, the reference's
R code and COCOA all add as modules and the tentative reimplementation did not. Not a discrepancy.

`coupling_exact/` matrices are now also copied after the writer finishes; the first ACR chain
raced it on iJB785.

### Retraction: COCOA does not retain terminal complexes in its kinetic modules

Earlier sections here report COCOA keeping 113 / 550 / 664 / 1027 globally terminal, non-balanced
complexes in its kinetic modules on iCN900, iAB_RBC_283, Chlorella_variabilis_iAJ526 and iAF1260b,
and treat that as an open defect. Those figures were produced by intersecting COCOA's module
membership, mapped into the Kimonu-JSON index space (zero complex dropped), with a terminality set
computed in CRN index space (zero complex present). Above the zero complex the two spaces are off
by one. The Kimonu and tentative-algorithm figures in the same tables were computed in CRN space
throughout and stand.

Redone in one space, COCOA's modules retain 8 / 0 / 10 / 4 globally terminal complexes -- exactly
Kimonu's 8 / 0 / 10 / 4, and all of them balanced linkage-class modules. COCOA's balanced set is
identical to Kimonu's exact one on all four models. And fed Kimonu's exact concordance partition
directly, COCOA's kinetic stage returns 249 non-singleton modules, largest 19, covering 589, with
8 terminal inside, at both `efficient` settings, against Kimonu's 250 / 19 / 591 / 8; the two
complexes' difference is within the two names the bijection could not match. There is no COCOA
kinetic defect to instrument.

This is the fourth cross-tool join error in this file. The three before it were name parsing.
The safeguard now: every cross-tool comparison converts to CRN index space before any set
operation, and the conversion is asserted against complex compositions.

### Gate and ACR containment, continued

iIS312_Trypomastigote on the fixed algorithm and the fixed bridge: balanced 228/228, classes
111/111 identical, and `lp_nonoptimal 0, lp_degraded 0` where the old bridge logged 4363 / 2640.
Amastigote is running.

iJB785's reference ACR stage: giant 33 against Kimonu's 22, inherited from its 297 false-positive
concordance edges; ACR 0 on both; 17 ACRR pairs, all 17 in Kimonu's `acrr` (14 in `acrr_free`).
Containment now holds on all six models with sound-or-adjudicated concordance.

### On identical input the two kinetic stages agree exactly

With the complex bijection completed (fractional coefficients rendered as COCOA writes them), COCOA's
kinetic stage fed Kimonu's exact concordance partition for iCN900 returns 250 non-singleton
modules, largest 19, covering 591 complexes, 8 globally terminal inside -- at both `efficient`
settings -- which is Kimonu's result to the complex. Every earlier Kimonu-vs-COCOA kinetic
difference therefore lives upstream in the concordance input, not in either kinetic stage.

`reference_sampled_seeds/` holds seeds 2026 and 31337 for eight of the nine models; iLJ478 has
only seed 1234. The staged final sampled rerun covers all 27 cells regardless.

### What the preliminary figures make visible about COCOA's seeds

COCOA `efficient=false` has occasional catastrophic seeds. On iIT341 one seed of three returns a
giant of 523 with 31 free ACR metabolites and 2089 free ACRR pairs; the other two return 185-190,
0 and 20-24, against Kimonu's 195, 0 and 24. On iHN637 likewise: one seed 825 / 25 / 2211, the
other two 211-213 / 0 / 116, against Kimonu's 243 / 0 / 121. In both cases a single sampled
concordance partition over-merged and everything downstream inflated with it. On iJN678 the
over-merge is stable across all three seeds (giant 1183, ACRR 2484, against 551 and 236), so it is
a property of that model under sampling rather than seed luck.

For the paper this argues for reporting the median over seeds with the individual seeds shown as
points, rather than mean with min-max bars: a mean of (523, 185, 190) is not a number anyone
computed. The preliminary figures still use mean and min-max.

### Gate passed; final reference reruns launched (2026-09-18)

iIS312_Amastigote on the fixed algorithm and the fixed bridge: balanced 315/315, classes 165/165
identical, zero differences, `lp_degraded 0` (`lp_nonoptimal 27` are the deliberately infeasible
negative branches, not failures). With Trypomastigote this makes two complete models reproducing
Kimonu's concordance exactly on the final pipeline.

Launched: nine exhaustive runs, jobs 549287-549295 (`ref_final.sbatch`, 14-day cap, output
`reference_final/`), and the 27-cell sampled array 549296 (`ref_final_sampled.sbatch`, 7-day cap,
output `reference_final_sampled/`). These replace every earlier reference concordance run.

The first submission (549251-549260) failed at once with no log: both scripts write their Slurm
log into the output directory, which the script itself creates only after Slurm has tried to open
the log. Directories created, resubmitted as 549287-549295 and array 549296.

Old-bridge exhaustive runs 548721 (iNF517) and 548722 (iSB619) cancelled at ~1.8 days; both models
added to the fixed rerun instead: exhaustive 549302 549303 (21-day cap) and six sampled cells as array 549304
(cells 27-32). The final reference set is therefore 11 models.

The reference jobs requested 16 GB each; a completed gate run's MaxRSS was 151 MB. With fifteen
jobs that booked 240 of node01's 257 GB while using 30 of 64 CPUs, so the added runs sat pending.
Pending jobs and array elements were reduced to 4 GB in place and both scripts now request 4 GB;
the nine already running keep their 16 GB until they finish.

### Thirteen sampled cells out of memory at 4 GB; resubmitted per model (2026-09-18)

The 4 GB right-sizing was taken from an *exhaustive* gate run (151 MB). Sampled mode is different:
`pre_concordant_sampling.m` line 84-85 forms `(sample+eps)./(sample_reshaped+eps)`, an m x 10 x m
dense double array (m = complexes), twice per round (`std` and `mean`), so the prefilter peaks at
about 7 x m^2 x 80 B. Measured: iIS312_Epimastigote 1738 complexes 1.5 GB, iCN900 3329 6.0 GB,
iJB785 6114 20.4 GB (old whole-node run 546238_15). Cells 17-23 and 27-32 (iIT341 seed 31337,
iJB785, iLJ478, iNF517, iSB619) were killed in rounds 1-4 of the prefilter; nothing was written.

Cells 15, 16 (iIT341) and 24-26 (Saha2011iRS1597) survived at 4 GB only because node01 has 16 GB
of swap: their prefilter ran at the cgroup cap (sstat MaxRSS 4.18 GB) and took 1160-2770 s. They
are now in the concordance phase at 200-430 MB RSS and their results are unaffected (swapping
changes time, not arithmetic); only their prefilter wall time is inflated, about 1-2% of the cell's
total runtime on the old-run scale (iIT341 sampled took 21 h whole-node).

Resubmitted with per-model memory (8 x m^2 x 80 B + 0.5 GB, rounded up), same script, same cells,
same seeds, `--mem` on the sbatch command line: 549880 (cell 17, 8 GB), 549881 (18-20 iJB785,
26 GB), 549882 (21-23 iLJ478, 9 GB), 549883 (27-29 iNF517, 12 GB), 549884 (30-32 iSB619, 12 GB).
The script default is now 8 GB with the scaling rule in a comment. iJB785 books 26 GB for a run that
needs it for the first ~45 min only; single-job simplicity was preferred over splitting the
prefilter into its own job. Node01 is memory-bound by the running 16 GB exhaustive jobs (booked
204 GB / 257 GB), so the 26 GB and 12 GB cells start as the iIS312-family exhaustive runs finish.
Running jobs cannot be shrunk (`scontrol`: "Job is no longer pending execution").

### Sampled reference arm complete; exhaustive landed on the iIS312 family (2026-09-21)

All 33 sampled cells finished (0.2-53 h each; the 13 resubmitted cells ran at their per-model
memory without incident, iJB785 peaking at 20.5 GB under its 26 GB request). Compared to Kimonu by
complex composition (`compare.py` partition): balanced sets identical on every cell; concordant
pairs a strict subset of Kimonu's on every cell, zero invented pairs, `lp_degraded 0` throughout.
Seeds agree with each other except iCN900 31337 (one pair fewer) and iSB619 31337 (2158 fewer).
Module counts are not monotone in the pairs: iIT341 178 vs Kimonu 173 and iLJ478 300 vs 267,
because a dropped pair splits a module. This is the shape the paper should show: the sampled
prefilter loses concordances (iSB619 -9504 to -11662 pairs of 177160), never gains them.

Exhaustive: iIS312, iIS312_Amastigote, iIS312_Epimastigote, iIS312_Trypomastigote landed (3-17 h),
each identical to Kimonu in balanced set, modules and pairs. The reference R stages (fork's fixed
kinetic code, exact coupling matrices) ran on their .mat files in `reference_final/acr/`: giants
9/9/9/11 equal Kimonu's; ACR 0 on all four; ACRR 1/1/1/3 pairs, every one in Kimonu's set by name
(Kimonu 9/9/9/11: kernel test vs the reference's sufficient combinatorial criterion). Seven
exhaustive runs continue. Exact coupling matrices now exist for all 11 final models.

`figures/gather_figdata.py` rebuilds `figdata.json` from the result directories; with `--check` it
reproduces the committed Kimonu and COCOA numbers exactly (only new COCOA efficient=true seeds
differ). Reference arms now come only from `reference_final*/`; figures F1-F5 re-rendered.

COCOA efficient=true is at 31/33 (30 at three seeds). It has the same seed-outlier behaviour as
efficient=false: iHN637 seed 2026 gives balanced 1150 (others 1149), giant 808 (213-216), ACRR
88019 (5616-5722), ACR 7 (0). iAF1260b efficient=true giant 40/52/40 over seeds, ACRR 766/953/789.
COCOA efficient=false iAF1260b cell 0 hit its 7-day cap: a non-finisher, not reported, not redone.

### Four-phase and two-phase upstream sets are identical on all 33; the giant difference is the merge rule (2026-09-21)

Measured seed by seed (`kinetic_variants/measure_phases.py`): for every concordance module ∪
balanced seed on every one of the 33 models, the four-phase algorithm of S.2.2 and the two-phase
algorithm of Remark S2-1 return the same upstream set. Zero complexes on either side, 33/33. The
final modules of both satisfy autonomy, feeding and Lemma S2-4 on every model once balanced-only
linkage classes are set aside. Remark S2-2's equivalence claim therefore holds on this benchmark.

**Retraction.** Earlier entries here and in the ledger called Remark S2-2 false ("on iCN900 they are
not identical", "giants agree on 25 of 33"). The comparison behind that mixed three bookkeeping
differences into the algorithm: singleton upstream sets counted on one side only, balanced-only
linkage classes added by Kimonu but not by the tentative script (iCN900 246 → 250 modules, 583 → 591
covered, the "8 terminal complexes"), and above all the merge rule. `kinetic_variants/giant_merge.py`
separates them: on all 33 models the four-phase and two-phase giants agree under overlap-only
merging, and the eight models where Kimonu's giant differs (iAM_Pb448/Pc455/Pf480/Pk459/Pv461 100–101
vs 96–97, iJR904 146 vs 144, iSB619 140 vs 139, iYO844 102 vs 100) differ only because Kimonu also
merges by Proposition S3-4 (modules whose representative composition difference lies in im(YΔ)).
The reference's R kinetic script merges by overlap only (`km_fixed.R`, `mdiff`/`intersect`); COCOA
applies the S3-4 merge with `efficient=false` and skips it with `efficient=true`
(`kinetic_analysis.jl` docstring). The "tentative vs S2-1" framing is withdrawn; the arm previously
labelled "Kimonu (tentative §S.2.2)*" is Kimonu with overlap-only merging, i.e. the reference's merge
rule, and is relabelled as such. `sweep_tentative/` logs remain valid as that arm's numbers.

Consequence for the paper: there is no algorithm choice to defend. Kimonu, COCOA and the reference
compute the same upstream sets; the published theory's S3-4 merge is implemented by Kimonu and
COCOA (`efficient=false`) and not by the reference code. The earlier finding that reduced-graph
terminality in Phase IV would change 32/33 partitions is unaffected (that variant is a third,
unpublished reading and stays unadopted).

### Decision: Proposition S3-4 merge stays Kimonu's default (2026-09-21)

Measured on the affected models: the S3-4 merge and overlap-only merging give identical ACR and
ACRR sets (iAM_Pb448 1422/1422, iJR904 2553/2553, iSB619 1841/1841, iYO844 1140/1140 ACRR; ACR
0/0 or 3/3), so the rule changes only how coupled complexes are grouped, never the reported
metabolites. Modules disconnected in the reaction graph occur under both rules (24-33 per model):
concordance is a flux property, not adjacency. Decision: S3-4 stays the default because it is the
closure of the coupling relation and the published theory; overlap-only giants are reported as the
starred comparability arm (F3) with one Methods sentence. `figdata.json` key renamed
`giant_4phase` -> `giant_overlap`.

### Roadmap additions (2026-09-21)

- F2C2 agreement figure for the Supplement: exact full coupling (proportional rows of K) vs the
  shipped F2C2 matrices, per model. Verified identical on iCN900; F2C2 misses 4 reactions on each
  iIS312 variant. Not yet drawn.
- Idea for StrainDesign: an F2C2-type exact flux-coupling function (full/partial/directional
  coupling from the exact nullspace). Community-facing; exactness is the selling point.

### Fifth reference concordance defect: the scaling variable's fixed bound (2026-09-22)

iCN900 exhaustive landed identical to Kimonu (744 balanced, 352 modules, 14618 pairs, 86 h).
iIT341 exhaustive landed as a strict subset: 429935 of Kimonu's 431235 pairs, 185 modules vs 173,
balanced identical, `lp_degraded 0`. The 1300 missing pairs split exactly two Kimonu modules (13
and 79 complexes) into small pieces; the exact M rows of the rejected pairs are proportional
(ratios -1, 8, -4, ...), and HiGHS bounds of the activity ratio are constant.

Mechanism: the Charnes-Cooper LPs normalise the denominator to A_j.w = 1 with w = t v and
t in [0, 999] (negative branch [-999, 0]). A complex whose activity never exceeds 1/999 cannot
reach the normalisation, both branches are infeasible, and the pair is reported non-concordant.
iIT341 has 198 such non-balanced complexes (max |A_c.v| = 1.18e-4, biomass-linked); all 1300
missing pairs involve one. The four landed iIS312-family models and iCN900 have none, which is why
they came out identical. `kinetic_variants/` is unaffected; this is concordance only.

Fix on the fork (commit after 28350a9): bound t per denominator, t <= 999 / max|A_j.v|, the
activity range computed once per denominator complex with two LPs. Validation through the CPLEX
bridge: the 8 hand-adjudicated iLJ478 pairs still 8/8; of the 1300 rejected iIT341 pairs, 1286 are
now accepted and 14 (all with denominator complex 2939, t ~ 8.5e6) stay rejected in CPLEX although
HiGHS solves them to constant ratios 1, -4, -1. Scaling the +-1e9 flux-variable bounds with t makes
CPLEX worse (382 rejected); free bounds break the bridge. Residual: solver precision at large t.

Consequences. The five exhaustive runs still in flight (iJB785, iLJ478, Saha, iNF517, iSB619) use the
fixed bound and will be subsets on any model with sub-1e-3 complexes; the sampled arm's subset
relation mixes this cap with the prefilter. Options: (a) report reference as subset with this cause;
(b) rerun the affected exhaustive runs with the fix; (c) after each run lands, re-test only the
pairs Kimonu accepts and the reference rejected, with the fix, and report "reference + fix = Kimonu
up to N solver-residual pairs" (acceptances are already verified as a subset of the exact answer).
Downstream: iIT341's reference R stage on the split concordance gives giant 245 vs Kimonu 195 and
ACRR 21 (free 8) vs 4048 (free 24); every reference pair is in Kimonu's set.

### Sixth reference defect, kinetic R stage: the balanced class is taken to be the largest class (2026-09-22)

On the retested iIT341 partition (identical to Kimonu's), the reference R stage first gave giant 245
vs Kimonu 195 and 7 of its 21 ACRR pairs outside Kimonu's set. Cause: `code_kineticModule_analysis.R`
picks the largest entry of `class_with_balanced` as the balanced class. iIT341's largest concordance
module has 913 complexes against 696 balanced (iLJ478: 777 vs 778, one short of the same failure),
so the module was seeded as "balanced" into every upstream set. The shipped .mat files carry `B`;
the fork now identifies the balanced class as the one equal to `B` (fork commit after e141a83). On
iCN900 and the iIS312 family the balanced set is the largest class, which is why they were fine.

### Kinetic-module counts: Kimonu seeds singleton concordance classes, the reference does not (2026-09-22)

Non-singleton kinetic modules on the landed models: Kimonu 49/33/49/26/250/103 against the
reference's 49/33/49/26/246/96 (iIS312 family, iCN900, iIT341). Set by set, every reference module
is a Kimonu module; the extra Kimonu modules (4 on iCN900, 7 on iIT341) are pairs of one balanced
and one non-balanced complex, the upstream sets of singleton concordance classes, which the
reference R script never seeds ("if there is more than one complex assigned"). COCOA's counts agree
with Kimonu's where its concordance is clean (iCN900 250 at both settings), so COCOA seeds them too.
Giants and ACR/ACRR are unaffected on every landed model. Whether a singleton class should seed a
kinetic module is a definitional choice for the paper; Kimonu keeps it for now. Reference counts
are read from the R stage's `final_lres` (dumped as `modules_<stem>.json` by `run_chain.sh`);
`Overview[4]-[5]` is not usable for this because the singleton count is taken before the late
singleton append.

COCOA's ACR/ACRR test, from `kinetic_analysis.jl`: `efficient=true` is a pairwise rule on
composition differences within a module (one differing species → ACR; exactly two with equal and
opposite differences → ACRR pair); `efficient=false` runs that rule and then the column-space test
`y ∈ im(YΔ)` by QR projection at tolerance 1e-8. The reference's rule (Corollary S2-10 as coded)
is the pairwise rule with a stricter pair pattern (each complex has exactly one further species in
common) and no sign test. Kimonu applying the reference's rule to its own modules
(`kinetic_variants/s210_criterion.py`) reproduces the reference's ACR and ACRR name sets exactly
on all six landed models; it is the starred arm in F4 and M2.

### Reference sampled mode gets the R stage; main text uses it (2026-09-22)

The reference's kinetic and ACR/ACRR scripts (fixed R code) now run on the 33 sampled-mode
.mat files (`reference_final_sampled/acr/`). On the six models with both modes: giants identical
(11/9/9/9/195/19), free ACRR pairs identical (1/1/1/1/4/1), non-singleton kinetic modules equal on
four and within 5% on two (iIT341 97 vs 96, iCN900 235 vs 246), where the prefilter's dropped pairs
split concordance modules. Decision: the main-text figures M1 and M2 carry COCOA efficient=true and
the reference in sampled mode (11 models, three seeds, whiskers), the supplement figures F3 and F4
carry both reference modes. iJB785's R kinetic stage is slow (hours) and the remaining sampled
stems land behind it; the figures are re-rendered when they do.

### Kinetic modules from a nullspace in complex-concentration space (2026-09-24, exploratory)

Mass action gives v = diag(k)·Sᵀψ. For rate constants built from a positive steady-state flux v*
(k_r = v*_r / ψ*_s(r)), every positive steady state satisfies N·diag(v*)·Sᵀ·u = 0 with u = ψ/ψ*, so
complexes whose rows of a basis of ker(N·diag(v*)·Sᵀ) are proportional have a fixed ψ ratio: a
coupling test with no graph traversal (`kinetic_variants/psi_nullspace.py`). Random rate constants
are the wrong notion of generic: they usually admit no positive steady state and force whole modules
to zero.

With uptake reactions detached from the shared zero complex (otherwise every inflow-driven module is
coupled to ψ_∅ = 1 and they all merge), on 33/33 models: the ψ partition is identical to the
upstream algorithm's overlap-merged modules on 3 models and strictly coarser on 30, never finer and
never crossing; two independent v* draws agree on 33/33; runtime comparable to Kimonu. Extra members
(`psi_extras.json`): 15 uptake-fed complexes, 50 in terminal strong components, 347 other, e.g.
enzyme-substrate complexes of fully coupled elementary steps (Langary's own remark that fully coupled
reactions make their substrate complexes coupled, which the upstream algorithm does not apply).
Against Kimonu's default (S3-4 merge) the relation is "neither" on the S3-4 models: the ψ test is
linear in ψ and does not see the log-space merges.

Status: the test is sound by construction; an independent numerical check (second steady state by
continuation) has not yet converged on any model and is not counted as validation.

### ψ-nullspace + log-space closure versus Kimonu's default, 33 models (2026-09-24)

`kinetic_variants/psi_closure.py`: ψ-nullspace coupling (uptakes detached), then Proposition S3-4 with
known-ACR augmentation iterated to a fixed point. No graph traversal. Against the current default
(upstream algorithm + S3-4 + S4-1 + S4-5): the new partition contains the default on 29/33 (identical
on iAB_RBC_283, iCN900, iLJ478; strictly coarser on 26); on 4 (iJR904, iNJ661, iSB619, iYO844) neither
contains the other, the default having a merge the closure does not make (iJR904 giant 146 vs 144,
iYO844 102 vs 101) -- the deficiency-one lemma, not yet part of the closure, is the likely source.
ACR and ACRR of the default are subsets of the new ones on 33/33. Gains include Saha2011iRS1597
free ACRR 531 -> 949 and ACR 2 -> 3, iJN678 ACR 0 -> 1, iMM904 free ACRR 22 -> 38, iJB785 26 -> 36.
Runtime 1-87 s, about twice the default (the extra LP for v* and the second exact nullspace).
