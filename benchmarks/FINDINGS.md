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
