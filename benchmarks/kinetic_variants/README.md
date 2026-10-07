Four-phase (S.2.2) vs two-phase (Remark S2-1) upstream sets, and the module merge rule.

`measure_phases.py` runs both algorithms seed by seed (concordance module ∪ balanced) on Kimonu's
topology and checks autonomy, feeding and Lemma S2-4 on the final modules with balanced-only
linkage classes excluded. `giant_merge.py` reports the largest module under three rules: four-phase
with overlap merge, two-phase with overlap merge, two-phase with the Proposition S3-4 coupling
merge (Kimonu's default). Logs are the 2026-09-21 runs over all 33 models.

Result: the two algorithms return identical upstream sets for every seed on all 33 models; every
giant difference between Kimonu and the reference is the S3-4 merge, which the reference's R
kinetic script does not perform.

## Two steady states, same rate constants (`two_point.py`, 2026-10-03)

Exact criterion: c . ln x is the same at every positive steady state for every k > 0 iff there is no r > 0 in
ker N and dxi with N (r o exp(G dxi)) = 0 and c . dxi != 0 (G = S^T Y^T; detached uptakes get a free rate
ratio). For fixed dxi this is an LP (`realise`, exact version `exact_realise` for integer G and rational x');
`search` looks for dxi by least squares; `certify` is the sign relaxation as a MILP, optionally with the
iterated route's fixed log-combinations (V dxi = 0) and one sign per kinetic module. The tangent space is the
first-order version of the criterion, the sign MILP its sign relaxation.

Small networks (tests in `test_tangent_one_shot.py`): former completeness toy, C and P disproven exactly (x' = (1, 1, 2, 1/2)),
and found by the search without the known steady states; A and B proven. Figure 1: G proven only with V;
D, E and the six ACRR metabolites individually disproven. EnvZ-OmpR Yp: no counterexample, proven only with V.
The iHN637 witness (7 reactions): all three metabolites proven by the MILP.

iHN637 gap (`two_point_gap.py`, results `two_point_iHN637.json`): 2 ACR + 12 ACRR links, all undecided.
The MILP found neither a sign pattern nor infeasibility in 3 h per relation, and the search does not work at
genome scale (it fails even on metabolites the tangent route leaves free, where counterexamples exist).
