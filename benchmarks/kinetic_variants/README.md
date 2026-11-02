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

## Second order on the gap models (`second_order_modp.py`, 2026-10-28)

Tests whether the tangent-only relations break at second order along the steady states through z* (J q =
-N (r o (G h)^2), relation v breaks iff v^T q != 0), exactly over GF(2^31 - 1) at a random r in ker N, 4 random
h in ker J. Results in `second_order_results/`. All tangent-only relations give v^T q = 0 on iCN718, iHN637 and
iNJ661 (ACR 4/2/0, ACRR 390/1198/231 pairs, complex pairs 6391/5057/1398), and all stay orthogonal to ker J at
the new r (the first-order claim at a second, independent flux). Controls: proven relations 0, random row-space
vectors nonzero 20/20. As expected from the argument in src/kimonu/tangent.py: a fixed pair relation
orthogonal to ker J(r) at generic r is orthogonal at every interior r, so it is constant along every path of
steady states, to every order. The gap can only close globally (separate steady states at the same k).
iHN637: rank J = 2919 < rank W = 2920 at generic r (structurally degenerate), the second-order system was still
consistent. The float version (`second_order.py`) is not usable: at r* scale 10^6 lsqr stalls and the proven
controls come out larger than the candidates.

## Attached against detached uptakes (`uptake_modes.py`, `uptake_lists.py`, 2026-10-30)

Both routes from src, both uptake modes, all 33 models (`uptake_modes/<model>.json`; the three slowest models
also as `<model>.detached.json`). The attached block reproduces `hierarchy/<model>.json` and the gap numbers of
the draft (relations 31/6/14, g_rel 21/3/7 on iCN718/iNJ661/iHN637). Detached gives the same three gap models
with the same relations and the same g_rel. `uptake_lists.py` (iterated route, lists in
/scratch/sadreyev_lab/pschneider/kimonu_uptake_modes/lists): over the 33 models, free-metabolite ACR drops from
265 (attached) to 16 (detached); 150 of the 249 lost are extracellular and 18 periplasmic; 132 of the 249 stay
in detached ACRR groups. No metabolite gains ACR by detaching (detached within attached, as the theory says:
the detached steady states for given rate constants are the union of the attached ones over all supplies).

## The gap hangs on a few root relations (`gap_roots.py`, 2026-10-30)

Assume some open relations hold (roots) and feed them back: an ACR metabolite is dropped from every complex
(its constant goes into the rate constants), an ACRR pair a/b is merged into one metabolite (balances summed).
The original steady states are then steady states of the reduced network, and the iterated route holds for
every k there, so whatever it proves on the reduced network holds on the original once the roots are proven.
Attached uptakes, relations of `two_point_gap.py prepare`.
- iHN637: roots glyc3p and ser__L (the two open ACR) close all 14 relations.
- iNJ661: any single one of the 6 relations, as root, closes all 6.
- iCN718: nadph, malACP, acACP, propacp (the four open ACR) close 30 of 31; adding E_89/E231_ACP_accoa closes
  31. Single roots close at most 8 (nadph), pairs with nadph at most 13.
- Controls: no root closes 0 (rebuild reproduces the original, all reactions kept); four random free-metabolite
  pairs of similar degree (iHN637) close 0 of 14.
The gap is therefore 8 root questions (2 + 1 + 5) rather than 51 relations; a certificate only has to settle
the roots.

## Local sign certificate with lazy rows (`two_point_local.py`, 2026-10-30)

The sign MILP of `two_point.certify` (with V and module classes) on a subset of balance rows, grown from the
target's rows by the reactions that move in each solution (hubs last), 900 s per round, at most 8 rounds or
1500 rows. Dropping rows only relaxes the MILP, so an infeasible subsystem would prove the relation. Toys behave
as expected (Figure 2 kinase X proven, C and P not; completeness toy A, B proven, C, P not; Figure 1 G proven
with V, D and E not). Gap models (attached, job 568337, outputs in
/scratch/sadreyev_lab/pschneider/kimonu_twopoint_local/): 0 of 51 proven. Every subsystem had a sign pattern
while small (up to 95-865 rows), and the MILP hit 900 s at 500-1434 rows; iCN718 link E_89/E231_ACP_accoa
reached the row cap with sign patterns throughout. The local sign patterns mostly extend to N d = 0 on the whole
network. No small local proof exists in these neighbourhoods; the obstruction, if any, is global.

## Relation values at given k (`relation_values.py`, 2026-11-02)

kimonu.relation_values on the 33 models, attached and detached (job 568428, outputs in
/scratch/sadreyev_lab/pschneider/kimonu_values/). Round trip: x random small rationals (1 on metabolites with a
fractional substrate coefficient), k = r* / ψ(x) by kimonu.rate_constants with the integer r* (scale 10^6). Every
ACR value came back as x_i, every ACRR ratio as x_j / x_i and every module ratio as ψ_d(x) / ψ_c(x), exactly, in
all 66 runs, and the relation counts equal the iterated route's (uptake_modes/). At most 57 s, 1 to 4 rounds;
every feedback ratio was rational. At prime rate constants (k_j = a permuted j-th prime), 20 models attached and
12 detached have no positive steady state (a complex forced to zero): random k rarely admits one at genome scale.

## Monomial formulas from flux coupling (`binomial_formulas.py`, 2026-11-02)

kimonu.binomial_formulas on the 33 models (job 568461, outputs in /scratch/sadreyev_lab/pschneider/kimonu_binomial/),
each formula evaluated at the k of `relation_values.py` and compared with the exact value there: 0 wrong, at most
78 s. Values with a monomial formula (each against the first member of its group): detached ACR 84/93, ACRR
3745/4945, module ratios 4253/24566; attached ACR 176/860, ACRR 3221/4502, module ratios 4390/25272. Up to 18 rate
constants in one formula. Flux coupling constrains k in 9 models attached and 8 detached (iSB619 and iYO844 5
constraints each); the k of rate_constants satisfies all of them.
