Four-phase (S.2.2) vs two-phase (Remark S2-1) upstream sets, and the module merge rule.

`measure_phases.py` runs both algorithms seed by seed (concordance module ∪ balanced) on Kimonu's
topology and checks autonomy, feeding and Lemma S2-4 on the final modules with balanced-only
linkage classes excluded. `giant_merge.py` reports the largest module under three rules: four-phase
with overlap merge, two-phase with overlap merge, two-phase with the Proposition S3-4 coupling
merge (Kimonu's default). Logs are the 2026-09-21 runs over all 33 models.

Result: the two algorithms return identical upstream sets for every seed on all 33 models; every
giant difference between Kimonu and the reference is the S3-4 merge, which the reference's R
kinetic script does not perform.
