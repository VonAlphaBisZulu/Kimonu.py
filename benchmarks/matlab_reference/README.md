# Running the MATLAB reference implementation without MATLAB

The reference implementation for kinetic modules is [ankueken/Upstream_Algorithm](https://github.com/ankueken/Upstream_Algorithm)
(Küken et al. 2021, 2022; Langary et al. 2023, 2025). It is written for MATLAB and calls `linprog`
from the Optimization Toolbox. This directory runs it under **GNU Octave** instead, so the
reference can be included in the benchmark without a licence.

The reference's own `.m` files are used unmodified. `compat/` supplies only what Octave lacks:

| file | why |
|---|---|
| `compat/linprog.m` | Octave has no `linprog`; this maps the MATLAB calling convention onto the bundled GLPK, including the `optimset('linprog')` defaults query |
| `compat/round.m` | Octave's built-in `round` takes no digit count, and the concordance test calls `round(x, 2)` |
| `compat/triplets_to_sparse.m`, `compat/write_results_json.m` | driver helpers, not part of the reference |

## Usage

```bash
export UPSTREAM_ALGORITHM=/path/to/Upstream_Algorithm     # defaults to ../../../Upstream_Algorithm
octave --no-gui benchmarks/matlab_reference/run_reference.m \
    benchmarks/processed/concordance/e_coli_core.crn.json \
    benchmarks/results/reference_concordance_i_e_coli_core.json
```

The input is the shared exported network (`benchmarks/export_crn.jl`), so the reference decides
concordance on the *same* complex decomposition as Kimonu and COCOA. What is being compared is
therefore the decision procedure, not three separate complex enumerations.

Two deliberate differences from the reference's own workflow, both of which make the comparison
stricter rather than looser:

- **No sampling prefilter.** `pre_concordant_sampling` is stochastic and only ever removes
  candidate pairs. The driver tests every pair, which gives the reference formulation's exact
  answer and makes the run reproducible.
- **No `run_preprocessing`.** The model arrives already preprocessed, so COBRA Toolbox and the R
  stage (`get_AY_matrix.r`, which needs a pinned legacy igraph/Matrix environment) are not needed.

## Diagnostics

- `diagnose_lp.m <crn.json> <numerator> <denominator>` — prints the activity-ratio bounds the
  reference derives for one ordered pair, in both the positive and negative branch.
- `audit_missed.m <crn.json> <missed_pairs.json>` — for a list of pairs, reports the sign class of
  both complexes and whether the reference rejected the pair because an LP was infeasible or
  because the bounds genuinely differed.

The adaptations these led to are in pull request 1 of the fork
(https://github.com/ankueken/Upstream_Algorithm/pull/1) and are described in Supplementary Text S3.
