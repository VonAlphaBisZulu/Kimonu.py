# Benchmarks — Kimonu vs COCOA.jl vs the MATLAB reference

Correctness and speed comparison against COCOA.jl and against the MATLAB reference implementation
(ankueken/Upstream_Algorithm), mirroring the benchmark COCOA.jl ran against that same predecessor —
with Kimonu in the challenger seat. The reference runs under Octave and needs no MATLAB licence;
see `matlab_reference/`.

## Design

**One preprocessed model, all tools.** The Julia side runs COCOA's own preprocessing pipeline and
writes the resulting network to `processed/<profile>/<model>.json`; the Python side reads exactly
that file and runs with `preprocess=False`. Complex counts depend on which preprocessing steps ran,
so tools that preprocess separately are not comparable — this removes the question entirely.

**Complexes are matched by composition, not by name.** Neither tool's complex identifiers mean
anything to the other. The comparison keys every complex by its `{metabolite: coefficient}` map.

**Partitions are compared, not counts.** Two runs can agree on the number of modules and still
disagree about which complexes belong together. `compare.py` compares the induced set partitions,
and when they differ it reports the *direction* — whether one tool's concordant pairs are a strict
subset of the other's (missed concordances) or a strict superset (invented ones) — plus the
offending complexes.

**Singletons are excluded from module and pair counts.** COCOA labels a complex concordant with
nothing as `-1`; Kimonu returns it as a one-element group. Counting them differently manufactures
disagreements that are not there.

## Running

```bash
# 1. Export the preprocessed network and COCOA's complex decomposition (A, Y, S, bounds)
julia --project=benchmarks/julia benchmarks/export_crn.jl

# 2. COCOA: run the LP sweep, record the partition
julia --project=benchmarks/julia benchmarks/run_cocoa.jl

# 3. Kimonu: same model, exact nullspace
python benchmarks/run_kimonu.py

# 4. MATLAB reference, under Octave (optional; see matlab_reference/README.md)
export UPSTREAM_ALGORITHM=/path/to/Upstream_Algorithm
octave --no-gui benchmarks/matlab_reference/run_reference.m \
    benchmarks/processed/concordance/e_coli_core.crn.json \
    benchmarks/results/reference_concordance_i_e_coli_core.json

# 5. Compare
python benchmarks/compare.py
```

All three tools consume the same exported network, and the reference and Kimonu additionally use
COCOA's complex decomposition (`export_crn.jl` writes A and Y), so complex enumeration is held
fixed and only the decision procedure varies.

Both runners take `--model M` (repeatable), `--profile`, and `--repeats`. Timings report the
fastest of `--repeats` runs, which the Julia side needs so that JIT compilation of the first call
is not charged to COCOA.

`benchmarks/julia` is a separate Julia project that `dev`s the COCOA checkout. Point it at whichever
COCOA tree should be measured — benchmark against upstream `main`, not a working branch.

### Profiles

| profile | preprocessing |
|---|---|
| `concordance` | `normalize_bounds`, `remove_blocked_reactions`, `remove_orphans`, `split_into_irreversible` |
| `kinetic` | the above plus `split_into_elementary` (COCOA's recommendation for kinetic-module analysis) |

The `kinetic` profile is what reproduces the complex counts in the Julia proof-of-concept. Kimonu
has no elementary-step splitter of its own yet (roadmap Phase 1b); until it does, the split model
has to come from the Julia side, which is fine for benchmarking but not for standalone use.

### Scenarios

`i` is the full steady-state cone. `ii` additionally bounds the objective to within a relative
tolerance of its optimum. Only the COCOA runner accepts `--scenario ii`; Kimonu's scenario (ii)
needs the kernel restricted to reactions left free under the bound, which is not implemented.

## Caveats

- COCOA's `seed` defaults to `rand(UInt)`, so its sampling-based candidate filter is not
  reproducible unless the seed is pinned. The runner pins it and records it.
- COCOA is run in-process (`workers = Int[]`). It is designed to solve LPs in parallel, so this
  measures single-worker performance on both sides rather than COCOA's best achievable wall-clock.
- The models here are small. Nothing in these numbers speaks to genome-scale behaviour, which is
  where COCOA's own scaling results live.
