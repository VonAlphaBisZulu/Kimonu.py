# Kimonu.py — roadmap to a publication

*A Python package for exact, deterministic kinetic-module analysis of metabolic networks —
a Python twin of COCOA.jl, re-engineered around exact sparse nullspace computation.*

Status: concept validated in a Julia proof-of-concept (`exact_nullspace.jl`,
`concordance_nullspace.jl`); Python package not yet started. This document is the plan.

---

## 1. The idea in one line

**Concordance is flux coupling lifted one matrix-multiply up.** Compute the exact nullspace
`K = ker(N)` of the stoichiometric matrix `N = Y·A`, then group the proportional rows of
`M = A·K`. That single exact linear-algebra step replaces COCOA's per-pair Charnes–Cooper LP
sweep with a rank-1 proportionality test — **no LP, no sampling, no CV threshold, fully
deterministic.**

```
N : m×n stoichiometry (= Y·A)     Y : m×p complex composition     A : p×n incidence
K = ker(N)  (n×ν, exact rational, SPARSE)      M = A·K  (p×ν)
concordant complexes ⇔ proportional rows of M      balanced ⇔ zero rows of M
```

The kinetic-module layer (autonomy + feeding refinement) and ACR/ACRR (column-span in the
composition-difference matrix) sit on top of this concordance core and reuse `Y`.

## 2. Why it is worth publishing (selling points)

1. **Exactness / determinism.** Concordance/balance are fixed by exact linear structure; the LP
   solver is irrelevant. Kimonu computes them exactly over the rationals — no seed, no CV
   prefilter, reproducible bit-for-bit. This *fixes* real defects in the LP route (see §3).
2. **Speed.** Julia PoC: e_coli_core **~6 s vs ~150 s** for the CC reference (~25×), *and* exact.
   The Python target rides the same sparse exact nullspace.
3. **Accessibility.** Brings the kinetic-module framework to the Python/cobrapy ecosystem
   (COCOA is Julia; the original reference is MATLAB+R).

## 3. Correctness leverage (already established, see COCOA_review_notes / NULLSPACE_RESULTS)

The nullspace route is not just faster — it is demonstrably *more correct* than the objective-
bounded LP route, from the LP route's own numbers:

- COCOA's objective-bounded run reports **123 modules / 2243–2245 pairs** on e_coli, deterministic
  across seeds {1234, 999, 7}. But the bounded flux set `S₂ ⊆ S₁` can only make concordance
  *coarser* ⇒ `pairs(ii) ≥ pairs(i) = 2436`. **2245 < 2436 is impossible for an exact method** ⇒
  ~190 false-negative concordances (a Charnes–Cooper objective-scaling bug: the bound is added as
  a raw constant, not scaled by the homogenizing `t`).
- `efficient=true` ACR path misses true ACR species (e.g. Yp in EnvZ-OmpR).

Kimonu is exact in both regimes, so it is the natural "ground truth" implementation.

## 4. What already exists

- **Julia PoC** (`exact_nullspace.jl` = dep-free sparse integer RREF ported from straindesign;
  `concordance_nullspace.jl` = row-proportionality grouping). `N·K == 0` verified exactly.
- **Validation** (NULLSPACE_RESULTS.md): EnvZ-OmpR full exact match (13 complexes / 5 balanced /
  3 modules / 8 pairs); e_coli_core exact scenario-(i) match (895 / 209 / **110 / 2436**), ~25×.
- **Shared models** (`shared_models/`): EnvZ-OmpR, e_coli_core, anaerobic e_coli_core, iAB_RBC_283.
- **MATLAB reference** (`Upstream_Algorithm/`, ankueken) — fragile (needs a pinned legacy R env +
  a `cbind` fix; its committed workflow doesn't call the elementary-step splitter). Runs, with work.

## 5. Package architecture (proposed)

Dependency: **the lifted sparse-exact-nullspace package** (working name straindesign — see the
separate nullspace roadmap). Kimonu is a thin, exact analysis layer on top of it.

```
kimonu/
  io.py            load cobra/SBML models; build N, Y, A (complex enumeration = graph nodes)
  nullspace.py     thin wrapper over straindesign: K = ker(N), exact & sparse
  concordance.py   M = A·K; group proportional rows (exact rational); balanced = zero rows;
                   ratios γ; trivial concordances read off Y's two-complex columns
  kinetic.py       kinetic modules: autonomy + feeding refinement
  acr.py           ACR / ACRR via column-span of the composition-difference matrix (rigorous path)
  scenarios.py     (i) irreversibility-only ker(N); (ii) near-optimal: restrict K to FVA-free reacs
  api.py / cli.py  user entry points; mirror COCOA's public surface where it aids comparison
```

Design rules: exact rationals end-to-end (no float thresholds in the core); reproducible;
cobrapy-native I/O; keep the concordance core independent of the kinetic/ACR layers.

## 6. Benchmark plan (the paper is benchmark-heavy)

- **(a) Exactness vs COCOA.jl & MATLAB reference** — complexes, balanced, modules, concordant
  pairs, ACR/ACRR on the shared models. Lead with an *exact* toy result (EnvZ 13/5/3/8). Confirm
  all three consume the *same* preprocessed (elementary-split) model before comparing counts.
- **(b) Reproducibility** — Kimonu identical across runs by construction; contrast COCOA's
  seed/CV sensitivity and the bounded-mode false negatives (§3).
- **(c) Speed, same machine** — head-to-head vs COCOA.jl **and** the MATLAB reference on the
  shared models. (The review notes flag that COCOA's own figures measure only its scaling, not a
  head-to-head — this is the gap Kimonu's paper fills.) *MATLAB license pending; note if deferred.*
- **(d) Scaling + peak memory** — genome-scale networks; the sparse exact nullspace is the enabler.
  This is the shared benchmark with the straindesign package (speed + sparsity of `K`).
- **(e) Full suite parity** — reassemble COCOA's 62-run suite (33 ordered / 29 random binding).

## 7. Publication

- **Venue:** BMC Bioinformatics application note (or similar app-note track).
- **Working title:** *Harnessing sparse exact nullspace computation for fully deterministic
  detection of concordance and kinetic modules.*
- **Selling points:** Python accessibility · significant speedup · exactness/determinism.
- **Suggested reviewers:** Anika Küken, Zoran Nikoloski, Venkata Satagopam.
- **Positioning vs COCOA.jl:** independent, exact Python implementation; cite and credit COCOA and
  the original kinetic-module work; the correctness findings (CC scaling bug, efficient-mode ACR)
  are context and a possible shared correction, not a takedown.

## 8. TODO — phased

**Phase 0 — foundation**
- [ ] Lift the sparse exact nullspace out of straindesign into a standalone package (straindesign);
      Kimonu and straindesign both depend on it. (See the nullspace-package roadmap.)

**Phase 1 — concordance core in Python**  *(core done; straindesign-backed, tested)*
- [x] `io.build_crn`: complex enumeration (substrate/product multisets, COCOA-style order) → A, Y, N;
      verified `N == Y·A == cobra stoichiometric matrix`.
- [x] `concordance.concordance`: FVA-unblock (the only, one-time LP) → K=ker(N) via **straindesign** →
      M=A·K → exact proportional-row grouping (balanced = zero rows). Handles the straindesign ExactCOO
      (>int64) path exactly. Tests: N==Y·A, internal consistency, determinism.
- **Phase 1b — preprocessing to match the PoC counts.** Exact agreement with the Julia PoC (EnvZ
  13/5/3/8, e_coli 895/209/110/2436) requires the **elementary-step splitting** the PoC ran; the
  shared SBML are un-split (raw EnvZ → 16 complexes, not 13). This is the preprocessing-consistency
  issue the review notes flagged. Implement `split_into_elementary_rxns`, then:
  - [ ] Reproduce EnvZ (13/5/3/8) and e_coli_core (895/209/110/2436).
  - [ ] Küken Fig.1 toy unit test: balanced {2A, D, A+E, F}; module {B, C, B+C, 2B}; γ(C, B+C) = −1.
  - [ ] Resolve the **896-vs-895** complex enumeration off-by-one (boundary/empty-complex edge case).
- [ ] concordance ratios γ_ij and the trivial-concordance shortcut from Y (in the Julia PoC; not yet ported).

**Phase 2 — kinetic modules + ACR**
- [ ] `kinetic.py`: autonomy + feeding refinement of concordance modules.
- [ ] `acr.py`: ACR/ACRR via the rigorous column-span path (avoid efficient-mode incompleteness).
- [ ] Scenario (ii): restrict K to FVA-free reactions; test on a *realistically bounded* model
      (COCOA's normalize_bounds wipes uptake limits, so (ii) ≡ (i) there — need a model where the
      bound actually fixes reactions).

**Phase 3 — benchmarks** (§6 a–e), incl. same-machine head-to-heads and genome-scale scaling/memory.

**Phase 4 — packaging & paper**
- [ ] PyPI package, docs, CI; reproducibility scripts for every figure.
- [ ] Draft the app note.

## 9. Open questions / risks
- MATLAB license for the reference head-to-head (currently unavailable).
- Scenario-(ii) genuinely-different test needs a realistically-bounded model.
- Complex-enumeration parity (896 vs 895) must be nailed before count comparisons are clean.
- Kinetic-module + ACR layers are less validated than the concordance core — need their own toy tests.
