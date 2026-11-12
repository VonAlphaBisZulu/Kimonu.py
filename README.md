# Kimonu.py
**KI**netic **MO**dules by **NU**llspace

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/VonAlphaBisZulu/Kimonu.py/main?labpath=examples%2Fkimonu_examples.ipynb)

Kimonu reads the steady-state invariants of a metabolic network from exact nullspaces. These are flux coupling,
balanced and concordant complexes, and, under mass-action kinetics, kinetic modules, absolute concentration robustness
(ACR) and absolute concentration ratio robustness (ACRR). Every nullspace is computed in exact rational arithmetic
with the sparse RREF of [StrainDesign](https://github.com/klamt-lab/straindesign), so zero, equal and proportional
rows are found without a tolerance.

The mass-action invariants come from two routes, which bound the true relations from both sides.

- The **iterated route** (the default) alternates the mass balance at a reference steady state with the comparison of
  complex compositions and feeds every proven coupling back. Every relation it reports holds at every positive
  steady state, but it can miss relations.
- The **tangent route** linearizes mass action at the reference steady state. It misses no relation, but it can
  report false couplings, which the linearization suggests but which do not hold at every positive steady state.

Where both routes agree, the list of relations is exact. Uptakes are detached by default, so a relation has to hold
for every supply, as in the published tools. Attached uptakes, which fix the medium, are an option.

## Install

Python 3.10 or later.

```
pip install git+https://github.com/VonAlphaBisZulu/Kimonu.py
```

This installs Kimonu from GitHub together with StrainDesign, COBRApy, NumPy and SciPy. From a clone of the
repository, `pip install .` does the same.

## Try it

[`examples/kimonu_examples.ipynb`](examples/kimonu_examples.ipynb) runs every call on the networks of Figures 1
and 2 of the paper: concordance, kinetic modules, ACR and ACRR by both routes, their exact values at given rate
constants, the rate constants for given values, and the monomial formulas in the rate constants that flux coupling
gives. The Binder badge opens it in the browser, with nothing to install.

## Use

The example runs every step on `e_coli_core`, which ships with COBRApy under the name `textbook`, so nothing has to
be downloaded.

```python
import cobra, kimonu
model = cobra.io.load_model("textbook")     # e_coli_core, or cobra.io.read_sbml_model("network.xml")
crn = kimonu.prepare(model)          # reversible reactions split, blocked ones removed, exact N = Y·A
conc = kimonu.concordance(crn)       # K, M, balanced complexes and concordance modules
it = kimonu.iterated_route(crn)      # proven kinetic modules, ACR and ACRR
tan = kimonu.tangent_route(crn)      # every true relation, false couplings possible
cert = kimonu.certify(crn)           # both routes at one r*: certified, proven lists, open relations
k = kimonu.rate_constants(crn, it.flux)     # rate constants at which x = 1 is a steady state, or your own
vals = kimonu.relation_values(crn, k)       # exact values of the relations at these rate constants

modules = [m for m in cert.modules if len(m) > 1]                           # a kinetic module has two or more complexes
gap = len(cert.open_pairs) + len(cert.open_acr) + len(cert.open_acrr)       # relations only the tangent route reports
share = sum(len(m) for m in modules) / len(crn.complexes)                   # the zero complex counts as a complex
print(f"certified {cert.certified}, {len(modules)} kinetic modules, gap {gap}, "
      f"{share:.1%} of {len(crn.complexes)} complexes in kinetic modules")
```

```
certified True, 4 kinetic modules, gap 0, 5.5% of 145 complexes in kinetic modules
```

A gap of 0 means that both routes report the same relations, so the list is complete. `kimonu.analyze(model)` is the
default in one call, the iterated route with detached uptakes, and `kimonu.analyze(model, route="both")` returns the
same certificate as `certify`.

For a network that already comes as matrices, `crn = kimonu.crn_from_matrices(A, Y)` builds the same `CRN`, with
N = Y·A computed exactly. Every function that reads a network takes the `CRN` first, `kimonu.iterated_route(crn)`,
`kimonu.concordance(crn)`, `kimonu.readout(crn, B)`. Every matrix Kimonu computes (K, M, W, B and Yᵀ B) is a list of
exact sparse rows, documented in `kimonu.linalg`.

## API

The functions follow the paper. Every one has a docstring that gives its purpose, its steps and the equation it
computes.

| paper | function | module |
|---|---|---|
| reversible reactions split, N = Y·A with the zero complex | `irreversible`, `build_crn` (returns a `CRN`) | `io` |
| blocked reactions (one LP), pinned rates | `prepare` | `io` |
| a network given as A and Y | `crn_from_matrices` | `io` |
| K = ker(N), M = A·K, balanced and concordant complexes (Fig. 1C, E) | `concordance` | `concordance` |
| ACR, ACRR and kinetic modules as zero and equal rows of B (Fig. 1F) | `readout`, `RouteResult` | `readout` |
| exact interior reference flux r* | `interior_flux` (an LP point by default, or draw "uniform" from 2^bits integers) | `steady_state` |
| W = N diag(r*) Sᵀ (Eq. 3) | `steady_state_map` | `steady_state` |
| iterated route, ker((YΔ)ᵀ) (Eq. 4, left) | `iterated_route`, `composition_differences`, `free_directions` | `iterated`, `readout` |
| tangent route, ker(W Yᵀ) (Eq. 4, right) | `tangent_route`, `tangent_space` | `tangent` |
| both routes, certified or open relations | `certify` (returns a `Certificate`) | `certify` |
| default mode from a cobra model | `analyze` (route "iterated", "tangent", "reference" or "both") | `analysis` |
| reference-equivalent mode (upstream sets, Shinar–Feinberg rule) | `reference_route` | `reference` |
| flux coupling under mass action, Λ = Sᵀ Yᵀ B | `kinetic_flux_coupling` | `coupling` |
| beyond the main text: values at given rate constants, rate constants for given values, monomial formulas | `relation_values`, `rate_constants`, `binomial_formulas` | `values` |

The three routes return one `RouteResult` with the modules, ACR, ACRR, the basis B they were read from, the
reference flux r* and the number of rounds. With detached uptakes the zero complex is a module of its own, listed
last, so the modules partition the complexes. `kimonu.linalg` holds the exact linear algebra they share.

## Tests

```
pip install ".[dev]"
pytest
```

## Paper and benchmarks

The paper is "Exact nullspaces find and certify kinetic modules and other mass-action invariants in metabolic
networks". [`benchmarks/README.md`](benchmarks/README.md) maps each of its figures, tables and numbers to the script
that made it, and the Zenodo archive holds the results (`benchmarks/ZENODO_README.md`).

## License

See [`LICENSE`](LICENSE).
