# Kimonu.py
**KI**netic **MO**dules by **NU**llspace

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/VonAlphaBisZulu/Kimonu.py/main?labpath=examples%2Fkimonu_examples.ipynb)

Kimonu reads the steady-state invariants of a metabolic network from exact nullspaces. These are flux coupling,
balanced and concordant complexes, and, under mass-action kinetics, kinetic modules, absolute concentration robustness
(ACR) and absolute concentration ratio robustness (ACRR), each with the exact ratio or value it fixes. Every nullspace is computed in exact rational arithmetic
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
be downloaded. It is [`examples/quickstart.py`](examples/quickstart.py).

```python
import cobra, kimonu
model = cobra.io.load_model("textbook")     # e_coli_core, or cobra.io.read_sbml_model("network.xml")
crn = kimonu.preprocess(model)       # reversible reactions split, blocked ones removed, exact N = Y·A
conc = kimonu.concordance(crn)       # K, M, balanced complexes, concordance modules and activity ratios
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

A gap of 0 means that both routes report the same relations, so the list is complete. Every function also takes the
cobra model itself and preprocesses it first, as in `kimonu.iterated_route(model)`. Each such call preprocesses the
model anew, so pass the `CRN` when calling several functions. Uptakes are detached by default, and
`detach_uptakes=False` attaches them in every call that takes this keyword.

For a network that already comes as matrices, `crn = kimonu.crn_from_matrices(A, Y)` builds the same `CRN`, with
N = Y·A computed exactly. Every function that reads a network takes it first, `kimonu.iterated_route(crn)`,
`kimonu.concordance(crn)`, `kimonu.readout(crn, B)`. Every matrix Kimonu computes (K, M, W, B and Yᵀ B) is a list of
exact sparse rows, documented in `kimonu.linalg`.

### Groups and their ratios

Each coupling comes with the ratio it fixes. Fully coupled reactions keep r_j = λ r_i at every steady state (the
proportional rows of K), and concordant complexes keep a_c = λ a_d for their activities a = A·r (the proportional rows
of M). An ACRR group of k metabolites holds k − 1 independent ratios. Continuing the example above:

```python
fc = kimonu.mass_balance_flux_coupling(crn)     # fully coupled reactions and their flux ratios λ
for g in fc.groups:
    if len(g) > 1:
        print("fully coupled:", {crn.reaction_ids[j]: str(fc.ratio[j]) for j in sorted(g)})
name = lambda c: " + ".join(m if v == 1 else f"{v} {m}" for m, v in crn.complexes[c]) or "0"
o2 = crn.complexes.index((("o2_c", 1),))        # the complex o2_c, formed by O2 transport
module = next(m for m in conc.modules if o2 in m)
print("concordant:", {name(c): str(conc.activity_ratio[c]) for c in module})
print("ACRR groups:", it.acrr_groups, "independent ratios:", sum(len(g) - 1 for g in it.acrr_groups))
```

```
fully coupled: {'EX_for_e': '1', 'PFL': '1'}
fully coupled: {'GND': '1', 'PGL': '1'}
fully coupled: {'ICL': '1', 'MALS': '1'}
concordant: {'2 h_c + 1/2 o2_c + q8h2_c': '1', 'h2o_c + 2 h_e + q8_c': '-1', 'o2_c': '-1/2'}
ACRR groups: [['acald_c', 'acald_e']] independent ratios: 1
```

The module of o2_c holds the substrate complex of cytochrome bd oxidase, which takes up half an O2, so the activity
of o2_c is −1/2 that of this complex.

### Relations at your own rate constants

The routes find the relations for generic rate constants. `kimonu.relation_values(crn, k)` reruns the iterated route
at given rate constants k itself and returns the kinetic modules at that k, each with the exact ratios ψ_c / ψ_d of
its complex concentrations, and the exact values of ACR (x_i) and ACRR (x_j / x_i). Values are rational, or a
rational times a root, and are exact either way. `kimonu.rate_constants(crn, flux, x)` gives the rate constants for
which chosen concentrations x are a steady state with a chosen interior flux, so every relation takes the value that
x gives it. Continuing the example above:

```python
x = [1] * len(crn.metabolite_ids)                   # the concentrations to prescribe, 1 unless set
x[crn.metabolite_ids.index("acald_c")] = 3
x[crn.metabolite_ids.index("for_e")] = 2
k = kimonu.rate_constants(crn, it.flux, x)          # x is a steady state for k, with the flux it.flux
vals = kimonu.relation_values(crn, k)               # the iterated route at k itself, with exact values
print("ACR:", vals.acr, " ACRR:", vals.acrr)
for first, ratios in vals.modules:                  # ψ_c / ψ_first within each kinetic module at k
    print(name(first), "|", {name(c): value for c, value in ratios.items()})
```

```
ACR: {}  ACRR: [('acald_c', {'acald_e': 1/3})]
acald_e | {'acald_c': 3}
icit_c | {'accoa_c + glx_c + h2o_c': 1}
for_e | {'coa_c + pyr_c': 1/2}
6pgc_c + nadp_c | {'6pgl_c + h2o_c': 1}
```

`kimonu.binomial_formulas(crn)` gives, where flux coupling fixes a relation, its value as a monomial formula in k.
The [notebook](examples/kimonu_examples.ipynb) runs all three on the networks of Figures 1 and 2.

## API

The functions follow the paper. Every one has a docstring that gives its purpose, its steps and the equation it
computes.

| paper | function | module |
|---|---|---|
| reversible reactions split, N = Y·A with the zero complex | `make_irreversible`, `build_crn` (returns a `CRN`) | `network` |
| blocked reactions (one LP), pinned rates | `preprocess` | `network` |
| a network given as A and Y | `crn_from_matrices` | `network` |
| a `CRN` or a cobra model, as every function takes it | `as_crn` | `network` |
| K = ker(N), fully coupled reactions and their flux ratios (Fig. 1C) | `mass_balance_flux_coupling` | `coupling` |
| M = A·K, balanced and concordant complexes and their activity ratios (Fig. 1E) | `concordance` | `concordance` |
| zero rows, groups of proportional rows and their ratios, as every readout of K, M and ker(W) uses them | `row_readout` | `readout` |
| ACR, ACRR and kinetic modules as zero and equal rows of B (Fig. 1F) | `readout`, `RouteResult` | `readout` |
| exact interior reference flux r* | `interior_flux` (an LP point by default, or draw "uniform" from 2^bits integers) | `steady_state` |
| W = N diag(r*) Sᵀ (Eq. 3) | `steady_state_map` | `steady_state` |
| iterated route, ker((YΔ)ᵀ) (Eq. 4, left) | `iterated_route`, `composition_differences`, `free_directions` | `iterated`, `readout` |
| tangent route, ker(W Yᵀ) (Eq. 4, right) | `tangent_route`, `tangent_space` | `tangent` |
| both routes, certified or open relations | `certify` (returns a `Certificate`) | `certify` |
| reference-equivalent mode (upstream sets, Shinar–Feinberg rule) | `reference_route` | `reference` |
| flux coupling under mass action, Λ = Sᵀ Yᵀ B | `kinetic_flux_coupling` | `coupling` |
| beyond the main text: values at given rate constants, rate constants for given values, monomial formulas | `relation_values`, `rate_constants`, `binomial_formulas` | `relation_values` |

The three routes return one `RouteResult` with the modules, ACR, ACRR and the ACRR groups, the basis B they were read
from, the reference flux r* and the number of rounds. With detached uptakes the zero complex is a module of its own,
listed last, so the modules partition the complexes. `kimonu.linalg` holds the exact linear algebra they share.

## Tests

```
pip install ".[dev]"
pytest
```

## Paper and scripts

The paper is "Exact nullspaces find and certify kinetic modules and other mass-action invariants in metabolic
networks". [`scripts/README.md`](scripts/README.md) maps each of its figures, tables and numbers to the script
that made it, and the Zenodo archive holds the results (`scripts/ZENODO_README.md`). No results are kept in git.
The scripts read one records folder in the layout of the archive, named by the environment variable `KIMONU_RECORDS`
(default `records/` at the top of the checkout), so the unpacked archive is all they need. They write to
`KIMONU_OUT` (default the same folder), which keeps the archive as it is when set to another folder.

## License

See [`LICENSE`](LICENSE).
