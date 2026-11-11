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
pip install "straindesign>=1.19.2"
pip install .
```

## Try it

[`examples/kimonu_examples.ipynb`](examples/kimonu_examples.ipynb) runs every call on the networks of Figures 1
and 2 of the paper: concordance, kinetic modules, ACR and ACRR by both routes, their exact values at given rate
constants, the rate constants for given values, and the values as formulas in the rate constants. The Binder badge
opens it in the browser, with nothing to install.

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
