"""Kimonu — exact, deterministic concordance and kinetic-module analysis of metabolic networks.

Pure Python. Depends on straindesign (exact sparse null space) and cobra.

Core idea: concordance is flux coupling lifted one matrix product up. With N = Y·A,
K = ker(N) (exact, sparse, via straindesign), and M = A·K, concordant complexes are the
proportional rows of M and balanced complexes are its zero rows — exact linear algebra,
no LP, no sampling.

Kinetic modules, ACR and ACRR are fixed log-combinations of the concentrations, read off one exact
basis by acr.readout. Two routes produce that basis: the iterated route (kinetic_modules, then
acr_acrr) proves relations at every steady state; the tangent route (tangent_route) reads the tangent
space of the steady states at one generic point, finds every relation that holds locally, and so
bounds the iterated route from above. relation_values gives the values those relations take at given
rate constants, exactly, and rate_constants the rate constants at which they take given values;
binomial_formulas and symbolic_formulas give the values as formulas in the rate constants.
"""
from kimonu.io import build_crn, irreversible
from kimonu.concordance import concordance, concordance_from_matrices, ConcordanceResult
from kimonu.kinetic import kinetic_modules, KineticResult
from kimonu.acr import acr_acrr, AcrResult, readout, Readout
from kimonu.tangent import tangent_route, TangentResult
from kimonu.values import relation_values, rate_constants, Values, Value
from kimonu.formulas import binomial_formulas, symbolic_formulas, Formulas, Monomial

__version__ = "0.0.1"
__all__ = ["build_crn", "irreversible", "concordance", "concordance_from_matrices", "ConcordanceResult",
           "kinetic_modules", "KineticResult", "acr_acrr", "AcrResult", "readout", "Readout",
           "tangent_route", "TangentResult", "relation_values", "rate_constants", "Values", "Value",
           "binomial_formulas", "symbolic_formulas", "Formulas", "Monomial", "__version__"]
