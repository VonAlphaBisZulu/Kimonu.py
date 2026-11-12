"""Kimonu (KInetic MOdules by NUllspace): exact concordance, kinetic modules, ACR and ACRR of metabolic networks.

Pure Python. Depends on straindesign (exact sparse nullspace), cobra, numpy and scipy. Every nullspace is computed
in rational arithmetic, so zero, equal and proportional rows are found by exact comparison, without tolerance.

The default mode, one call from a cobra model:

    import kimonu
    crn = kimonu.prepare(model)                 # split, unblocked, N = Y·A
    res = kimonu.analyze(crn)                   # iterated route, detached uptakes, seed 0
    res.modules, res.acr, res.acrr              # complex indices into crn.complexes, metabolite ids

Paper concept -> function, in the order of the paper (Schneider, Exact nullspaces find and certify
kinetic modules and other mass-action invariants in metabolic networks).

Networks and preprocessing
    split reversible reactions                     irreversible(model)
    N = Y·A, the zero complex as a complex         build_crn(model) -> CRN (A, Y, N exact, N_float, ids)
    blocked reactions (one LP), pinned rates       prepare(model) -> CRN
Flux balance (Fig. 1C, E)
    K = ker(N), fully coupled reactions            concordance(...).K, kimonu.linalg.proportional_classes
    M = A·K, balanced and concordant complexes     concordance(model), concordance_from_matrices(A, N)
Mass action, the basis B (Fig. 1F)
    ACR, ACRR, kinetic modules from B              readout(Y, B) -> Readout, and RouteResult of every route
    Yᵀ B                                           kimonu.readout.complex_directions(Y, B)
Reference steady state (Eq. 3, Fig. 1G)
    exact interior flux r*, generic k = r*         interior_flux(N, seed, draw="lp", bits=64)
    W = N diag(r*) Sᵀ                              steady_state_map(A, Y, N, r*)
Iterated route (Eq. 4, left)
    coupled sets of ker(W), YΔ, ker((YΔ)ᵀ), feedback   iterated_route(A, Y, N) -> RouteResult, .B
    ker((YΔ)ᵀ) of given modules                    free_directions(Y, modules)
Tangent route (Eq. 4, right)
    ker(W Yᵀ), false couplings possible            tangent_route(A, Y, N) -> RouteResult, .B
    the full basis, with source coordinates        tangent_space(A, Y, N, r*)
Both routes (Eq. 4)
    agreement certifies the list, open relations   certify(A, Y, N) -> Certificate
Default mode and alternatives
    one call from a cobra model                    analyze(model, route="iterated" | "tangent" | "reference" | "both")
Reference-equivalent mode (Fig. 4E)
    upstream sets, overlap merge, Shinar–Feinberg rule   reference_route(A, Y, N)
Flux coupling under mass action
    Λ = Sᵀ Yᵀ B                                    kinetic_flux_coupling(A, Y, N, route=...)
Values of the relations
    exact values at given k                        relation_values(A, Y, N, k)
    k for given concentrations and fluxes          rate_constants(A, Y, r, x)
    formulas in k                                  binomial_formulas(A, Y, N), symbolic_formulas(A, Y, N)

The three routes return one result type, :class:`RouteResult`, with the modules, ACR and ACRR (named when
metabolite_ids is given), the basis B they were read from, the flux r* and the route-specific extras.
"""
from kimonu.io import CRN, build_crn, irreversible, prepare
from kimonu.concordance import ConcordanceResult, concordance, concordance_from_matrices
from kimonu.readout import Readout, RouteResult, free_directions, readout
from kimonu.steady_state import interior_flux, steady_state_map
from kimonu.iterated import iterated_route
from kimonu.tangent import tangent_route, tangent_space
from kimonu.certify import Certificate, certify
from kimonu.analysis import analyze
from kimonu.reference import reference_route
from kimonu.coupling import FluxCouplingResult, kinetic_flux_coupling
from kimonu.values import Value, Values, rate_constants, relation_values
from kimonu.formulas import Formulas, Monomial, binomial_formulas, symbolic_formulas

__version__ = "0.0.1"
__all__ = [
    # networks and preprocessing
    "irreversible", "build_crn", "prepare", "CRN",
    # flux balance
    "concordance", "concordance_from_matrices", "ConcordanceResult",
    # mass action, the basis B and its readout
    "readout", "Readout", "RouteResult", "free_directions",
    # reference steady state
    "interior_flux", "steady_state_map",
    # iterated route
    "iterated_route",
    # tangent route
    "tangent_route", "tangent_space",
    # both routes
    "certify", "Certificate",
    # default mode and alternatives
    "analyze",
    # reference-equivalent mode
    "reference_route",
    # flux coupling under mass action
    "kinetic_flux_coupling", "FluxCouplingResult",
    # values and formulas
    "relation_values", "rate_constants", "Values", "Value",
    "binomial_formulas", "symbolic_formulas", "Formulas", "Monomial",
    "__version__",
]
