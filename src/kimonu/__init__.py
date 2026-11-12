"""Kimonu (KInetic MOdules by NUllspace): exact concordance, kinetic modules, ACR and ACRR of metabolic networks.

Pure Python. Depends on straindesign (exact sparse nullspace), cobra, numpy and scipy. Every nullspace is computed
in rational arithmetic, so zero, equal and proportional rows are found by exact comparison, without tolerance.

The default mode, one call from a cobra model:

    import kimonu
    crn = kimonu.prepare(model)                 # split, unblocked, N = Y·A
    res = kimonu.analyze(crn)                   # iterated route, detached uptakes, seed 0
    res.modules, res.acr, res.acrr              # complex indices into crn.complexes, metabolite ids

Conventions
    The network.  A network is a :class:`CRN`, from :func:`prepare` for a cobra model or from
        :func:`crn_from_matrices` for a network given as A and Y. Every function that reads a network takes the CRN
        as its first argument and nothing else of it, ``iterated_route(crn)``, ``concordance(crn)``,
        ``readout(crn, B)``. :func:`analyze` and :func:`concordance` also take a cobra model and prepare it first.
    Matrices.  Every matrix Kimonu computes, K, M, W, Yᵀ, (YΔ)ᵀ, B and Yᵀ B, is in the row form of
        :mod:`kimonu.linalg`, a list with one row per matrix row, each row a sorted tuple of (column, value) pairs
        with exact values. Equal rows are equal tuples, so every readout groups rows in a dictionary. A is a scipy
        integer matrix and N a straindesign RationalMatrix, the input of every nullspace of N.
    Results.  Complexes and reactions are indices into crn.complexes and crn.reaction_ids. Metabolites in results
        are named by crn.metabolite_ids, which are indices for a CRN built from matrices without names. Arguments
        that list metabolites, such as ``known_acr``, take indices.

Paper concept -> function, in the order of the paper (Schneider, Exact nullspaces find and certify
kinetic modules and other mass-action invariants in metabolic networks).

Networks and preprocessing
    split reversible reactions                     irreversible(model)
    N = Y·A, the zero complex as a complex         build_crn(model) -> CRN (A, Yt, N exact, ids)
    blocked reactions (one LP), pinned rates       prepare(model) -> CRN
    a network given as matrices                    crn_from_matrices(A, Y) -> CRN
Flux balance (Fig. 1C, E)
    K = ker(N), fully coupled reactions            concordance(crn).K, kimonu.linalg.proportional_classes
    M = A·K, balanced and concordant complexes     concordance(crn) -> ConcordanceResult
Mass action, the basis B (Fig. 1F)
    ACR, ACRR, kinetic modules from B              readout(crn, B) -> Readout, and RouteResult of every route
    Yᵀ B                                           complex_directions(crn, B)
Reference steady state (Eq. 3, Fig. 1G)
    exact interior flux r*, generic k = r*         interior_flux(crn, seed, draw="lp", bits=64)
    W = N diag(r*) Sᵀ                              steady_state_map(crn, r*)
Iterated route (Eq. 4, left)
    coupled sets of ker(W), YΔ, ker((YΔ)ᵀ), feedback   iterated_route(crn) -> RouteResult, .B
    YΔ and ker((YΔ)ᵀ) of given modules             composition_differences(crn, modules), free_directions(crn, modules)
Tangent route (Eq. 4, right)
    ker(W Yᵀ), false couplings possible            tangent_route(crn) -> RouteResult, .B
    the full basis, with source coordinates        tangent_space(crn, r*)
Both routes (Eq. 4)
    agreement certifies the list, open relations   certify(crn) -> Certificate
Default mode and alternatives
    one call from a cobra model                    analyze(model, route="iterated" | "tangent" | "reference" | "both")
Reference-equivalent mode (Fig. 4E)
    upstream sets, overlap merge, Shinar–Feinberg rule   reference_route(crn)
Flux coupling under mass action
    Λ = Sᵀ Yᵀ B                                    kinetic_flux_coupling(crn, route=...)
Values of the relations (beyond the main text, see the Supplementary Text)
    exact values at given k                        relation_values(crn, k)
    k for given concentrations and fluxes          rate_constants(crn, r, x)
    monomial formulas in k from flux coupling      binomial_formulas(crn)

The three routes return one result type, :class:`RouteResult`, with the modules, ACR and ACRR, the basis B they
were read from, the flux r* and the route-specific extras.
"""
from kimonu.io import CRN, build_crn, crn_from_matrices, irreversible, prepare
from kimonu.concordance import ConcordanceResult, concordance
from kimonu.readout import (Readout, RouteResult, complex_directions, composition_differences, free_directions,
                            readout)
from kimonu.steady_state import interior_flux, steady_state_map
from kimonu.iterated import iterated_route
from kimonu.tangent import tangent_route, tangent_space
from kimonu.certify import Certificate, certify
from kimonu.analysis import analyze
from kimonu.reference import reference_route
from kimonu.coupling import FluxCouplingResult, kinetic_flux_coupling
from kimonu.values import Formulas, Monomial, Value, Values, binomial_formulas, rate_constants, relation_values

__version__ = "0.0.1"
__all__ = [
    # networks and preprocessing
    "irreversible", "build_crn", "prepare", "crn_from_matrices", "CRN",
    # flux balance
    "concordance", "ConcordanceResult",
    # mass action, the basis B and its readout
    "readout", "Readout", "RouteResult", "complex_directions", "composition_differences", "free_directions",
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
    # values of the relations
    "relation_values", "rate_constants", "binomial_formulas", "Values", "Value", "Formulas", "Monomial",
    "__version__",
]
