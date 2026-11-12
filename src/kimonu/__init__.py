"""Kimonu (KInetic MOdules by NUllspace): exact concordance, kinetic modules, ACR and ACRR of metabolic networks.

Pure Python. Depends on straindesign (exact sparse nullspace), cobra, numpy and scipy. Every nullspace is computed
in rational arithmetic, so zero, equal and proportional rows are found by exact comparison, without tolerance.

The default analysis, from a cobra model:

    import kimonu
    crn = kimonu.preprocess(model)              # split, unblocked, N = Y·A
    res = kimonu.iterated_route(crn)            # detached uptakes, seed 0
    res.modules, res.acr, res.acrr              # complex indices into crn.complexes, metabolite ids

Conventions
    The network.  A network is a :class:`CRN`, from :func:`preprocess` for a cobra model or from
        :func:`crn_from_matrices` for a network given as A and Y. Every function that reads a network takes it as
        its first argument, never its matrices one by one, as in ``iterated_route(crn)``, ``concordance(crn)`` and
        ``readout(crn, B)``. Each also takes a cobra model and preprocesses it first (:func:`as_crn`). Pass the CRN
        when calling several functions, so that preprocessing runs once.
    Matrices.  Every matrix Kimonu computes, K, M, W, Yᵀ, (YΔ)ᵀ, B and Yᵀ B, is in the row form of
        :mod:`kimonu.linalg`, a list with one row per matrix row, each row a sorted tuple of (column, value) pairs
        with exact values. Equal rows are equal tuples, so every readout groups rows in a dictionary
        (:func:`row_readout`, :func:`readout`). A is a scipy integer matrix and N a straindesign RationalMatrix, the
        input of every nullspace of N. Yᵀ B is ``kimonu.linalg.matrix_product(crn.Yt, B)``.
    Results.  Complexes and reactions are indices into crn.complexes and crn.reaction_ids. Metabolites in results
        are named by crn.metabolite_ids, which are indices for a CRN built from matrices without names. Arguments
        that list metabolites, such as ``known_acr``, take indices. Groups of reactions, complexes or metabolites
        are called groups throughout, as in ``acrr_groups``, ``FluxCouplingResult.groups`` and
        ``RowReadout.groups``.

Paper concept -> function, in the order of the paper (Schneider, Exact nullspaces find and certify
kinetic modules and other mass-action invariants in metabolic networks).

Networks and preprocessing (kimonu.network)
    split reversible reactions                     make_irreversible(model)
    N = Y·A, the zero complex as a complex         build_crn(model) -> CRN (A, Yt, N exact, ids)
    blocked reactions (one LP), pinned rates       preprocess(model) -> CRN
    a network given as matrices                    crn_from_matrices(A, Y) -> CRN
    a CRN or a cobra model, as every call takes it as_crn(network) -> CRN
Flux balance (Fig. 1C, E; kimonu.coupling, kimonu.concordance)
    K = ker(N), fully coupled reactions, r_j = λ r_i     mass_balance_flux_coupling(crn) -> MassBalanceCouplingResult
    M = A·K, balanced and concordant complexes, a_c = λ a_d   concordance(crn) -> ConcordanceResult
    zero rows, proportional rows and their ratio λ       row_readout(rows) -> RowReadout
Mass action, the basis B (Fig. 1F; kimonu.readout)
    ACR, ACRR, kinetic modules from B              readout(crn, B) -> Readout, and RouteResult of every route
Reference steady state (Eq. 3, Fig. 1G; kimonu.steady_state)
    exact interior flux r*, generic k = r*         interior_flux(crn, seed, draw="lp", bits=64)
    W = N diag(r*) Sᵀ                              steady_state_map(crn, r*)
Iterated route (Eq. 4, left; kimonu.iterated)
    coupled sets of ker(W), YΔ, ker((YΔ)ᵀ), feedback   iterated_route(crn) -> RouteResult, .B
    YΔ and ker((YΔ)ᵀ) of given modules             composition_differences(crn, modules), free_directions(crn, modules)
Tangent route (Eq. 4, right; kimonu.tangent)
    ker(W Yᵀ), false couplings possible            tangent_route(crn) -> RouteResult, .B
    the full basis, with source coordinates        tangent_space(crn, r*)
Both routes (Eq. 4; kimonu.certify)
    agreement certifies the list, open relations   certify(crn) -> Certificate
Reference-equivalent mode (Fig. 4E; kimonu.reference)
    upstream sets, overlap merge, Shinar–Feinberg rule   reference_route(crn)
Flux coupling under mass action (kimonu.coupling)
    Λ = Sᵀ Yᵀ B                                    kinetic_flux_coupling(crn, route=...) -> FluxCouplingResult
Values of the relations (beyond the main text, see the Supplementary Text; kimonu.relation_values)
    exact values at given k                        relation_values(crn, k)
    k for given concentrations and fluxes          rate_constants(crn, r, x)
    monomial formulas in k from flux coupling      binomial_formulas(crn)

The three routes return one result type, :class:`RouteResult`, with the modules, ACR and ACRR, the ACRR groups, the
basis B they were read from, the flux r* and the route-specific extras.
"""
from kimonu.network import CRN, as_crn, build_crn, crn_from_matrices, make_irreversible, preprocess
from kimonu.concordance import ConcordanceResult, concordance
from kimonu.readout import (Readout, RouteResult, RowReadout, composition_differences, free_directions, readout,
                            row_readout)
from kimonu.steady_state import interior_flux, steady_state_map
from kimonu.iterated import iterated_route
from kimonu.tangent import tangent_route, tangent_space
from kimonu.certify import Certificate, certify
from kimonu.reference import reference_route
from kimonu.coupling import (FluxCouplingResult, MassBalanceCouplingResult, kinetic_flux_coupling,
                             mass_balance_flux_coupling)
from kimonu.relation_values import (Formulas, Monomial, Value, Values, binomial_formulas, rate_constants,
                                    relation_values)

__version__ = "0.0.1"
__all__ = [
    # networks and preprocessing
    "make_irreversible", "build_crn", "preprocess", "crn_from_matrices", "as_crn", "CRN",
    # flux balance
    "mass_balance_flux_coupling", "MassBalanceCouplingResult", "concordance", "ConcordanceResult",
    "row_readout", "RowReadout",
    # mass action, the basis B and its readout
    "readout", "Readout", "RouteResult", "composition_differences", "free_directions",
    # reference steady state
    "interior_flux", "steady_state_map",
    # iterated route
    "iterated_route",
    # tangent route
    "tangent_route", "tangent_space",
    # both routes
    "certify", "Certificate",
    # reference-equivalent mode
    "reference_route",
    # flux coupling under mass action
    "kinetic_flux_coupling", "FluxCouplingResult",
    # values of the relations
    "relation_values", "rate_constants", "binomial_formulas", "Values", "Value", "Formulas", "Monomial",
    "__version__",
]
