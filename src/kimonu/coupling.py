"""Flux coupling: fully coupled reactions under mass balance, and kinetic modules lifted into flux space.

Mass balance (:func:`mass_balance_flux_coupling`). Every steady-state flux is a combination r = K·α of the columns
of a basis K of ker(N). When rows i and j of K are proportional, row j = λ · row i, reactions i and j are fully
coupled, and r_j = λ r_i at every steady state (Burgard et al., 2004). A zero row is a blocked reaction, which no
steady state runs. K is the basis that concordance reads (:func:`kimonu.concordance`), with the pinned rates as rows
of N, so the groups hold over the feasible steady states.

Mass action (:func:`kinetic_flux_coupling`). Under flux balance, two reactions are fully coupled when their rates
keep a fixed ratio over every steady-state flux. Under mass action with fixed rate constants only some of these
fluxes are steady states. With the rate constants fixed by an interior flux r* at x* = 1, as for both routes
(:mod:`kimonu.steady_state`), reaction j runs at

    r_j = r*_j u_c(j),

where c(j) is its substrate complex, u_c = ψ_c / ψ*_c, or, for a detached uptake, a source of its own
(:func:`kimonu.network.uptake_columns`). Two rates keep a fixed ratio exactly when the u of their columns
do. This happens for

- two reactions that consume the same complex, whose ratio is that of their rate constants,
- two reactions that consume complexes of one kinetic module,
- with detached uptakes, an uptake whose source keeps a fixed ratio to a complex or to another source.

A rate is fixed when the u of its column is. That holds for an uptake from an attached zero complex, and for every
reaction that consumes a complex whose concentration is fixed, as that of a complex made of ACR metabolites is.

Both routes read these groups off their basis B, in the way concordance reads M = A K. A reaction runs at its rate
constant times the concentration of its substrate complex, so ln r − ln r* = Sᵀ Yᵀ (ln x − ln x*). A detached
source enters as one more species with a free concentration (Ŷ, the composition of the complexes and the
sources, :func:`kimonu.network.with_sources`), and

    Λ = Sᵀ Ŷᵀ B

lifts B into flux space. Reactions with equal rows of Λ keep a fixed flux ratio, and reactions with a zero row have
a fixed rate. Since ker(W Yᵀ) ⊆ B ⊆ ker((YΔ)ᵀ) (Eq. 4), the groups inherit iterated ⊆ true ⊆ tangent.

Iterated route. B = ker((ŶΔ)ᵀ), with Δ built from the final kinetic modules, the ACR metabolites, and each
source coupled to a complex. That coupling is linear: columns with proportional rows in ker(W), with the
rows u_c = u_d of the final modules added, keep a fixed ratio at every positive steady state. On the
complexes these groups are the kinetic modules, so a source is the only thing they add. Every coupling and
fixed rate found holds at every positive steady state.

Tangent route. B = ker(W Yᵀ), whose rows already include the sources. No true coupling or fixed rate is missed,
and false couplings are possible, as for the kinetic modules.

Every full coupling under flux balance holds under mass action, since diag(r*) Sᵀ maps ker(W) into
ker(N). Within a group, the ratio of two rates at every steady state is their ratio at r*. Directional and
partial coupling have no counterpart, since every rate is positive at a positive steady state.
"""
from collections import namedtuple

import numpy as np

from kimonu.concordance import _with_pins
from kimonu.iterated import coupled_columns, iterated_route
from kimonu.linalg import matrix_product, nullspace
from kimonu.network import as_crn, uptake_columns, with_sources, zero_complexes
from kimonu.readout import free_directions, row_readout
from kimonu.steady_state import interior_flux
from kimonu.tangent import tangent_space

MassBalanceCouplingResult = namedtuple(
    "MassBalanceCouplingResult",
    ["groups", "group_of", "n_groups", "largest_group", "ratio", "blocked", "K"],
)
MassBalanceCouplingResult.__doc__ = """Fully coupled reactions under mass balance, from the proportional rows of K.

groups         sets of fully coupled reactions, partitioning the reactions that are not blocked, singletons
               included, in the order of their first reaction
group_of       numpy array, the group of each reaction, counted from 1, and 0 for a blocked reaction
n_groups       number of groups with more than one reaction
largest_group  size of the largest group
ratio          {reaction: λ} for every reaction that is not blocked, with r_j = λ r_i at every steady state for the
               first reaction i of its group (row j of K = λ · row i), so λ = 1 for the first reaction
blocked        reactions with a zero row of K, which no steady state runs, none after :func:`kimonu.preprocess`
K              the RREF basis of ker(N) in row form, one row per reaction, as in :class:`kimonu.ConcordanceResult`
"""

FluxCouplingResult = namedtuple(
    "FluxCouplingResult",
    ["groups", "group_of", "n_groups", "largest_group", "fixed_rate", "route", "flux"],
)
FluxCouplingResult.__doc__ = """Reactions whose rates keep a fixed ratio under mass action, and the rates that are fixed.

groups         sets of reactions with a fixed ratio of rates, partitioning the reactions, singletons included, in the
               order of their first reaction
group_of       numpy array, the group of each reaction, counted from 1
n_groups       number of groups with more than one reaction
largest_group  size of the largest group
fixed_rate     reactions whose rate is the same at every steady state, proven by the iterated route and reported by
               the tangent route, which can include false couplings
route          "iterated" or "tangent"
flux           the interior flux r* the route ran at. For the rate constants k = r*, two coupled rates keep the ratio
               r*_i / r*_j at every steady state, and a fixed rate keeps the value r*_j.
"""


def mass_balance_flux_coupling(network, preprocess=True, solver=None):
    """Fully coupled reactions under mass balance, with the fixed ratio of their rates, and blocked reactions.

    Args:
        network: a :class:`kimonu.CRN`, or a cobra model, which is preprocessed first (:func:`kimonu.as_crn`). Pass
            the CRN when calling several functions, so that preprocessing runs once.
        preprocess, solver: for a cobra model, as for :func:`kimonu.concordance`.

    The pins of the CRN enter N as rows that fix the ratio of every two pinned rates, as for concordance.

    Returns a :class:`MassBalanceCouplingResult`.
    """
    crn = as_crn(network, preprocess, solver)
    # K = ker(N): every steady-state flux is K·α. Proportional rows are fully coupled reactions, zero rows blocked.
    K = nullspace(_with_pins(crn))
    blocked, groups, ratio = row_readout(K)
    groups = [set(g) for g in groups]
    group_of = np.zeros(len(K), dtype=np.int64)
    for idx, g in enumerate(groups, start=1):
        for j in g:
            group_of[j] = idx
    return MassBalanceCouplingResult(
        groups=groups,
        group_of=group_of,
        n_groups=sum(1 for g in groups if len(g) > 1),
        largest_group=max((len(g) for g in groups), default=0),
        ratio=ratio,
        blocked=blocked,
        K=K,
    )


def _labels(extended, B):
    """Group label of every column of W (complexes, then sources) from the equal rows of Ŷᵀ B, and the label of the
    zero row, the fixed columns, if there is one."""
    keys = {}
    label = {c: keys.setdefault(row, len(keys)) for c, row in enumerate(matrix_product(extended.Yt, B))}   # Ŷᵀ B
    return label, {keys[()]} if () in keys else set()


def _iterated(crn, detach_uptakes, flux):
    """Group label of every column and the labels of the fixed groups."""
    found = iterated_route(crn, detach_uptakes, flux=flux)
    equal = [(min(g), c) for g in found.modules for c in g if c != min(g)]
    groups = coupled_columns(crn, flux, detach_uptakes, equal)
    skip = zero_complexes(crn) if detach_uptakes else set()
    extended = with_sources(crn, detach_uptakes)
    index = {name: i for i, name in enumerate(crn.metabolite_ids)}
    B = free_directions(extended, [set(g) - skip for g in groups], [index[a] for a in found.acr])
    label, fixed = _labels(extended, B)
    for c in skip:                                        # no reaction consumes a detached zero complex
        label[c] = -1 - c
    return label, fixed


def _tangent(crn, detach_uptakes, flux):
    """Group label of every column and the labels of the fixed groups."""
    extended = with_sources(crn, detach_uptakes)
    label, fixed = _labels(extended, tangent_space(crn, flux, detach_uptakes))
    for c in (zero_complexes(crn) if detach_uptakes else ()):
        label[c] = -1 - c
    return label, fixed


def kinetic_flux_coupling(network, route="iterated", detach_uptakes=True, seed=0, draw="lp", bits=64, flux=None):
    """Reactions whose rates keep a fixed ratio under mass action, and the rates that are fixed.

    Args:
        network: a :class:`kimonu.CRN` without blocked reactions, or a cobra model, which is preprocessed first
            (:func:`kimonu.as_crn`). Pass the CRN when calling several functions, so that preprocessing runs once.
        route: "iterated" (default) for couplings proven at every positive steady state, or "tangent" for every
            true coupling and possibly false couplings (see :func:`kimonu.tangent_route`).
        detach_uptakes: give every reaction out of the zero complex its own source (default), as for
            :func:`kimonu.iterated_route`. False attaches the uptakes to the zero complex, so they run at
            constant rate.
        seed, draw, bits: of the reference flux r* (:func:`kimonu.interior_flux`). The results hold for generic rate
            constants.
        flux: a reference flux r* to use instead of drawing one.

    Returns a :class:`FluxCouplingResult`.
    """
    if route not in ("iterated", "tangent"):
        raise ValueError(f"route must be 'iterated' or 'tangent', not {route!r}")
    crn = as_crn(network)
    if flux is None:
        flux = interior_flux(crn, seed, draw, bits)
    label, fixed = (_iterated if route == "iterated" else _tangent)(crn, detach_uptakes, flux)
    column, _ = uptake_columns(crn, detach_uptakes)
    q = crn.A.shape[1]
    by_label = {}
    for j in range(q):
        by_label.setdefault(label[column[j]], []).append(j)
    groups = [set(g) for g in by_label.values()]
    group_of = np.zeros(q, dtype=np.int64)
    for idx, g in enumerate(groups, start=1):
        for j in g:
            group_of[j] = idx
    return FluxCouplingResult(
        groups=groups,
        group_of=group_of,
        n_groups=sum(1 for g in groups if len(g) > 1),
        largest_group=max((len(g) for g in groups), default=0),
        fixed_rate=sorted(j for j in range(q) if label[column[j]] in fixed),
        route=route,
        flux=flux,
    )
