"""Flux coupling under mass action: kinetic modules lifted into flux space.

Under flux balance, two reactions are fully coupled when their rates keep a fixed ratio over every
steady-state flux, which makes their rows of K = ker(N) proportional. Under mass action with fixed rate
constants only some of these fluxes are steady states. With the rate constants fixed by an interior flux
r* at x* = 1, as in kinetic.py, reaction j runs at

    r_j = r*_j u_c(j),

where c(j) is its substrate complex, u_c = ψ_c / ψ*_c, or, for a detached uptake, a source of its own
(:func:`kimonu.kinetic.rate_columns`). Two rates keep a fixed ratio exactly when the u of their columns
do. This happens for

- two reactions that leave the same complex, whose ratio is that of their rate constants,
- two reactions that leave complexes of one kinetic module,
- with detached uptakes, an uptake whose source keeps a fixed ratio to a complex or to another source.

A rate is fixed when the u of its column is: an uptake from an attached zero complex, and every
reaction that leaves a complex whose monomial is fixed, as a complex made of ACR metabolites is.

Both routes read these classes off their basis B of free log-directions, in the way concordance reads
M = A K. A detached source enters as one more species with a free concentration (Ŷ, the composition of the
complexes and the sources), so that ln r_j − ln r*_j = ŷ_c(j)ᵀ ln x̂ for every reaction, and

    Λ = Sᵀ Ŷᵀ B.

Reactions with equal rows of Λ keep a fixed flux ratio, and reactions with a zero row have a fixed rate.
Since ker(W Yᵀ) ⊆ B ⊆ ker((YΔ)ᵀ), the classes inherit iterated ⊆ true ⊆ tangent.

Iterated route. B = ker((ŶΔ)ᵀ), with Δ built from the final kinetic modules, the ACR metabolites, and each
source coupled to a complex. That coupling is linear: columns with proportional rows in ker(W), with the
rows u_c = u_d of the final modules added, keep a fixed ratio at every positive steady state. On the
complexes these classes are the kinetic modules, so a source is the only thing they add. Every coupling and
fixed rate found holds at every positive steady state.

Tangent route. B = T = ker(W Yᵀ), whose rows already include the sources. Couplings and fixed rates hold
near the steady state, and none is missed there.

Every full coupling under flux balance holds under mass action, since diag(r*) Sᵀ maps ker(W) into
ker(N). Within a class, the ratio of two rates at every steady state is their ratio at r*. Directional and
partial coupling have no counterpart, since every rate is positive at a positive steady state.
"""
from collections import namedtuple

import numpy as np
from scipy import sparse

from kimonu.acr import _combined, free_directions
from kimonu.io import zero_complexes
from kimonu.kinetic import coupled_columns, kinetic_modules, rate_columns
from kimonu.linalg import positive_flux
from kimonu.tangent import tangent_space

FluxCouplingResult = namedtuple(
    "FluxCouplingResult",
    ["classes", "class_of", "n_classes", "largest_class", "fixed_rate", "route", "flux"],
)


def _with_sources(Y, n_columns):
    """Ŷ: the composition of the complexes, and one species of its own for each detached source."""
    return sparse.block_diag((sparse.csc_matrix(Y), sparse.identity(n_columns - Y.shape[1])), format="csc")


def _read(Yhat, rows, n_columns):
    """Class label of every column from the equal rows of Ŷᵀ B, and the label of the zero row, if any."""
    keys = {}
    label = {c: keys.setdefault(_combined(Yhat, rows, c), len(keys)) for c in range(n_columns)}
    return label, {keys[()]} if () in keys else set()


def _iterated(A, Y, N, detach_uptakes, seed):
    """Class label of every column, the labels of the fixed classes, and r*."""
    n = A.shape[0]
    found = kinetic_modules(A, Y, N, detach_uptakes, seed)
    equal = [(min(g), c) for g in found.modules for c in g if c != min(g)]
    classes = coupled_columns(A, Y, N, found.flux, detach_uptakes, equal)
    n_columns = sum(len(g) for g in classes)
    skip = zero_complexes(Y) if detach_uptakes else set()
    Yhat = _with_sources(Y, n_columns)
    rows = free_directions(Yhat, [set(g) - skip for g in classes], found.known_acr)
    label, fixed = _read(Yhat, rows, n_columns)
    for c in skip:                                        # no reaction leaves a detached zero complex
        label[c] = -1 - c
    return label, fixed, found.flux


def _tangent(A, Y, N, detach_uptakes, seed):
    """Class label of every column, the labels of the locally fixed classes, and r*."""
    flux = positive_flux(N, seed, scale=10**6)            # as in tangent_route
    rows = tangent_space(A, Y, N, flux, detach_uptakes)
    _, n_columns = rate_columns(A, Y, detach_uptakes)
    label, fixed = _read(_with_sources(Y, n_columns), rows, n_columns)
    for c in (zero_complexes(Y) if detach_uptakes else ()):
        label[c] = -1 - c
    return label, fixed, flux


def kinetic_flux_coupling(A, Y, N, detach_uptakes=True, seed=0, route="iterated"):
    """Reactions whose rates keep a fixed ratio under mass action, and the rates that are fixed.

    Args:
        A: complexes x reactions incidence matrix, one -1 and one +1 per column.
        Y: metabolites x complexes composition matrix.
        N: metabolites x reactions stoichiometric matrix, exact (straindesign RationalMatrix) or
            integer; N = Y A. The network must have no blocked reactions.
        detach_uptakes: give every reaction out of the zero complex its own source (default), as in
            :func:`kimonu.kinetic_modules`; False attaches the uptakes to the zero complex, so they run
            at constant rate.
        seed: seed of the random interior flux; results hold for generic rate constants.
        route: "iterated" (default) for couplings that hold at every positive steady state, or
            "tangent" for every coupling that holds near the steady state, which can include couplings
            that differ between separate steady states (see :func:`kimonu.tangent_route`).

    Returns a :class:`FluxCouplingResult`. ``classes`` partition the reactions, singletons included,
    ordered by their first reaction, and ``class_of`` gives each reaction's class, counted from 1.
    ``n_classes`` counts the classes with more than one reaction. ``fixed_rate`` are the reactions whose
    rate is the same at every steady state (iterated) or near the steady state (tangent). ``flux`` is the
    interior flux r* the route ran at: for the rate constants k = r*, two coupled rates keep the ratio
    r*_i / r*_j at every steady state, and a fixed rate keeps the value r*_j.
    """
    if route == "iterated":
        label, fixed, flux = _iterated(A, Y, N, detach_uptakes, seed)
    elif route == "tangent":
        label, fixed, flux = _tangent(A, Y, N, detach_uptakes, seed)
    else:
        raise ValueError(f"route must be 'iterated' or 'tangent', not {route!r}")
    column, _ = rate_columns(A, Y, detach_uptakes)
    groups = {}
    for j in range(A.shape[1]):
        groups.setdefault(label[column[j]], []).append(j)
    classes = [set(g) for g in groups.values()]
    class_of = np.zeros(A.shape[1], dtype=np.int64)
    for idx, g in enumerate(classes, start=1):
        for j in g:
            class_of[j] = idx
    return FluxCouplingResult(
        classes=classes,
        class_of=class_of,
        n_classes=sum(1 for g in classes if len(g) > 1),
        largest_class=max((len(g) for g in classes), default=0),
        fixed_rate=sorted(j for j in range(A.shape[1]) if label[column[j]] in fixed),
        route=route,
        flux=flux,
    )
