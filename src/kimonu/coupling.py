"""Flux coupling under mass action: kinetic modules lifted into flux space.

Under flux balance, two reactions are fully coupled when their rates keep a fixed ratio over every
steady-state flux, which makes their rows of K = ker(N) proportional. Under mass action with fixed rate
constants only some of these fluxes are steady states. With the rate constants fixed by an interior flux
r* at x* = 1, as for both routes (:mod:`kimonu.steady_state`), reaction j runs at

    r_j = r*_j u_c(j),

where c(j) is its substrate complex, u_c = ψ_c / ψ*_c, or, for a detached uptake, a source of its own
(:func:`kimonu.steady_state.rate_columns`). Two rates keep a fixed ratio exactly when the u of their columns
do. This happens for

- two reactions that consume the same complex, whose ratio is that of their rate constants,
- two reactions that consume complexes of one kinetic module,
- with detached uptakes, an uptake whose source keeps a fixed ratio to a complex or to another source.

A rate is fixed when the u of its column is. That holds for an uptake from an attached zero complex, and for every
reaction that consumes a complex whose concentration is fixed, as that of a complex made of ACR metabolites is.

Both routes read these classes off their basis B, in the way concordance reads M = A K. A reaction runs at its rate
constant times the concentration of its substrate complex, so ln r − ln r* = Sᵀ Yᵀ (ln x − ln x*). A detached
source enters as one more metabolite with a free concentration (Ŷ, the composition of the complexes and the
sources), and

    Λ = Sᵀ Ŷᵀ B

lifts B into flux space. Reactions with equal rows of Λ keep a fixed flux ratio, and reactions with a zero row have
a fixed rate. Since ker(W Yᵀ) ⊆ B ⊆ ker((YΔ)ᵀ) (Eq. 4), the classes inherit iterated ⊆ true ⊆ tangent.

Iterated route. B = ker((ŶΔ)ᵀ), with Δ built from the final kinetic modules, the ACR metabolites, and each
source coupled to a complex. That coupling is linear: columns with proportional rows in ker(W), with the
rows u_c = u_d of the final modules added, keep a fixed ratio at every positive steady state. On the
complexes these classes are the kinetic modules, so a source is the only thing they add. Every coupling and
fixed rate found holds at every positive steady state.

Tangent route. B = ker(W Yᵀ), whose rows already include the sources. No true coupling or fixed rate is missed,
and false couplings are possible, as for the kinetic modules.

Every full coupling under flux balance holds under mass action, since diag(r*) Sᵀ maps ker(W) into
ker(N). Within a class, the ratio of two rates at every steady state is their ratio at r*. Directional and
partial coupling have no counterpart, since every rate is positive at a positive steady state.
"""
from collections import namedtuple

import numpy as np
from scipy import sparse

from kimonu.io import zero_complexes
from kimonu.iterated import coupled_columns, iterated_route
from kimonu.readout import _combined, free_directions
from kimonu.steady_state import interior_flux, rate_columns
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


def _iterated(A, Y, N, detach_uptakes, flux):
    """Class label of every column and the labels of the fixed classes."""
    found = iterated_route(A, Y, N, detach_uptakes, flux=flux)
    equal = [(min(g), c) for g in found.modules for c in g if c != min(g)]
    classes = coupled_columns(A, Y, N, flux, detach_uptakes, equal)
    n_columns = sum(len(g) for g in classes)
    skip = zero_complexes(Y) if detach_uptakes else set()
    Yhat = _with_sources(Y, n_columns)
    rows = free_directions(Yhat, [set(g) - skip for g in classes], found.acr)
    label, fixed = _read(Yhat, rows, n_columns)
    for c in skip:                                        # no reaction consumes a detached zero complex
        label[c] = -1 - c
    return label, fixed


def _tangent(A, Y, N, detach_uptakes, flux):
    """Class label of every column and the labels of the fixed classes."""
    rows = tangent_space(A, Y, N, flux, detach_uptakes)
    _, n_columns = rate_columns(A, Y, detach_uptakes)
    label, fixed = _read(_with_sources(Y, n_columns), rows, n_columns)
    for c in (zero_complexes(Y) if detach_uptakes else ()):
        label[c] = -1 - c
    return label, fixed


def kinetic_flux_coupling(A, Y, N, detach_uptakes=True, seed=0, route="iterated", draw="lp", bits=64,
                          flux=None):
    """Reactions whose rates keep a fixed ratio under mass action, and the rates that are fixed.

    Args:
        A: complexes x reactions incidence matrix, one -1 and one +1 per column.
        Y: metabolites x complexes composition matrix.
        N: metabolites x reactions stoichiometric matrix, exact (straindesign RationalMatrix) or
            integer; N = Y A. The network must have no blocked reactions.
        detach_uptakes: give every reaction out of the zero complex its own source (default), as for
            :func:`kimonu.iterated_route`. False attaches the uptakes to the zero complex, so they run at
            constant rate.
        seed, draw, bits: of the reference flux r* (:func:`kimonu.interior_flux`). The results hold for generic rate
            constants.
        route: "iterated" (default) for couplings proven at every positive steady state, or "tangent" for every
            true coupling and possibly false couplings (see :func:`kimonu.tangent_route`).
        flux: a reference flux r* to use instead of drawing one.

    Returns a :class:`FluxCouplingResult`. ``classes`` partition the reactions, singletons included,
    ordered by their first reaction, and ``class_of`` gives each reaction's class, counted from 1.
    ``n_classes`` counts the classes with more than one reaction. ``fixed_rate`` are the reactions whose
    rate is the same at every steady state (proven by the iterated route, reported by the tangent route, which
    can include false couplings). ``flux`` is the
    interior flux r* the route ran at: for the rate constants k = r*, two coupled rates keep the ratio
    r*_i / r*_j at every steady state, and a fixed rate keeps the value r*_j.
    """
    if route not in ("iterated", "tangent"):
        raise ValueError(f"route must be 'iterated' or 'tangent', not {route!r}")
    if flux is None:
        flux = interior_flux(N, seed, draw, bits)
    label, fixed = (_iterated if route == "iterated" else _tangent)(A, Y, N, detach_uptakes, flux)
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
