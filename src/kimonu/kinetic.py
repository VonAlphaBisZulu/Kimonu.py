"""Kinetic modules by the iterated route: a null space in complex space, closed in log space.

Under mass action, reaction j runs at r_j = k_j ψ_s(j), where ψ_c = exp(y_cᵀ ln x) is the monomial of
complex c and s(j) the substrate complex of j. Any positive steady state (r*, ψ*) fixes rate constants
k_j = r*_j / ψ*_s(j), and every steady state for those constants satisfies

    W u = 0,   W = N diag(r*) Sᵀ,   u = ψ / ψ*,

with S the substrate part of the incidence matrix. The equations are linear in u; all the nonlinearity
sits in u = exp(Yᵀ(ln x - ln x*)). Every rate-constant set that admits a positive steady state arises
this way, and r* is drawn at random from the interior of the steady-state cone, so what follows holds
for generic rate constants.

Two kinds of relation are combined until neither yields anything new, following the scheme of
Langary et al. (2025, Supplement S.4.1):

1. Linear: complexes whose rows in a basis of ker(W) are proportional keep a fixed ψ ratio at every
   positive steady state; they are kinetically coupled.
2. Log: within a coupled set, ln ψ_c - ln ψ_d = (y_c - y_d)ᵀ ln x is constant, so the composition
   differences span a space of fixed log-combinations. The readout of acr.py groups the complexes
   whose compositions differ by such a combination (Proposition S3-4) and finds the ACR metabolites.
   One readout suffices: the differences within a merged group and the unit vectors of the ACR
   metabolites (Remark S3-6) already lie in the span, so repeating it returns the same readout.
3. Feedback: a merge found in log space is a new linear fact, u_c = u_d for the merged complexes,
   since their ratio is fixed and equals one at u = 1. Added to W, it can make further rows of ker(W)
   proportional, which the proportional-row readout alone misses because it ignores linear relations
   among three or more coupled sets; it plays the part of the concordance merge of Proposition S4-1.
   On EnvZ-OmpR this is what couples the complexes that differ by Yp.

Every relation found holds at every positive steady state, but relations may be missed: the route is
sound, not complete. The tangent route (tangent.py) reads the same quantities off the tangent space of
the steady states, misses nothing locally, and so bounds this route from above.

Complexes that are no reaction's substrate have a zero column in W, so their u is free in ker(W) and
they reach the log step as singletons. Their ψ is still a monomial of x, and the readout joins them to a
module exactly when their composition differs from a member's by a fixed log-combination; this is the
coupling of Langary et al. (2025, Supplement S.1), which is not restricted to substrates. Joining adds no
new vector to the span, so it changes no ACR or ACRR, and the u_c = u_d it feeds back only fixes u_c.

The zero complex has ψ = 1. Detached (default, as in the published tools), each uptake reaction is given
its own source, so only couplings that hold whatever the uptake rates remain, and the zero complex is a
module of its own. Attached, its uptake reactions run at constant rate and complexes coupled to it share its
module: a fixed medium.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np

from kimonu.acr import free_directions, readout
from kimonu.io import substrate_complexes, zero_complexes
from kimonu.linalg import entries, exact, null_space_rows, positive_flux, proportional_classes

KineticResult = namedtuple(
    "KineticResult",
    ["modules", "module_of", "n_modules", "giant_module", "known_acr", "rounds", "flux"],
)


def rate_columns(A, Y, detach_uptakes):
    """The column of u that carries each reaction's rate, {reaction: column}, and the number of columns.

    Reaction j runs at r*_j u_c for its column c: its substrate complex, or, for a detached uptake, a
    source of its own after the n complexes, numbered in the order of the reactions.
    """
    n = A.shape[0]
    zero = zero_complexes(Y) if detach_uptakes else set()
    column, n_columns = {}, n
    for j, c in substrate_complexes(A).items():
        if c in zero:
            column[j], n_columns = n_columns, n_columns + 1
        else:
            column[j] = c
    return column, n_columns


def steady_state_map(A, Y, N, flux, detach_uptakes):
    """W = N diag(r*) Sᵀ as {(metabolite, column): value}, and its shape: dx/dt = W u.

    Column c collects, for every reaction consuming complex c, its stoichiometry times its rate at r*.
    A detached uptake gets a column of its own after the n complexes (:func:`rate_columns`), a free source
    rather than the constant ψ_0 = 1. tangent.py builds its Jacobian from the same map.
    """
    column, n_columns = rate_columns(A, Y, detach_uptakes)
    by_reaction = {}
    triples, (m, _) = entries(N)
    for i, j, v in triples:
        by_reaction.setdefault(j, []).append((i, v))
    W = {}
    for j, col in column.items():
        for i, v in by_reaction.get(j, ()):
            W[(i, col)] = W.get((i, col), Fraction(0)) + flux[j] * v
    return W, (m, n_columns)


def coupled_columns(A, Y, N, flux, detach_uptakes, equal=()):
    """Classes of the columns of W (complexes, then detached sources) whose u keep a fixed ratio.

    ``equal`` lists pairs (c, d) already known to satisfy u_c = u_d at every steady state; they enter W
    as extra rows. Every steady state solves W u = 0, so columns with proportional rows in ker(W) keep a
    fixed ratio at every steady state, and since u = 1 is one of them, that ratio is 1. No row is zero,
    so every column lands in exactly one class.
    """
    W, (m, n_columns) = steady_state_map(A, Y, N, flux, detach_uptakes)
    for k, (c, d) in enumerate(equal):                    # u_c - u_d = 0
        W[(m + k, c)], W[(m + k, d)] = Fraction(1), Fraction(-1)
    W = exact(((i, j, v) for (i, j), v in W.items()), (m + len(equal), n_columns))
    _, classes = proportional_classes(null_space_rows(W))
    return classes


def _psi_coupled(A, Y, N, flux, detach_uptakes, equal=()):
    """Complexes whose ψ ratios are fixed by the steady-state equations W u = 0.

    ``equal`` lists pairs (c, d) already known to satisfy u_c = u_d. Returns the coupled sets and the
    substrate complexes (the only complexes that appear in a rate law).
    """
    n = A.shape[0]
    classes = coupled_columns(A, Y, N, flux, detach_uptakes, equal)
    drop = zero_complexes(Y) if detach_uptakes else set()
    return [{c for c in g if c < n and c not in drop} for g in classes], set(substrate_complexes(A).values()) - drop


def _closure(Y, sets, skip):
    """Merge sets whose compositions differ by a fixed log-combination (Proposition S3-4) and find the
    ACR metabolites, in one readout at L = ker(YΔᵀ).

    Within a set every composition difference is a column of YΔ, so all its members have the same row
    of Yᵀ L and the readout only merges whole sets. One readout is already closed: complexes it groups
    differ by a vector of im(YΔ), and the unit vector of an ACR metabolite lies in im(YΔ), so YΔ built
    from the merged sets and those unit vectors (Remark S3-6) has the same span, the same L and the same
    readout.
    """
    found = readout(Y, free_directions(Y, sets), skip)
    return found.modules, set(found.acr)


def kinetic_modules(A, Y, N, detach_uptakes=True, seed=0):
    """Kinetic modules of a network with no blocked reactions, by the iterated route.

    Args:
        A: complexes x reactions incidence matrix, one -1 and one +1 per column.
        Y: metabolites x complexes composition matrix.
        N: metabolites x reactions stoichiometric matrix, exact (straindesign RationalMatrix) or
            integer; N = Y A.
        detach_uptakes: give every reaction out of the zero complex its own source (default), so
            couplings through the constant zero-complex monomial are not used; False attaches the
            uptakes to the zero complex, a fixed medium.
        seed: seed of the random interior flux. Results hold for generic rate constants, so the
            seed changes them only with probability zero.

    Returns a :class:`KineticResult`; every complex belongs to exactly one module, singletons
    included, and modules are ordered by their first complex. ``known_acr`` are the metabolites the
    route proves to have ACR; :func:`kimonu.acr_acrr` reports them again together with the ACRR pairs.
    ``flux`` is the interior flux r* the route ran at: the rate constants k = r* have the steady state
    x = 1, and :func:`kimonu.rate_constants` moves it to any other x.
    """
    n = A.shape[0]
    skip = zero_complexes(Y) if detach_uptakes else set()
    # 1) Generic rate constants: an interior steady-state flux r* fixes k_j = r*_j at ψ* = 1.
    flux = positive_flux(N, seed)
    equal, rounds = [], 0
    while True:
        rounds += 1
        # 2) Linear part: every steady state for these k solves W u = 0; proportional rows of ker(W)
        #    are complexes with fixed ψ ratios.
        coupled, _ = _psi_coupled(A, Y, N, flux, detach_uptakes, equal)
        coupled = [s for s in coupled if s]
        # 3) Log part: merge sets whose compositions differ by a fixed log-combination, and fix the
        #    metabolites that become ACR on the way.
        sets, known = _closure(Y, coupled, skip)
        # 4) Feedback: a merge found in log space is the linear fact u_c = u_d; add it to W and repeat
        #    until nothing merges.
        if len(sets) == len(coupled):
            break
        equal = [(min(s), c) for s in sets for c in s if c != min(s)]

    # 5) A detached zero complex is coupled to nothing and forms a module of its own.
    modules = sets + [{c} for c in sorted(skip)]
    module_of = np.zeros(n, dtype=np.int64)
    for idx, group in enumerate(modules, start=1):
        for i in group:
            module_of[i] = idx

    return KineticResult(
        modules=modules,
        module_of=module_of,
        n_modules=sum(1 for m in modules if len(m) > 1),
        giant_module=max((len(m) for m in modules), default=0),
        known_acr=sorted(known),
        rounds=rounds,
        flux=flux,
    )
