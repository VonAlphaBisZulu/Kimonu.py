"""The iterated route: kinetic modules, ACR and ACRR proven for generic rate constants, possibly not all of them.

A positive steady state satisfies two conditions, each linear on its own scale (Eq. 1). The mass balance
W u = 0 (Eq. 3, :mod:`kimonu.steady_state`) is linear in the complex concentrations u = ψ / ψ*, measured in units of
a reference steady state with flux r*, and the composition ln ψ = Yᵀ ln x is linear in the logs. The iterated route
alternates between the two, following the scheme of Langary et al. (2025), with exact nullspaces in place of their
linear programs and coupled groups read as equal rows.

1. Mass balance. If rows c and d of a basis of ker(W) are proportional, u_c / u_d is the same in every solution.
   Since u = 1 is a solution, the ratio is 1, so the rows are in fact equal and u_c = u_d. ψ_c / ψ_d is then fixed,
   and c and d are coupled. Groups of such complexes are the coupled sets (:func:`coupled_columns`).
2. Composition. Within a coupled set, ln ψ_c - ln ψ_d = (y_c - y_d)ᵀ ln x is fixed (Eq. 2). These composition
   differences are the columns of YΔ, with Δ (complexes x pairs) holding one column e_c - e_d for each member c of a
   set other than its first member d. Every column of YΔ is a combination already proven fixed, so
   ker((YΔ)ᵀ) contains B (Eq. 4), and its zero rows and equal rows are relations of B
   (:func:`kimonu.readout.free_directions`, :func:`kimonu.readout.readout`). Equal rows of Yᵀ ker((YΔ)ᵀ) can merge
   coupled sets, for instance two complexes that share a metabolite and differ by an ACRR pair.
3. Feedback. The mass balance treats every u_c as an unknown of its own, so a merge found on the log scale is new
   information for it. Every complex c of a merged set adds the row u_c = u_d against the first member d of its new
   set to W. The nullspace shrinks, and the smaller nullspace can couple further complexes. On the EnvZ–OmpR network of
   Fig. 1 this joins the set that contains G to the rest and so gives G its ACR.

The rounds end when the readout merges no sets, and that last readout gives the modules, ACR and ACRR. The test
compares the two partitions of the complexes, the coupled sets and the modules of the readout, which are equal
exactly when nothing merged (see the comment at the test). The loop always ends, since every round that continues
merges at least two sets. One readout per round is enough. Complexes it groups differ by a vector in the span of
YΔ, and the unit vector of an ACR metabolite lies there too, so YΔ built from the merged sets and those unit vectors
has the same span and gives the same readout.

Every relation reported is proven for generic rate constants, but the list can be incomplete, because the
composition enters only for complexes that are already coupled. The self-activating kinase of Fig. 2 has ACR in X,
which this route misses. :func:`kimonu.certify` runs the tangent route at the same r* and tells whether the list is
also complete.

Complexes that are no reaction's substrate have a zero column in W, so their u is free in ker(W) and they reach the
composition step as singletons. Their ψ is still a product of metabolite concentrations, and the readout joins them
to a module exactly when their composition differs from a member's by a fixed combination. This coupling is not
restricted to substrates (Langary et al., 2025). Joining adds no new vector to the span of YΔ, so it changes no ACR
or ACRR, and the row u_c = u_d it feeds back only fixes u_c.

The zero complex has ψ_∅ = 1. Detached (the default, as in the published tools), each uptake reaction draws from a
source of its own, so only couplings that hold for every supply remain. The zero complex then consumes nothing in W,
is a coupled set of its own, and is a module of its own, listed last. Attached, its uptake reactions run at constant
rate, and complexes coupled to it share its module, a fixed medium.
"""
from kimonu.linalg import nullspace
from kimonu.network import as_crn
from kimonu.readout import free_directions, readout, route_result, row_readout
from kimonu.steady_state import interior_flux, steady_state_map


def coupled_columns(crn, flux, detach_uptakes=True, equal=()):
    """Groups of the columns of W (complexes, then detached sources) whose u keep a fixed ratio, from ker(W).

    ``equal`` lists pairs (c, d) already known to satisfy u_c = u_d at every steady state, the feedback rows. They
    enter W as extra rows. Every steady state solves W u = 0, so columns with proportional rows in ker(W) keep a
    fixed ratio at every steady state, and since u = 1 is one of them, that ratio is 1. For the same reason no row of
    ker(W) is zero, so every column lands in exactly one group. Restricted to the complexes, the groups are the
    coupled sets.
    """
    W, n_columns = steady_state_map(crn, flux, detach_uptakes)
    W += [tuple(sorted(((c, 1), (d, -1)))) for c, d in equal]                 # u_c - u_d = 0
    return row_readout(nullspace(W, n_columns)).groups


def iterated_route(network, detach_uptakes=True, seed=0, draw="lp", bits=64, flux=None):
    """Kinetic modules, ACR and ACRR by the iterated route, each proven for generic rate constants.

    The default analysis of Kimonu. :func:`kimonu.certify` adds the tangent route at the same r* and tells whether
    the list is complete.

    Args:
        network: a :class:`kimonu.CRN` without blocked reactions, or a cobra model, which is preprocessed first
            (:func:`kimonu.as_crn`). Pass the CRN when calling several functions, so that preprocessing runs once.
        detach_uptakes: give every reaction out of the zero complex its own source (default), so only relations
            that hold for every supply are reported. False attaches the uptakes to the zero complex, a fixed medium.
        seed, draw, bits: of the reference flux r* (:func:`kimonu.interior_flux`). The results do not depend on
            them, except for finely tuned reference fluxes, which a random draw does not hit.
        flux: a reference flux r* to use instead of drawing one, as :func:`kimonu.certify` does.

    Returns a :class:`kimonu.RouteResult` with ``rounds``, the number of rounds of mass balance and composition, and
    ``dimension``, dim ker((YΔ)ᵀ) of the last round.
    """
    crn = as_crn(network)
    n = len(crn.Yt)
    # 1) Generic rate constants: an interior steady-state flux r* fixes k = r* with the steady state x* = 1.
    if flux is None:
        flux = interior_flux(crn, seed, draw, bits)
    equal, rounds = [], 0
    while True:
        rounds += 1
        # 2) Mass balance: equal rows of ker(W), with the feedback rows u_c = u_d, are the coupled sets.
        groups = coupled_columns(crn, flux, detach_uptakes, equal)
        coupled = [s for s in ({c for c in g if c < n} for g in groups) if s]
        # 3) Composition: B = ker((YΔ)ᵀ) of the coupled sets. Its readout merges sets whose compositions differ by
        #    a fixed combination and finds ACR and ACRR.
        B = free_directions(crn, coupled)
        found = readout(crn, B, detach_uptakes)
        # 4) Stop when the readout merged no sets. Both lists partition the complexes. Every complex is in exactly
        #    one coupled set, since no row of ker(W) is zero, and in exactly one module. Each module is a union of
        #    coupled sets, since the composition differences within a set are columns of YΔ and so give equal rows
        #    of Yᵀ B. The partitions are therefore equal exactly when nothing merged. Otherwise every merge is the
        #    linear fact u_c = u_d, a feedback row of W, and the next round starts with fewer sets.
        if {frozenset(s) for s in found.modules} == {frozenset(s) for s in coupled}:
            break
        equal = [(min(s), c) for s in found.modules for c in s if c != min(s)]
    return route_result("iterated", crn, found, B=B, dimension=len({j for row in B for j, _ in row}), flux=flux,
                        rounds=rounds)
