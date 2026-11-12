"""The tangent route: one nullspace at one generic steady state, which misses no relation but can add false ones.

The iterated route (:mod:`kimonu.iterated`) builds the fixed combinations up from relations it can prove. The
tangent route applies both conditions of Eq. 1 at once, but only to first order around the reference steady state,
as a tangent line approximates a curve near one point. Pascual-Escudero & Feliu (2022) introduced this
linearization as a test for local ACR. Kimonu extends it to ACRR and kinetic modules and evaluates it as one exact
sparse nullspace at a random interior point, which runs at genome scale.

With the rate constants fixed by an interior flux r* at x* = 1, as for the iterated route, let ln x move away from
ln x* by a small step h. Then ln ψ moves by Yᵀ h, so u = exp(Yᵀ h) ≈ 1 + Yᵀ h. Since W·1 = 0 at the reference,
staying at steady state requires, to first order,

    W Yᵀ h = 0,    W = N diag(r*) Sᵀ (Eq. 3),

so ker(W Yᵀ) holds the directions in which ln x can move without leaving the steady states, to first order. It is
read out like B (:func:`kimonu.readout.readout`), without iteration. Equal rows of Yᵀ ker(W Yᵀ) are kinetic
modules, zero rows of ker(W Yᵀ) are ACR metabolites and equal nonzero rows are ACRR groups.

Why no relation is missed. A steady state z near z* = 0 (z = ln x) is itself the point z = 0 for the rate constants
fixed by its own flux r* ∘ exp(Sᵀ Yᵀ z), so the tangent spaces at nearby steady states are the ker(W Yᵀ) of nearby
interior fluxes. For generic r* the rank of W Yᵀ does not change between nearby interior fluxes, and the
constant-rank theorem then makes the steady states near x* a manifold with tangent space ker(W Yᵀ). A combination
that is fixed at every steady state is constant on this manifold, so it is orthogonal to ker(W Yᵀ), and the route
reports it. ker(W Yᵀ) is therefore contained in B (Eq. 4, right). The argument does not need the steady state to
be nondegenerate.

False couplings. The route checks only whether a combination changes when the steady state moves a little. Steady
states that no continuous path of steady states connects to the reference are never compared, and a combination can
take different values there. The route then reports a relation that does not hold, a false coupling. The
self-activating kinase of Fig. 2, 0 <-> X, 2X -> C, C + 2P -> 3P, C -> 0, P -> 0, has the two steady states X = 10,
(C, P) = (1, 4) and (4, 1) for one choice of rate constants. Neither can move without leaving the steady states, so
ker(W Yᵀ) = {0} at both, and the route reports ACR for X, C and P, C and P falsely. The iterated route claims
neither C nor P and misses X. Where the two routes agree, the list is correct and complete (:func:`kimonu.certify`).

The zero complex has ψ_∅ = 1. Detached (default), each uptake reaction draws from a source of its own, a species
with a free concentration (:func:`kimonu.steady_state.with_sources`), so W Yᵀ becomes W Ŷᵀ with one more coordinate
per source. Only the metabolite rows are read, and the zero complex is a module of its own, listed last. Attached,
its uptake reactions run at constant rate and add nothing to W Yᵀ.
"""
from kimonu.linalg import matrix_product, null_space_rows
from kimonu.readout import readout, route_result
from kimonu.steady_state import interior_flux, steady_state_map, with_sources


def tangent_space(crn, flux, detach_uptakes=True):
    """An exact basis of ker(W Yᵀ) at the steady state x = 1 for the rate constants k = ``flux``, in row form.

    One row per coordinate. First come the m metabolites (ln x), then, with detached uptakes, the sources (ln u) in
    the order of :func:`kimonu.steady_state.rate_columns`, the coordinates of :func:`kimonu.steady_state.with_sources`.
    The first m rows are the B of the tangent route.
    """
    W, _ = steady_state_map(crn, flux, detach_uptakes)
    species = with_sources(crn, detach_uptakes)
    # W Ŷᵀ: the column of complex c spreads over c's metabolites, and a detached source keeps a column of its own
    return null_space_rows(matrix_product(W, species.Yt), len(species.metabolite_ids))


def tangent_route(crn, detach_uptakes=True, seed=0, draw="lp", bits=64, flux=None):
    """Kinetic modules, ACR and ACRR by the tangent route, which reports every true relation and possibly false ones.

    Every relation that holds at every positive steady state (for generic rate constants) is reported. A reported
    relation can be a false coupling, one that differs between steady states that no path of steady states
    connects (see the module docstring). :func:`kimonu.certify` compares with the iterated route at the same r*.

    Args:
        crn: the network, as for :func:`kimonu.iterated_route`, without blocked reactions.
        detach_uptakes: give every reaction out of the zero complex its own source (default), as for
            :func:`kimonu.iterated_route`.
        seed, draw, bits: of the reference flux r* (:func:`kimonu.interior_flux`). The default LP point keeps the
            integers small. draw="uniform" adds the probability bound and makes this route several times slower.
        flux: a reference flux r* to use instead of drawing one.

    Returns a :class:`kimonu.RouteResult`, with ``dimension`` the dimension of ker(W Yᵀ), sources included.
    """
    # 1) Generic rate constants from an interior flux r*, as for the iterated route
    if flux is None:
        flux = interior_flux(crn, seed, draw, bits)
    # 2) ker(W Yᵀ), one row per coordinate. The metabolite rows are B.
    rows = tangent_space(crn, flux, detach_uptakes)
    B = rows[:len(crn.metabolite_ids)]
    # 3) The readout. With detached uptakes the zero complex is a module of its own, listed last.
    return route_result("tangent", crn, readout(crn, B, detach_uptakes), B=B,
                        dimension=len({j for r in rows for j, _ in r}), flux=flux)
