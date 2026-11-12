"""The reference steady state: the exact interior flux r* of both draws, and W = N diag(r*) Sᵀ of Eq. 3.

W is checked by hand on the open network, and against its defining property W·1 = N r* = 0 on EnvZ-OmpR.
"""
from fractions import Fraction

import pytest

import kimonu
from kimonu.network import uptake_columns, with_sources
from networks import ENVZ, ENVZ_SPECIES, OPEN, network


def _steady(crn, flux):
    Nd = crn.N.to_numpy()
    return all(sum(Fraction(int(Nd[i, j])) * flux[j] for j in range(Nd.shape[1])) == 0 for i in range(Nd.shape[0]))


@pytest.mark.parametrize("draw", ["uniform", "lp"])
def test_interior_flux_is_an_exact_integer_interior_steady_state(draw):
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    flux = kimonu.interior_flux(crn, draw=draw)
    assert all(isinstance(x, Fraction) and x > 0 and x.denominator == 1 for x in flux)
    assert _steady(crn, flux)


def test_the_uniform_draw_reaches_beyond_the_lp_point():
    """The uniform draw adds integers up to 2^64 to the free coordinates, so its entries are far larger, and it
    depends on the seed."""
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    lp, uniform = kimonu.interior_flux(crn, draw="lp"), kimonu.interior_flux(crn, draw="uniform")
    assert max(uniform) > 2**64 > max(lp)
    assert kimonu.interior_flux(crn, seed=1, draw="uniform") != uniform
    assert kimonu.interior_flux(crn, draw="uniform") == uniform
    assert kimonu.interior_flux(crn) == lp


def test_bits_set_the_size_of_the_uniform_draw():
    """r* = K (T c + s) with s below 2^bits and T = ceil(2 · 2^bits · reach / min r_c), so r* has at most about bits
    more bits than the LP point r_c (T is 1 while 2^bits is small against min r_c). bits=64 is the default of the
    uniform draw, and the LP point (the default draw) ignores bits."""
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    lp_flux = kimonu.interior_flux(crn, draw="lp")
    lp = max(lp_flux).numerator.bit_length()
    for bits in (8, 20, 64):
        flux = kimonu.interior_flux(crn, draw="uniform", bits=bits)
        assert all(x > 0 and x.denominator == 1 for x in flux) and _steady(crn, flux) and flux != lp_flux
        assert max(flux).numerator.bit_length() <= lp + bits + 12
    assert max(kimonu.interior_flux(crn, draw="uniform", bits=64)).numerator.bit_length() > lp + 40
    assert kimonu.interior_flux(crn, draw="uniform", bits=64) == kimonu.interior_flux(crn, draw="uniform")
    assert kimonu.interior_flux(crn, draw="lp", bits=8) == kimonu.interior_flux(crn, draw="lp")
    with pytest.raises(ValueError):
        kimonu.interior_flux(crn, bits=0)


def test_bits_reach_every_route_and_certify():
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    flux = kimonu.interior_flux(crn, bits=16)
    assert kimonu.iterated_route(crn, bits=16).flux == flux == kimonu.tangent_route(crn, bits=16).flux
    cert = kimonu.certify(crn, bits=16)
    assert cert.flux == flux and cert.certified


def test_unknown_draw():
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    with pytest.raises(ValueError):
        kimonu.interior_flux(crn, draw="sampled")


def test_blocked_network_is_rejected():
    """a -> b with nothing consuming b has no positive steady state."""
    crn, cx = network([("a", "b")], ["a", "b"])
    with pytest.raises(ValueError):
        kimonu.interior_flux(crn)


def test_uptake_columns_give_a_detached_uptake_its_own_source():
    """0 -> a, a <-> b, b -> 0. Attached, the uptake runs from the zero complex. Detached, from a source after the
    three complexes."""
    crn, cx = network(OPEN, ["a", "b"])
    zero, a, b = cx.index("0"), cx.index("a"), cx.index("b")
    assert uptake_columns(crn, detach_uptakes=False) == ({0: zero, 1: a, 2: b, 3: b}, 3)
    assert uptake_columns(crn) == ({0: 3, 1: a, 2: b, 3: b}, 4)


def test_with_sources_adds_one_species_per_detached_uptake():
    """Detached, the uptake of a draws from a source, a species of its own after a and b with composition e_2 and a
    row of its own after the three complexes. Attached, the network stays as it is."""
    crn, cx = network(OPEN, ["a", "b"])
    species = with_sources(crn)
    assert species.metabolite_ids == ["a", "b", "source of 0"] and species.Yt == crn.Yt + [((2, 1),)]
    assert with_sources(crn, detach_uptakes=False) is crn


def test_steady_state_map_by_hand():
    """Column c of W sums the columns of N of the reactions that consume c, weighted by r*."""
    crn, cx = network(OPEN, ["a", "b"])
    r = [Fraction(v) for v in (2, 5, 3, 2)]
    W, n_columns = kimonu.steady_state_map(crn, r, detach_uptakes=False)
    zero, a, b = cx.index("0"), cx.index("a"), cx.index("b")
    assert n_columns == 3 and len(W) == 2
    assert W == [tuple(sorted([(zero, 2), (a, -5), (b, 3)])), tuple(sorted([(a, 5), (b, -5)]))]
    W, n_columns = kimonu.steady_state_map(crn, r)
    assert n_columns == 4 and (3, 2) in W[0] and zero not in dict(W[0])


def test_w_annihilates_the_reference():
    """u = 1 is the reference steady state, so W·1 = N r* = 0."""
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    W, n_columns = kimonu.steady_state_map(crn, kimonu.interior_flux(crn))
    assert len(W) == len(ENVZ_SPECIES) and all(sum(v for _, v in row) == 0 for row in W)
