"""The reference steady state: the exact interior flux r* of both draws, and W = N diag(r*) Sᵀ of Eq. 3.

W is checked by hand on the open network, and against its defining property W·1 = N r* = 0 on EnvZ-OmpR.
"""
from fractions import Fraction

import pytest

import kimonu
from kimonu.steady_state import rate_columns
from networks import ENVZ, ENVZ_SPECIES, OPEN, matrices


def _steady(N, flux):
    Nd = N.toarray()
    return all(sum(Fraction(int(Nd[i, j])) * flux[j] for j in range(Nd.shape[1])) == 0 for i in range(Nd.shape[0]))


@pytest.mark.parametrize("draw", ["uniform", "lp"])
def test_interior_flux_is_an_exact_integer_interior_steady_state(draw):
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    flux = kimonu.interior_flux(N, draw=draw)
    assert all(isinstance(x, Fraction) and x > 0 and x.denominator == 1 for x in flux)
    assert _steady(N, flux)


def test_the_uniform_draw_reaches_beyond_the_lp_point():
    """The uniform draw adds integers up to 2^64 to the free coordinates, so its entries are far larger, and it
    depends on the seed."""
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    lp, uniform = kimonu.interior_flux(N, draw="lp"), kimonu.interior_flux(N)
    assert max(uniform) > 2**64 > max(lp)
    assert kimonu.interior_flux(N, seed=1) != uniform
    assert kimonu.interior_flux(N) == uniform


def test_bits_set_the_size_of_the_uniform_draw():
    """r* = K (T c + s) with s below 2^bits and T = ceil(2 · 2^bits · reach / min r_c), so r* has at most about bits
    more bits than the LP point r_c (T is 1 while 2^bits is small against min r_c). bits=64 is the default, and the
    LP point ignores bits."""
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    lp_flux = kimonu.interior_flux(N, draw="lp")
    lp = max(lp_flux).numerator.bit_length()
    for bits in (8, 20, 64):
        flux = kimonu.interior_flux(N, bits=bits)
        assert all(x > 0 and x.denominator == 1 for x in flux) and _steady(N, flux) and flux != lp_flux
        assert max(flux).numerator.bit_length() <= lp + bits + 12
    assert max(kimonu.interior_flux(N, bits=64)).numerator.bit_length() > lp + 40
    assert kimonu.interior_flux(N, bits=64) == kimonu.interior_flux(N)
    assert kimonu.interior_flux(N, draw="lp", bits=8) == kimonu.interior_flux(N, draw="lp")
    with pytest.raises(ValueError):
        kimonu.interior_flux(N, bits=0)


def test_bits_reach_every_route_and_certify():
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    flux = kimonu.interior_flux(N, bits=16)
    assert kimonu.iterated_route(A, Y, N, bits=16).flux == flux == kimonu.tangent_route(A, Y, N, bits=16).flux
    cert = kimonu.certify(A, Y, N, bits=16)
    assert cert.flux == flux and cert.certified


def test_unknown_draw():
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    with pytest.raises(ValueError):
        kimonu.interior_flux(N, draw="sampled")


def test_blocked_network_is_rejected():
    """a -> b with nothing consuming b has no positive steady state."""
    A, Y, N, cx = matrices([("a", "b")], ["a", "b"])
    with pytest.raises(ValueError):
        kimonu.interior_flux(N)


def test_rate_columns_give_a_detached_uptake_its_own_source():
    """0 -> a, a <-> b, b -> 0. Attached, the uptake runs from the zero complex. Detached, from a source after the
    three complexes."""
    A, Y, N, cx = matrices(OPEN, ["a", "b"])
    zero, a, b = cx.index("0"), cx.index("a"), cx.index("b")
    assert rate_columns(A, Y, detach_uptakes=False) == ({0: zero, 1: a, 2: b, 3: b}, 3)
    assert rate_columns(A, Y) == ({0: 3, 1: a, 2: b, 3: b}, 4)


def test_steady_state_map_by_hand():
    """Column c of W sums the columns of N of the reactions that consume c, weighted by r*."""
    A, Y, N, cx = matrices(OPEN, ["a", "b"])
    r = [Fraction(v) for v in (2, 5, 3, 2)]
    W, shape = kimonu.steady_state_map(A, Y, N, r, detach_uptakes=False)
    zero, a, b = cx.index("0"), cx.index("a"), cx.index("b")
    assert shape == (2, 3)
    assert {k: v for k, v in W.items() if v} == {(0, zero): 2, (0, a): -5, (1, a): 5, (0, b): 3, (1, b): -5}
    W, shape = kimonu.steady_state_map(A, Y, N, r)
    assert shape == (2, 4) and W[(0, 3)] == 2 and (0, zero) not in W


def test_w_annihilates_the_reference():
    """u = 1 is the reference steady state, so W·1 = N r* = 0."""
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    W, (m, n_columns) = kimonu.steady_state_map(A, Y, N, kimonu.interior_flux(N))
    assert all(sum(W.get((i, c), 0) for c in range(n_columns)) == 0 for i in range(m))
