"""Value tests: the exact values of the relations at given rate constants, the rate constants at given values, and the
monomial formulas in k from flux coupling.

Closed forms pin the open network (rational values) and 0 -> X, 2X -> 0 (a square root). The round trip
checks every relation on the toy networks: at k = r* / ψ(x) the steady state is x, so each ACR value
must come out as x_i, each ACRR ratio as x_j / x_i and each module ratio as ψ_d(x) / ψ_c(x), exactly.
EnvZ-OmpR, whose Yp value needs the feedback, is also checked against a simulation. The open network has one
binomial value (b = k_0 / k_3) and one that is not (a), and two uptakes with coupled fluxes constrain k when attached.
Every monomial, evaluated at k = r* / ψ(x), must give the exact value relation_values gives there.
"""
from fractions import Fraction
import random

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import kimonu
from networks import ENVZ, ENVZ_SPECIES, FIG1, KINASE, KINASE_SPECIES, OPEN, matrices, network

SQRT = [("0", "X"), ("X+X", "0")]
TOYS = {"envz": (ENVZ, ENVZ_SPECIES), "open": (OPEN, ["a", "b"]), "kinase": (KINASE, KINASE_SPECIES),
        "fig1": (FIG1, list("ABCDEFGHI"))}


def _psi(crn, c, x):
    out = Fraction(1)
    for i, y in crn.Yt[c]:
        out *= x[i] ** int(y)
    return out


def test_open_network_closed_form():
    """0 -> a, a <-> b, b -> 0: b = k1 / k4 and a = k1 (k3 + k4) / (k2 k4); detached only a / b is left."""
    crn, _ = network(OPEN, ["a", "b"])
    k = [2, 3, 5, 7]
    v = kimonu.relation_values(crn, k, detach_uptakes=False)
    assert v.acr == {"a": Fraction(2 * 12, 3 * 7), "b": Fraction(2, 7)}
    v = kimonu.relation_values(crn, k, detach_uptakes=True)
    assert v.acr == {} and v.acrr == [("a", {"b": Fraction(3, 12)})]


def test_square_root_value():
    """dx/dt = k1 - 2 k2 x^2: x = (k1 / 2 k2)^(1/2), irrational for k = (3, 5) and exact as a Value."""
    crn, _ = network(SQRT, ["X"])
    x = kimonu.relation_values(crn, [3, 5], detach_uptakes=False).acr["X"]
    assert x.rational is None and (x.r, x.b, x.q) == (1, Fraction(3, 10), 2)
    assert float(x) == pytest.approx((3 / 10) ** 0.5, rel=1e-15)
    assert kimonu.relation_values(crn, [8, 1], detach_uptakes=False).acr["X"] == 2


@pytest.mark.parametrize("name", sorted(TOYS))
@pytest.mark.parametrize("detach", [False, True])
def test_round_trip(name, detach):
    """At k = r* / ψ(x) every value is the one x gives it, and the relations are those of the iterated route."""
    reactions, species = TOYS[name]
    crn, _ = network(reactions, species, names=False)
    kin = kimonu.iterated_route(crn, detach_uptakes=detach)
    rng = random.Random(0)
    x = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in species]
    k = kimonu.rate_constants(crn, kin.flux, x)
    v = kimonu.relation_values(crn, k, detach_uptakes=detach)
    assert sorted(v.acr) == kin.acr
    assert v.acr == {i: x[i] for i in v.acr}
    for ref, ratios in v.acrr:
        assert ratios == {j: x[j] / x[ref] for j in ratios}
    assert {frozenset([c, *r]) for c, r in v.modules} == {frozenset(m) for m in kin.modules if len(m) > 1}
    for c, ratios in v.modules:
        assert ratios == {d: _psi(crn, d, x) / _psi(crn, c, x) for d in ratios}


def test_envz_yp_matches_a_simulation():
    """Yp from the exact value and from integrating the ODE to steady state, for two different totals."""
    A, Y, N, _ = matrices(ENVZ, ENVZ_SPECIES)
    k = [3, 2, 5, 2, 8, 8, 8, 7, 4, 2, 8, 1, 7, 7]
    v = kimonu.relation_values(network(ENVZ, ENVZ_SPECIES)[0], k)
    assert v.rounds == 2 and list(v.acr) == ["Yp"]
    Nd, S = N.toarray(), -np.minimum(A.toarray(), 0)
    Yd = Y.toarray()

    def rhs(_, z):
        x = np.exp(z)
        return Nd @ (np.array(k, float) * np.prod(x[:, None] ** (Yd @ S), axis=0)) / x

    for start in (np.ones(9), np.linspace(0.5, 3, 9)):
        z = solve_ivp(rhs, (0, 400), np.log(start), method="LSODA", rtol=1e-10, atol=1e-12).y[:, -1]
        assert np.exp(z)[ENVZ_SPECIES.index("Yp")] == pytest.approx(float(v.acr["Yp"]), rel=1e-6)


def test_rate_constants_default_to_the_flux():
    crn, _ = network(OPEN, ["a", "b"])
    flux = kimonu.iterated_route(crn).flux
    assert kimonu.rate_constants(crn, flux) == flux


def test_wrong_k_is_rejected():
    crn, _ = network(OPEN, ["a", "b"])
    with pytest.raises(ValueError):
        kimonu.relation_values(crn, [1, 2, 3])
    with pytest.raises(ValueError):
        kimonu.relation_values(crn, [1, 2, 0, 4])


COUPLED = [("0", "A"), ("0", "B"), ("A+B", "C"), ("C", "0")]


def test_binomial_formulas_of_the_open_network():
    crn, _ = network(OPEN, ["a", "b"])
    b = kimonu.binomial_formulas(crn, detach_uptakes=False)
    assert b.acr["a"] is None and b.acr["b"].exponents == {0: 1, 3: -1} and b.acr["b"].constant == 1
    assert b.acr["b"].format(["up", "ab", "ba", "out"]) == "k_up·k_out^-1"


def test_coupled_uptakes_constrain_k():
    """Attached, both uptakes run at k0 and k1 with equal fluxes: k0 = k1. Detached, nothing is constrained."""
    crn, _ = network(COUPLED, ["A", "B", "C"])
    b = kimonu.binomial_formulas(crn, detach_uptakes=False)
    assert [c.exponents for c in b.constraints] in ([{0: 1, 1: -1}], [{0: -1, 1: 1}])
    assert kimonu.binomial_formulas(crn, detach_uptakes=True).constraints == []


@pytest.mark.parametrize("name", sorted(TOYS))
@pytest.mark.parametrize("detach", [False, True])
def test_monomials_give_the_exact_values(name, detach):
    reactions, species = TOYS[name]
    crn, _ = network(reactions, species, names=False)
    flux = kimonu.iterated_route(crn, detach_uptakes=detach).flux
    rng = random.Random(1)
    k = kimonu.rate_constants(crn, flux, [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in species])
    v = kimonu.relation_values(crn, k, detach_uptakes=detach)
    b = kimonu.binomial_formulas(crn, detach_uptakes=detach)
    assert all(c(k) == 1 for c in b.constraints)
    assert v.acr.keys() == b.acr.keys() and dict(v.acrr).keys() == dict(b.acrr).keys()
    assert dict(v.modules).keys() == dict(b.modules).keys()
    pairs = [(v.acr[i], b.acr[i]) for i in v.acr]
    for values, monomials in ((v.acrr, dict(b.acrr)), (v.modules, dict(b.modules))):
        pairs += [(x, monomials[r][j]) for r, ratios in values for j, x in ratios.items()]
    for exact, monomial in pairs:
        if monomial is not None:
            assert monomial(k) == exact
