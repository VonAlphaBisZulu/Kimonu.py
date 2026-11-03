"""Value tests: the exact values of the relations at given rate constants, and the rate constants at given
values.

Closed forms pin the open network (rational values) and 0 -> X, 2X -> 0 (a square root). The round trip
checks every relation on the toy networks: at k = r* / ψ(x) the steady state is x, so each ACR value
must come out as x_i, each ACRR ratio as x_j / x_i and each module ratio as ψ_d(x) / ψ_c(x), exactly.
EnvZ-OmpR, whose Yp value needs the feedback, is also checked against a simulation.
"""
from fractions import Fraction
import random

import numpy as np
import pytest
from scipy import sparse
from scipy.integrate import solve_ivp

import kimonu
from test_kinetic import _crn, ENVZ, ENVZ_SPECIES, OPEN
from test_tangent import FIG1, KINASE, KINASE_SPECIES

SQRT = [("0", "X"), ("X+X", "0")]
TOYS = {"envz": (ENVZ, ENVZ_SPECIES), "open": (OPEN, ["a", "b"]), "kinase": (KINASE, KINASE_SPECIES),
        "fig1": (FIG1, list("ABCDEFGHI"))}


def _psi(Y, c, x):
    col = sparse.csc_matrix(Y)[:, c]
    out = Fraction(1)
    for i, y in zip(col.indices, col.data):
        out *= x[i] ** int(y)
    return out


def test_open_network_closed_form():
    """0 -> a, a <-> b, b -> 0: b = k1 / k4 and a = k1 (k3 + k4) / (k2 k4); detached only a / b is left."""
    A, Y, N, _ = _crn(OPEN, ["a", "b"])
    k = [2, 3, 5, 7]
    v = kimonu.relation_values(A, Y, N, k, metabolite_ids=["a", "b"], detach_uptakes=False)
    assert v.acr == {"a": Fraction(2 * 12, 3 * 7), "b": Fraction(2, 7)}
    v = kimonu.relation_values(A, Y, N, k, detach_uptakes=True, metabolite_ids=["a", "b"])
    assert v.acr == {} and v.acrr == [("a", {"b": Fraction(3, 12)})]


def test_square_root_value():
    """dx/dt = k1 - 2 k2 x^2: x = (k1 / 2 k2)^(1/2), irrational for k = (3, 5) and exact as a Value."""
    A, Y, N, _ = _crn(SQRT, ["X"])
    x = kimonu.relation_values(A, Y, N, [3, 5], metabolite_ids=["X"], detach_uptakes=False).acr["X"]
    assert x.rational is None and (x.r, x.b, x.q) == (1, Fraction(3, 10), 2)
    assert float(x) == pytest.approx((3 / 10) ** 0.5, rel=1e-15)
    assert kimonu.relation_values(A, Y, N, [8, 1], metabolite_ids=["X"], detach_uptakes=False).acr["X"] == 2


@pytest.mark.parametrize("name", sorted(TOYS))
@pytest.mark.parametrize("detach", [False, True])
def test_round_trip(name, detach):
    """At k = r* / ψ(x) every value is the one x gives it, and the relations are those of the iterated route."""
    reactions, species = TOYS[name]
    A, Y, N, _ = _crn(reactions, species)
    kin = kimonu.kinetic_modules(A, Y, N, detach_uptakes=detach)
    rng = random.Random(0)
    x = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in species]
    k = kimonu.rate_constants(A, Y, kin.flux, x)
    v = kimonu.relation_values(A, Y, N, k, detach_uptakes=detach)
    assert sorted(v.acr) == kin.known_acr
    assert v.acr == {i: x[i] for i in v.acr}
    for ref, ratios in v.acrr:
        assert ratios == {j: x[j] / x[ref] for j in ratios}
    assert {frozenset([c, *r]) for c, r in v.modules} == {frozenset(m) for m in kin.modules if len(m) > 1}
    for c, ratios in v.modules:
        assert ratios == {d: _psi(Y, d, x) / _psi(Y, c, x) for d in ratios}


def test_envz_yp_matches_a_simulation():
    """Yp from the exact value and from integrating the ODE to steady state, for two different totals."""
    A, Y, N, _ = _crn(ENVZ, ENVZ_SPECIES)
    k = [3, 2, 5, 2, 8, 8, 8, 7, 4, 2, 8, 1, 7, 7]
    v = kimonu.relation_values(A, Y, N, k, metabolite_ids=ENVZ_SPECIES)
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
    A, Y, N, _ = _crn(OPEN, ["a", "b"])
    flux = kimonu.kinetic_modules(A, Y, N).flux
    assert kimonu.rate_constants(A, Y, flux) == flux


def test_wrong_k_is_rejected():
    A, Y, N, _ = _crn(OPEN, ["a", "b"])
    with pytest.raises(ValueError):
        kimonu.relation_values(A, Y, N, [1, 2, 3])
    with pytest.raises(ValueError):
        kimonu.relation_values(A, Y, N, [1, 2, 0, 4])
