"""Formula tests: monomials from flux coupling, symbolic values, and both against the exact values.

The open network has one binomial value (b = k0 / k3) and one that is not (a); two uptakes with coupled
fluxes constrain k when attached. On the toy networks every formula, evaluated at k = r* / ψ(x), must
give the exact value relation_values gives there.
"""
from fractions import Fraction
import random

import pytest
import sympy as sp

import kimonu
from networks import ENVZ, ENVZ_SPECIES, OPEN, matrices
from test_values import TOYS

COUPLED = [("0", "A"), ("0", "B"), ("A+B", "C"), ("C", "0")]


def test_open_network():
    A, Y, N, _ = matrices(OPEN, ["a", "b"])
    b = kimonu.binomial_formulas(A, Y, N, metabolite_ids=["a", "b"], detach_uptakes=False)
    assert b.acr["a"] is None and b.acr["b"].exponents == {0: 1, 3: -1} and b.acr["b"].constant == 1
    s = kimonu.symbolic_formulas(A, Y, N, metabolite_ids=["a", "b"], detach_uptakes=False)
    k0, k1, k2, k3 = sp.symbols("k0:4", positive=True)
    assert sp.simplify(s.acr["a"] - k0 * (k2 + k3) / (k1 * k3)) == 0


def test_coupled_uptakes_constrain_k():
    """Attached, both uptakes run at k0 and k1 with equal fluxes: k0 = k1. Detached, nothing is constrained."""
    A, Y, N, _ = matrices(COUPLED, ["A", "B", "C"])
    b = kimonu.binomial_formulas(A, Y, N, detach_uptakes=False)
    assert [c.exponents for c in b.constraints] in ([{0: 1, 1: -1}], [{0: -1, 1: 1}])
    s = kimonu.symbolic_formulas(A, Y, N, detach_uptakes=False)
    assert len(s.constraints) == 1
    assert kimonu.binomial_formulas(A, Y, N, detach_uptakes=True).constraints == []


def test_envz_yp_formula():
    A, Y, N, _ = matrices(ENVZ, ENVZ_SPECIES)
    yp = kimonu.symbolic_formulas(A, Y, N, metabolite_ids=ENVZ_SPECIES).acr["Yp"]
    k = [3, 2, 5, 2, 8, 8, 8, 7, 4, 2, 8, 1, 7, 7]
    assert yp.subs({sp.Symbol(f"k{j}", positive=True): v for j, v in enumerate(k)}) == sp.Rational(60, 29)


@pytest.mark.parametrize("name", sorted(TOYS))
@pytest.mark.parametrize("detach", [False, True])
def test_formulas_give_the_exact_values(name, detach):
    reactions, species = TOYS[name]
    A, Y, N, _ = matrices(reactions, species)
    flux = kimonu.iterated_route(A, Y, N, detach_uptakes=detach).flux
    rng = random.Random(1)
    k = kimonu.rate_constants(A, Y, flux, [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in species])
    v = kimonu.relation_values(A, Y, N, k, detach_uptakes=detach)
    b = kimonu.binomial_formulas(A, Y, N, detach_uptakes=detach)
    s = kimonu.symbolic_formulas(A, Y, N, detach_uptakes=detach)
    at = {sp.Symbol(f"k{j}", positive=True): sp.Rational(x.numerator, x.denominator) for j, x in enumerate(k)}
    assert all(c(k) == 1 for c in b.constraints)
    for got, sym in ((v.acr, s.acr), (dict(v.acrr), dict(s.acrr)), (dict(v.modules), dict(s.modules))):
        assert got.keys() == sym.keys()
    pairs = [(v.acr[i], b.acr[i], s.acr[i]) for i in v.acr]
    for values, bins, syms in ((v.acrr, b.acrr, s.acrr), (v.modules, b.modules, s.modules)):
        bins, syms = dict(bins), dict(syms)
        pairs += [(x, bins[r][j], syms[r][j]) for r, ratios in values for j, x in ratios.items()]
    for exact, monomial, expression in pairs:
        if monomial is not None:
            assert monomial(k) == exact
        assert float(expression.subs(at)) == pytest.approx(float(exact), rel=1e-12)
