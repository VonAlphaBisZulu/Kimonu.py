"""Tangent-route tests.

EnvZ-OmpR has one steady state per rate-constant set up to its conservation laws, so both routes must
agree on it exactly. The cooperative kinase of Figure 2 is the case where they part: two separate
steady states at the same rate constants, X equal in both and C, P not. The tangent route calls all
three fixed (C and P falsely), the iterated route none, and the truth (X alone) sits between them.
"""
from fractions import Fraction

import kimonu
from networks import CUBIC, ENVZ, ENVZ_SPECIES, K, KINASE, KINASE_SPECIES, OPEN, TOYS, network


def _partition(modules):
    return {frozenset(m) for m in modules}


def _pairs(groups):
    return {(a, b) for g in groups for a in g for b in g if a < b}


def test_envz_both_routes_agree():
    """One steady state per compatibility class: modules, ACR (Yp) and the 6-clique of ACRR are the iterated ones."""
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    tan = kimonu.tangent_route(crn)
    it = kimonu.iterated_route(crn)
    assert tan.route == "tangent" and _partition(tan.modules) == _partition(it.modules)
    assert tan.acr == ["Yp"] and len(tan.acrr) == 15 and tan.acrr == it.acrr
    assert tan.dimension == 2
    assert tan.module_of[cx.index("XD")] == tan.module_of[cx.index("XD+Yp")]


def test_kinase_tangent_route_adds_false_relations():
    """At k = K the kinase has the steady states X = 10, (C, P) = (1, 4) and (4, 1). The tangent route calls
    X, C and P fixed (C and P falsely), the iterated route none: it misses X, and rightly does not claim C or P."""
    crn, cx = network(KINASE, KINASE_SPECIES)
    for state in ({"X": 10, "C": 1, "P": 4}, {"X": 10, "C": 4, "P": 1}):
        dx = {s: Fraction(0) for s in KINASE_SPECIES}
        for (sub, prod), k in zip(KINASE, K):
            rate = Fraction(k)
            for s in ([] if sub == "0" else sub.split("+")):
                rate *= state[s]
            for s in ([] if sub == "0" else sub.split("+")):
                dx[s] -= rate
            for s in ([] if prod == "0" else prod.split("+")):
                dx[s] += rate
        assert not any(dx.values()), state
    tan = kimonu.tangent_route(crn, detach_uptakes=False)
    assert sorted(tan.acr) == ["C", "P", "X"] and tan.dimension == 0
    assert kimonu.iterated_route(crn, detach_uptakes=False).acr == []


def test_iterated_route_is_contained_in_the_tangent_route():
    """iterated ⊆ true ⊆ tangent: every coupled pair, ACR metabolite and ACRR pair the iterated route
    proves is also in the tangent space readout, attached and detached."""
    for name, (R, S) in TOYS.items():
        crn, cx = network(R, S)
        for detach in (False, True):
            it = kimonu.iterated_route(crn, detach_uptakes=detach)
            tan = kimonu.tangent_route(crn, detach_uptakes=detach)
            assert _pairs(it.modules) <= _pairs(tan.modules), (name, detach)
            assert set(it.acr) <= set(tan.acr), (name, detach)
            assert set(it.acrr) <= set(tan.acrr), (name, detach)


def test_tangent_route_detached_uptakes():
    """A detached uptake is a free source: a loses its ACR, and the zero complex is a module of its own."""
    crn, cx = network(OPEN, ["a", "b"])
    attached = kimonu.tangent_route(crn, detach_uptakes=False)
    detached = kimonu.tangent_route(crn, detach_uptakes=True)
    assert attached.acr == ["a", "b"] and detached.acr == []
    assert detached.modules[-1] == {cx.index("0")}
    assert _partition(attached.modules) == {frozenset(range(len(cx)))}


def test_cubic_is_fixed_on_each_branch_only():
    """Each of the three steady states is isolated (T = {0}), so the tangent route calls A fixed; it takes three
    values, and the iterated route reports nothing."""
    crn, cx = network(CUBIC, ["A"], names=False)
    assert kimonu.tangent_route(crn, detach_uptakes=False).acr == [0]
    assert kimonu.iterated_route(crn, detach_uptakes=False).acr == []


def test_tangent_space_has_the_source_coordinates_after_the_metabolites():
    """Detached, the open network's uptake gets a source coordinate after a and b. The steady states move along
    x = t·1 with source t, so ker(W Yᵀ) is one-dimensional, equal in all three coordinates, and B is its first two
    rows."""
    crn, cx = network(OPEN, ["a", "b"])
    flux = kimonu.interior_flux(crn)
    rows = kimonu.tangent_space(crn, flux)
    assert len(rows) == 3 and len({r for r in rows}) == 1 and rows[0]
    tan = kimonu.tangent_route(crn, flux=flux)
    assert tan.B == rows[:2] and tan.dimension == 1
    assert len(kimonu.tangent_space(crn, flux, detach_uptakes=False)) == 2


def test_tangent_route_does_not_depend_on_the_draw():
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    first = kimonu.tangent_route(crn)
    for draw in ("uniform", "lp"):
        for seed in (1, 2):
            other = kimonu.tangent_route(crn, seed=seed, draw=draw)
            assert _partition(other.modules) == _partition(first.modules)
            assert (other.acr, other.acrr) == (first.acr, first.acrr)
