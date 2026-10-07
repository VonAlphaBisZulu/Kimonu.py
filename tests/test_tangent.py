"""Tangent-route tests, and the readout both routes share.

EnvZ-OmpR has one steady state per rate-constant set up to its conservation laws, so both routes must
agree on it exactly. The cooperative kinase of Figure 3 is the case where they part: two separate
steady states at the same rate constants, X equal in both and C, P not. The tangent route calls all
three fixed (local ACR), the iterated route none, and the truth (X alone) sits between them.
"""
from fractions import Fraction

import numpy as np
from scipy import sparse

import kimonu
from kimonu.acr import coupling_signatures, free_directions, readout
from kimonu.io import zero_complexes
from kimonu.linalg import positive_flux
from test_kinetic import _crn, ENVZ, ENVZ_SPECIES, OPEN

# Figure 3: protein X is made and degraded, dimerises into an inactive kinase C, and active P
# activates C cooperatively. Reactions in the order of the rate constants K below.
KINASE = [("0", "X"), ("X", "0"), ("X+X", "C"), ("C+P+P", "P+P+P"), ("C", "0"), ("P", "0")]
KINASE_SPECIES = ["X", "C", "P"]
K = [20, 1, Fraction(1, 20), Fraction(1, 4), 1, 1]

FIG1 = [("A", "B"), ("B", "A"), ("B", "C"), ("C", "B"), ("C", "D"), ("D+E", "F"), ("F", "D+E"), ("F", "B+G"),
        ("C+G", "H"), ("H", "C+G"), ("H", "C+E"), ("A+G", "I"), ("I", "A+G"), ("I", "A+E")]

# 0 <-> A, 2A <-> 3A: three steady states x_A = 1, 2, 3 for suitable rate constants
CUBIC = [("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")]

TOYS = {"envz": (ENVZ, ENVZ_SPECIES), "open": (OPEN, ["a", "b"]), "kinase": (KINASE, KINASE_SPECIES),
        "fig1": (FIG1, list("ABCDEFGHI")), "cubic": (CUBIC, ["A"])}


def _partition(modules):
    return {frozenset(m) for m in modules}


def _pairs(groups):
    return {(a, b) for g in groups for a in g for b in g if a < b}


def test_envz_both_routes_agree():
    """One steady state per class: modules, ACR (Yp) and the 6-clique of ACRR are the iterated ones."""
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    tan = kimonu.tangent_route(A, Y, N, metabolite_ids=ENVZ_SPECIES)
    it = kimonu.kinetic_modules(A, Y, N)
    assert _partition(tan.modules) == _partition(it.modules)
    assert tan.acr == ["Yp"] and len(tan.acrr) == 15
    assert tan.dimension == 2
    assert tan.module_of[cx.index("XD")] == tan.module_of[cx.index("XD+Yp")]


def test_kinase_tangent_route_claims_what_holds_only_locally():
    """At k = K the kinase has the steady states X = 10, (C, P) = (1, 4) and (4, 1). The tangent route calls
    X, C and P fixed, the iterated route none: it misses X, and rightly does not claim C or P."""
    A, Y, N, cx = _crn(KINASE, KINASE_SPECIES)
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
    tan = kimonu.tangent_route(A, Y, N, metabolite_ids=KINASE_SPECIES)
    assert sorted(tan.acr) == ["C", "P", "X"] and tan.dimension == 0
    it = kimonu.kinetic_modules(A, Y, N)
    assert kimonu.acr_acrr(Y, it.modules, metabolite_ids=KINASE_SPECIES, known_acr=it.known_acr).acr == []


def test_iterated_route_is_contained_in_the_tangent_route():
    """iterated ⊆ true ⊆ tangent: every coupled pair, ACR metabolite and ACRR pair the iterated route
    proves is also in the tangent space readout, attached and detached."""
    for name, (R, S) in TOYS.items():
        A, Y, N, cx = _crn(R, S)
        for detach in (False, True):
            it = kimonu.kinetic_modules(A, Y, N, detach_uptakes=detach)
            found = kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr)
            tan = kimonu.tangent_route(A, Y, N, detach_uptakes=detach)
            assert _pairs(it.modules) <= _pairs(tan.modules), (name, detach)
            assert set(found.acr) <= set(tan.acr), (name, detach)
            assert set(found.acrr) <= set(tan.acrr), (name, detach)


def test_tangent_route_detached_uptakes():
    """A detached uptake is a free source: a loses its ACR, and the zero complex is a module of its own."""
    A, Y, N, cx = _crn(OPEN, ["a", "b"])
    attached = kimonu.tangent_route(A, Y, N, metabolite_ids=["a", "b"])
    detached = kimonu.tangent_route(A, Y, N, detach_uptakes=True, metabolite_ids=["a", "b"])
    assert attached.acr == ["a", "b"] and detached.acr == []
    assert {cx.index("0")} in detached.modules
    assert _partition(attached.modules) == {frozenset(range(len(cx)))}


def test_cubic_is_fixed_on_each_branch_only():
    """Each of the three steady states is isolated (T = {0}), so the tangent route calls A fixed; it takes three
    values, and the iterated route reports nothing."""
    A, Y, N, cx = _crn(CUBIC, ["A"])
    assert kimonu.tangent_route(A, Y, N).acr == [0]
    it = kimonu.kinetic_modules(A, Y, N)
    assert kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr).acr == []


def test_readout_reproduces_the_iterated_modules():
    """At the iterated route's final L = ker(YΔᵀ), equal rows of Yᵀ L over all complexes are its modules, the
    zero rows of L its ACR metabolites, and the grouping agrees with merging whole modules by the signature
    of a representative (Proposition S3-4)."""
    for name, (R, S) in TOYS.items():
        A, Y, N, cx = _crn(R, S)
        for detach in (False, True):
            it = kimonu.kinetic_modules(A, Y, N, detach_uptakes=detach)
            skip = zero_complexes(Y) if detach else set()
            found = readout(Y, free_directions(Y, it.modules, it.known_acr), skip)
            assert _partition(found.modules) | {frozenset({c}) for c in skip} == _partition(it.modules), (name, detach)
            assert found.acr == it.known_acr, (name, detach)
            kept = [m for m in it.modules if not (set(m) & skip)]
            signatures = coupling_signatures(Y, kept, it.known_acr)
            if signatures is not None:
                assert len(set(signatures)) == len(kept), (name, detach)


def test_integer_flux_is_an_exact_interior_steady_state():
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    flux = positive_flux(N, scale=10**6)
    assert all(x > 0 and x.denominator == 1 for x in flux)
    Nd = N.toarray()
    for i in range(Nd.shape[0]):
        assert sum(int(Nd[i, j]) * flux[j] for j in range(Nd.shape[1])) == 0


def test_free_directions_are_exact_for_long_decimal_coefficients():
    """Biomass-like coefficients with 15 significant digits: their difference does not fit in a float, and a
    basis built from the rounded difference would separate the two complexes of the module."""
    Y = sparse.csr_matrix(np.array([[9.86915907870449e-05, 4.11265742135554e-05],
                                    [0.000655436261176449, 1.65224841134181e-05],
                                    [1.0, 1.0]]))
    assert readout(Y, free_directions(Y, [{0, 1}])).modules == [{0, 1}]
