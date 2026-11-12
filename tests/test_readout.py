"""The shared readout: ACR, ACRR and kinetic modules as zero rows and equal rows of a basis B.

EnvZ-OmpR is the canonical Shinar-Feinberg example. The phosphorylated response regulator Yp has absolute
concentration robustness and nothing else does. The tests supply the kinetic modules directly, so they pin the
readout independently of how the modules are derived.
"""
from fractions import Fraction

import numpy as np
from scipy import sparse

import kimonu
from kimonu.readout import complex_directions, composition_differences
from networks import ENVZ, ENVZ_SPECIES, TOYS, matrices

# EnvZ-OmpR as COCOA constructs it: 9 species, 13 complexes.
METABOLITES = ["X", "XD", "XDYp", "XT", "XTYp", "Xp", "XpY", "Y", "Yp"]

# complex index -> composition, in the model's complex order
COMPLEXES = [
    {"XD": 1},                 # 0
    {"X": 1},                  # 1
    {"XTYp": 1},               # 2
    {"XT": 1, "Yp": 1},        # 3
    {"XT": 1, "Y": 1},         # 4
    {"XD": 1, "Yp": 1},        # 5
    {"XDYp": 1},               # 6
    {"XD": 1, "Y": 1},         # 7
    {"XT": 1},                 # 8
    {"Xp": 1},                 # 9
    {"Xp": 1, "Y": 1},         # 10
    {"XpY": 1},                # 11
    {"X": 1, "Yp": 1},         # 12
]

# The kinetic modules of this network; the large one is the giant module.
KINETIC_MODULES = [{0, 1, 2, 3, 5, 6, 8, 10, 11}, {9}, {4}, {7}, {12}]


def _Y():
    index = {m: i for i, m in enumerate(METABOLITES)}
    rows, cols, vals = [], [], []
    for j, comp in enumerate(COMPLEXES):
        for met, coeff in comp.items():
            rows.append(index[met])
            cols.append(j)
            vals.append(coeff)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(len(METABOLITES), len(COMPLEXES)), dtype=np.int64)


def _named(found):
    return [METABOLITES[i] for i in found.acr], [[METABOLITES[i] for i in g] for g in found.acrr_groups]


def test_acr_is_exactly_yp():
    """Yp has ACR and no other species does."""
    acr, _ = _named(kimonu.readout(_Y(), kimonu.free_directions(_Y(), KINETIC_MODULES)))
    assert acr == ["Yp"]


def test_acrr_is_one_group_of_six():
    """The ratio-robust metabolites form one group of six, 15 pairs. ACR metabolites are left out of the groups:
    a constant ratio between two constants is not an independent finding."""
    _, groups = _named(kimonu.readout(_Y(), kimonu.free_directions(_Y(), KINETIC_MODULES)))
    assert [sorted(g) for g in groups] == [sorted(["X", "XD", "XDYp", "XT", "XTYp", "XpY"])]


def test_no_modules_means_no_relation():
    """With no coupled complexes there is nothing to span, B is the identity, and nothing is fixed."""
    B = kimonu.free_directions(_Y(), [{i} for i in range(len(COMPLEXES))])
    assert B == [((i, Fraction(1)),) for i in range(len(METABOLITES))]
    found = kimonu.readout(_Y(), B)
    assert found.acr == [] and found.acrr_groups == [] and len(found.modules) == len(COMPLEXES)


def test_known_acr_enters_as_a_unit_vector():
    """A metabolite known to have ACR gets a zero row of B even without a module that shows it."""
    B = kimonu.free_directions(_Y(), [], known_acr=[METABOLITES.index("Y")])
    assert kimonu.readout(_Y(), B).acr == [METABOLITES.index("Y")]


def test_composition_differences_against_the_first_member():
    """YΔ holds y_c - y_d against the first member d of each module, and a unit vector per known ACR metabolite."""
    Y = _Y()
    YD = composition_differences(Y, [{1, 12}, {4}], known_acr=[0])
    assert YD.dtype == np.int64 and YD.shape == (len(METABOLITES), 2)
    assert {METABOLITES[i]: int(v) for i, v in zip(YD[:, 0].nonzero()[0], YD[:, 0].data)} == {"Yp": -1}
    assert YD[:, 1].toarray().ravel().tolist() == [1] + [0] * (len(METABOLITES) - 1)


def test_complex_directions_are_the_rows_of_yt_b():
    """Row c of Yᵀ B adds the rows of B of complex c's metabolites, and equal rows are the readout's modules."""
    Y = _Y()
    B = kimonu.free_directions(Y, KINETIC_MODULES)
    rows = complex_directions(Y, B)
    assert len(rows) == len(COMPLEXES)
    groups = {}
    for c, row in enumerate(rows):
        groups.setdefault(row, set()).add(c)
    assert sorted(map(sorted, groups.values())) == sorted(map(sorted, kimonu.readout(Y, B).modules))


def test_skipped_complexes_are_left_out_of_the_modules():
    found = kimonu.readout(_Y(), kimonu.free_directions(_Y(), KINETIC_MODULES), skip={9})
    assert all(9 not in m for m in found.modules)


def test_readout_reproduces_the_iterated_route():
    """At the iterated route's final B = ker((YΔ)ᵀ), equal rows of Yᵀ B over all complexes are its modules and the zero
    rows of B its ACR metabolites, also when the modules and ACR are given back to free_directions."""
    for name, (R, S) in TOYS.items():
        A, Y, N, cx = matrices(R, S)
        for detach in (False, True):
            it = kimonu.iterated_route(A, Y, N, detach_uptakes=detach)
            skip = kimonu.io.zero_complexes(Y) if detach else set()
            again = kimonu.readout(Y, kimonu.free_directions(Y, it.modules, it.acr), skip)
            assert {frozenset(m) for m in again.modules} | {frozenset({c}) for c in skip} == \
                {frozenset(m) for m in it.modules}, (name, detach)
            assert again.acr == it.acr, (name, detach)


def test_route_result_layout_is_shared_by_the_three_routes():
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    results = [kimonu.iterated_route(A, Y, N), kimonu.tangent_route(A, Y, N), kimonu.reference_route(A, Y, N)]
    assert [r.route for r in results] == ["iterated", "tangent", "reference"]
    assert all(type(r) is kimonu.RouteResult for r in results)
    it, tan, ref = results
    assert it.flux and tan.flux and ref.flux is None and ref.B is None and tan.rounds is None
    assert ref.concordance.n_balanced == 5


def test_free_directions_are_exact_for_long_decimal_coefficients():
    """Biomass-like coefficients with 15 significant digits: their difference does not fit in a float, and a basis
    built from the rounded difference would separate the two complexes of the module."""
    Y = sparse.csr_matrix(np.array([[9.86915907870449e-05, 4.11265742135554e-05],
                                    [0.000655436261176449, 1.65224841134181e-05],
                                    [1.0, 1.0]]))
    assert kimonu.readout(Y, kimonu.free_directions(Y, [{0, 1}])).modules == [{0, 1}]
