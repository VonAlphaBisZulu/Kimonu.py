"""The shared readout: ACR, ACRR and kinetic modules as zero rows and equal rows of a basis B.

EnvZ-OmpR is the canonical Shinar-Feinberg example. The phosphorylated response regulator Yp has absolute
concentration robustness and nothing else does. The tests supply the kinetic modules directly, so they pin the
readout independently of how the modules are derived.
"""
from fractions import Fraction

import numpy as np
import kimonu
from kimonu.linalg import matrix_product
from kimonu.readout import Readout, composition_differences, route_result
from networks import ENVZ, ENVZ_SPECIES, OPEN, TOYS, network

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

# The kinetic modules of this network; the large one is the largest module.
KINETIC_MODULES = [{0, 1, 2, 3, 5, 6, 8, 10, 11}, {9}, {4}, {7}, {12}]


def _crn(Y=None):
    """The CRN of the compositions alone, without reactions, which is all the readout reads."""
    if Y is None:
        index = {m: i for i, m in enumerate(METABOLITES)}
        Y = np.zeros((len(METABOLITES), len(COMPLEXES)), dtype=np.int64)
        for j, comp in enumerate(COMPLEXES):
            for met, coeff in comp.items():
                Y[index[met], j] = coeff
    return kimonu.crn_from_matrices(np.zeros((Y.shape[1], 0), dtype=np.int64), Y)


def _named(found):
    return [METABOLITES[i] for i in found.acr], [[METABOLITES[i] for i in g] for g in found.acrr_groups]


def test_acr_is_exactly_yp():
    """Yp has ACR and no other species does."""
    acr, _ = _named(kimonu.readout(_crn(), kimonu.free_directions(_crn(), KINETIC_MODULES)))
    assert acr == ["Yp"]


def test_acrr_is_one_group_of_six():
    """The ratio-robust metabolites form one group of six, 15 pairs. ACR metabolites are left out of the groups:
    a constant ratio between two constants is not an independent finding."""
    _, groups = _named(kimonu.readout(_crn(), kimonu.free_directions(_crn(), KINETIC_MODULES)))
    assert [sorted(g) for g in groups] == [sorted(["X", "XD", "XDYp", "XT", "XTYp", "XpY"])]


def test_no_modules_means_no_relation():
    """With no coupled complexes there is nothing to span, B is the identity, and nothing is fixed."""
    B = kimonu.free_directions(_crn(), [{i} for i in range(len(COMPLEXES))])
    assert B == [((i, Fraction(1)),) for i in range(len(METABOLITES))]
    found = kimonu.readout(_crn(), B)
    assert found.acr == [] and found.acrr_groups == [] and len(found.modules) == len(COMPLEXES)


def test_known_acr_enters_as_a_unit_vector():
    """A metabolite known to have ACR gets a zero row of B even without a module that shows it."""
    B = kimonu.free_directions(_crn(), [], known_acr=[METABOLITES.index("Y")])
    assert kimonu.readout(_crn(), B).acr == [METABOLITES.index("Y")]


def test_composition_differences_against_the_first_member():
    """(YΔ)ᵀ has one row y_c - y_d against the first member d of each module, and a unit row per known ACR
    metabolite."""
    rows = composition_differences(_crn(), [{1, 12}, {4}, set()], known_acr=[0])
    assert rows == [((METABOLITES.index("Yp"), -1),), ((0, 1),)]


def test_equal_rows_of_yt_b_are_the_modules():
    """Row c of Yᵀ B adds the rows of B of complex c's metabolites, and equal rows are the readout's modules."""
    crn = _crn()
    B = kimonu.free_directions(crn, KINETIC_MODULES)
    rows = matrix_product(crn.Yt, B)
    assert len(rows) == len(COMPLEXES)
    groups = {}
    for c, row in enumerate(rows):
        groups.setdefault(row, set()).add(c)
    assert sorted(map(sorted, groups.values())) == sorted(map(sorted, kimonu.readout(crn, B).modules))


def test_detached_readout_lists_the_zero_complex_last():
    """0 -> a, a <-> b, b -> 0 with all three complexes coupled, as attached. a and b have ACR, so the rows of Yᵀ B
    of a, b and the zero complex are all zero. Detached, the zero complex is not grouped and comes last."""
    crn, cx = network(OPEN, ["a", "b"], names=False)
    B = kimonu.free_directions(crn, [{0, 1, 2}])
    assert kimonu.readout(crn, B, detach_uptakes=False).modules == [{0, 1, 2}]
    assert kimonu.readout(crn, B).modules == [{1, 2}, {0}] and kimonu.readout(crn, B).acr == [0, 1]


def test_route_result_names_metabolites_and_gives_every_acrr_pair():
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    found = Readout(modules=[{0, 1}] + [{c} for c in range(2, len(cx))], acr=[6], acrr_groups=[[3, 0, 1]])
    res = route_result("iterated", crn, found, dimension=1)
    assert res.acr == ["Yp"] and res.acrr == [("XD", "X"), ("Xp", "X"), ("Xp", "XD")]
    assert res.acrr_groups == [["XD", "X", "Xp"]]
    assert res.n_modules == 1 and res.largest_module == 2 and list(res.module_of[:3]) == [1, 1, 2]


def test_acrr_groups_join_pairs_that_share_a_metabolite():
    """The reference route's Shinar–Feinberg pairs (0, 1) and (1, 3) give one group of three, while acrr keeps the
    pairs as they are. EnvZ-OmpR's iterated route has one group of six, the 15 pairs of acrr."""
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    found = Readout(modules=[{c} for c in range(len(cx))], acr=[], acrr_groups=[[0, 1], [1, 3], [4, 5]])
    res = route_result("reference", crn, found)
    assert res.acrr == [("X", "Xp"), ("XD", "X"), ("Y", "XpY")]
    assert res.acrr_groups == [["XD", "X", "Xp"], ["Y", "XpY"]]
    it = kimonu.iterated_route(crn)
    assert len(it.acrr_groups) == 1 and len(it.acrr_groups[0]) == 6 and len(it.acrr) == 15
    assert sum(len(g) - 1 for g in it.acrr_groups) == 5


def test_readout_reproduces_the_iterated_route():
    """At the iterated route's final B = ker((YΔ)ᵀ), equal rows of Yᵀ B over all complexes are its modules and the zero
    rows of B its ACR metabolites, also when the modules and ACR are given back to free_directions."""
    for name, (R, S) in TOYS.items():
        crn, cx = network(R, S, names=False)
        for detach in (False, True):
            it = kimonu.iterated_route(crn, detach_uptakes=detach)
            again = kimonu.readout(crn, kimonu.free_directions(crn, it.modules, it.acr), detach)
            assert again.modules == it.modules, (name, detach)
            assert again.acr == it.acr, (name, detach)


def test_route_result_layout_is_shared_by_the_three_routes():
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    results = [kimonu.iterated_route(crn), kimonu.tangent_route(crn), kimonu.reference_route(crn)]
    assert [r.route for r in results] == ["iterated", "tangent", "reference"]
    assert all(type(r) is kimonu.RouteResult for r in results)
    it, tan, ref = results
    assert it.flux and tan.flux and ref.flux is None and ref.B is None and tan.rounds is None
    assert ref.concordance.n_balanced == 5


def test_free_directions_are_exact_for_long_decimal_coefficients():
    """Biomass-like coefficients with 15 significant digits: their difference does not fit in a float, and a basis
    built from the rounded difference would separate the two complexes of the module."""
    crn = _crn(np.array([[9.86915907870449e-05, 4.11265742135554e-05],
                         [0.000655436261176449, 1.65224841134181e-05],
                         [1.0, 1.0]]))
    assert kimonu.readout(crn, kimonu.free_directions(crn, [{0, 1}])).modules == [{0, 1}]
