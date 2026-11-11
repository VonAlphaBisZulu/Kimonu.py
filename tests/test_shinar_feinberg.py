"""Networks of Shinar and Feinberg with known answers.

Shinar & Feinberg (2010) give two biochemical networks with ACR for every choice of rate constants, by their
deficiency-one criterion (two non-terminal complexes that differ in one species). EnvZ-OmpR, here in their form with
seven species (test_kinetic.py has the form with ADP-bound EnvZ of Figure 1), has ACR in the phosphorylated response
regulator Yp. The IDHKP-IDH core, with the bifunctional kinase-phosphatase E (X1), phosphorylated IDH (X2) and IDH (X4),
has ACR in X4. Networks 11 and 12 are the concordant network and the EnvZ-OmpR variant that the injectivity test of
benchmarks/kinetic_variants/test_tangent_one_shot.py uses. Network 12 has ACR in F by the same criterion (B and B+F),
network 11 has none.

On all four both routes agree, so the lists are exact, and the tests also pin the ACRR pairs and the kinetic modules.
"""
import pytest

import kimonu
from test_kinetic import _crn


def _rev(a, b):
    return [(a, b), (b, a)]


NETWORKS = {
    "envz_ompr": (_rev("X", "XT") + [("XT", "Xp")] + _rev("Xp+Y", "XpY") + [("XpY", "X+Yp")]
                  + _rev("XT+Yp", "XTYp") + [("XTYp", "XT+Y")],
                  ["X", "XT", "Xp", "Y", "Yp", "XpY", "XTYp"]),
    "idhkp_idh": (_rev("X1+X2", "X3") + [("X3", "X1+X4")] + _rev("X3+X4", "X5") + [("X5", "X3+X2")],
                  ["X1", "X2", "X3", "X4", "X5"]),
    "network_11": ([("A", "B"), ("B", "A"), ("C+D+D+D+D", "E"), ("E", "C+D+D+D+D"), ("B+E", "F"), ("F", "A+E"),
                    ("A+E", "F"), ("B", "D+B"), ("D", "0")], list("ABCDEF")),
    "network_12": ([("A", "B"), ("B", "A"), ("B", "C"), ("C+D", "E"), ("E", "C+D"), ("E", "A+F"), ("B+F", "G"),
                    ("G", "B+F"), ("G", "B+D")], list("ABCDEFG")),
}

# ACR, ACRR groups (all pairs within a group have ACRR) and kinetic modules of two or more complexes
EXPECTED = {
    "envz_ompr": (["Yp"], [{"X", "XT", "XpY", "XTYp"}],
                  [{"X", "X+Yp", "XT", "XT+Yp", "XTYp", "Xp+Y", "XpY"}]),
    "idhkp_idh": (["X4"], [{"X3", "X5"}], [{"X1+X2", "X3", "X3+X4", "X5"}]),
    "network_11": ([], [{"B", "D"}], [{"B", "D"}, {"C+D+D+D+D", "E"}]),
    "network_12": (["F"], [{"A", "B", "E", "G"}], [{"A", "A+F", "B", "B+F", "C+D", "E", "G"}]),
}


def _pairs(groups):
    return sorted(tuple(sorted((a, b))) for g in groups for a in g for b in g if a < b)


@pytest.mark.parametrize("name", sorted(NETWORKS))
@pytest.mark.parametrize("detach", [False, True])
def test_both_routes_give_the_known_answer(name, detach):
    reactions, species = NETWORKS[name]
    acr, acrr_groups, modules = EXPECTED[name]
    A, Y, N, cx = _crn(reactions, species)
    it = kimonu.kinetic_modules(A, Y, N, detach_uptakes=detach)
    res = kimonu.acr_acrr(Y, it.modules, metabolite_ids=species, known_acr=it.known_acr)
    tan = kimonu.tangent_route(A, Y, N, metabolite_ids=species, detach_uptakes=detach)
    named = lambda mods: sorted(sorted(cx[i] for i in m) for m in mods if len(m) > 1)
    assert sorted(res.acr) == acr and sorted(tan.acr) == acr, name
    assert sorted(tuple(sorted(p)) for p in res.acrr) == _pairs(acrr_groups), name
    assert sorted(tuple(sorted(p)) for p in tan.acrr) == _pairs(acrr_groups), name
    assert named(it.modules) == named(tan.modules) == sorted(sorted(m) for m in modules), name
