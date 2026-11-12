"""Flux coupling: fully coupled reactions under mass balance, and the kinetic modules lifted into flux space.

The small open network is checked by hand in both uptake modes. EnvZ-OmpR has all nine substrate
complexes in one kinetic module, so at fixed rate constants every steady-state flux is a multiple of r*,
where flux balance couples only two pairs. The cooperative kinase has two steady states at the same rate
constants, which shows what the iterated route proves, what the tangent route claims only near one steady
state, and that the truth sits between them.
"""
from fractions import Fraction

import pytest

import kimonu
from kimonu.linalg import nullspace
from kimonu.readout import row_readout
from networks import ENVZ, ENVZ_SPECIES, K, KINASE, KINASE_SPECIES, OPEN, TOYS, network


def _pairs(groups):
    return {(a, b) for g in groups for a in g for b in g if a < b}


def _flux_balance(crn):
    """Full coupling under flux balance: the proportional rows of K = ker(N)."""
    return row_readout(nullspace(crn.N)).groups


def _rates(state, k=K, reactions=KINASE):
    """Mass-action rates of the kinase at a concentration state."""
    out = []
    for (sub, _), kj in zip(reactions, k):
        rate = Fraction(kj)
        for s in ([] if sub == "0" else sub.split("+")):
            rate *= state[s]
        out.append(rate)
    return out


def test_open_network_by_hand():
    """0 -> a <-> b -> 0. Flux balance couples only the uptake and the export (r1 = r4). Under mass action
    the steady state fixes x_a / x_b = (k3 + k4) / k2, so a and b form a module and all four rates are
    coupled. Detached, the supply is free and no rate is fixed; attached, a and b have ACR and every rate
    is fixed."""
    crn, cx = network(OPEN, ["a", "b"])
    assert _pairs(_flux_balance(crn)) == {(0, 3)}
    for route in ("iterated", "tangent"):
        detached = kimonu.kinetic_flux_coupling(crn, route=route)
        assert detached.groups == [{0, 1, 2, 3}] and detached.fixed_rate == [], route
        attached = kimonu.kinetic_flux_coupling(crn, detach_uptakes=False, route=route)
        assert attached.groups == [{0, 1, 2, 3}] and attached.fixed_rate == [0, 1, 2, 3], route
        assert detached.n_groups == attached.n_groups == 1 and detached.largest_group == 4


def test_envz_every_rate_in_one_group():
    """All substrate complexes of EnvZ-OmpR lie in one kinetic module, so every rate keeps a fixed ratio to
    every other, while flux balance couples only XD <-> X and XT -> Xp with XpY -> X + Yp."""
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    assert sorted(sorted(g) for g in _flux_balance(crn) if len(g) > 1) == [[0, 1], [4, 7]]
    for route in ("iterated", "tangent"):
        found = kimonu.kinetic_flux_coupling(crn, route=route)
        assert found.n_groups == 1 and found.largest_group == len(ENVZ), route
        assert list(found.group_of) == [1] * len(ENVZ) and found.fixed_rate == [], route


def test_flux_balance_and_both_routes_nest():
    """flux balance ⊆ iterated ⊆ tangent for the coupled pairs, iterated ⊆ tangent for the fixed rates,
    and a fixed rate fixes its whole group, on every toy network in both uptake modes."""
    for name, (R, S) in TOYS.items():
        crn, cx = network(R, S)
        fb = _pairs(_flux_balance(crn))
        for detach in (True, False):
            it = kimonu.kinetic_flux_coupling(crn, detach_uptakes=detach)
            tan = kimonu.kinetic_flux_coupling(crn, detach_uptakes=detach, route="tangent")
            assert fb <= _pairs(it.groups) <= _pairs(tan.groups), (name, detach)
            assert set(it.fixed_rate) <= set(tan.fixed_rate), (name, detach)
            for found in (it, tan):
                assert sorted(j for g in found.groups for j in g) == list(range(len(R)))
                for g in found.groups:
                    assert not (g & set(found.fixed_rate)) or g <= set(found.fixed_rate), (name, detach, g)


def test_groups_lift_the_kinetic_modules():
    """Two reactions whose substrate complexes share a kinetic module share a group, and two reactions that
    consume complexes (not detached sources) share a group only when their complexes share a module."""
    for name, (R, S) in TOYS.items():
        crn, cx = network(R, S)
        for detach in (True, False):
            for route, modules in (("iterated", kimonu.iterated_route(crn, detach_uptakes=detach).modules),
                                   ("tangent", kimonu.tangent_route(crn, detach_uptakes=detach).modules)):
                found = kimonu.kinetic_flux_coupling(crn, detach_uptakes=detach, route=route)
                module_of = {c: k for k, g in enumerate(modules) for c in g}
                substrate = {j: cx.index(s) for j, (s, _) in enumerate(R) if not (detach and s == "0")}
                for i in substrate:
                    for j in substrate:
                        same_module = module_of[substrate[i]] == module_of[substrate[j]]
                        assert same_module == (found.group_of[i] == found.group_of[j]), (name, detach, route, i, j)


def test_kinase_two_steady_states():
    """At k = K the attached kinase has the steady states X = 10 with (C, P) = (1, 4) and (4, 1). The
    iterated route couples r3 and r5 and fixes the uptake, and both hold at both states. The tangent route
    calls every rate fixed, which holds near each state but not across them: r3, r4 and r5 differ. The true
    fixed rates are r0, r1 and r2, since X has ACR."""
    crn, cx = network(KINASE, KINASE_SPECIES)
    r1, r2 = (_rates({"X": 10, "C": c, "P": p}) for c, p in ((1, 4), (4, 1)))
    it = kimonu.kinetic_flux_coupling(crn, detach_uptakes=False)
    assert [sorted(g) for g in it.groups if len(g) > 1] == [[3, 5]] and it.fixed_rate == [0]
    assert r1[3] / r1[5] == r2[3] / r2[5]
    assert all(r1[j] == r2[j] for j in it.fixed_rate)
    tan = kimonu.kinetic_flux_coupling(crn, detach_uptakes=False, route="tangent")
    assert tan.fixed_rate == list(range(len(KINASE)))
    assert [j for j in range(len(KINASE)) if r1[j] == r2[j]] == [0, 1, 2]


def test_ratio_is_the_ratio_at_the_reference_flux():
    """For the rate constants k = r*, the open network detached has the steady states x = t·1 with source
    u = t, and every rate is t·r*_j, so coupled rates keep the ratio r*_i / r*_j."""
    crn, cx = network(OPEN, ["a", "b"])
    found = kimonu.kinetic_flux_coupling(crn)
    r = found.flux
    for t in (Fraction(1, 3), Fraction(7, 2)):
        x = {"a": t, "b": t}
        rates = [r[0] * t] + _rates(x, k=r[1:], reactions=OPEN[1:])
        dx = {"a": rates[0] - rates[1] + rates[2], "b": rates[1] - rates[2] - rates[3]}
        assert not any(dx.values())
        for g in found.groups:
            i, *rest = sorted(g)
            assert all(rates[j] / rates[i] == r[j] / r[i] for j in rest)


def test_mass_balance_coupling_gives_the_flux_ratios():
    """EnvZ-OmpR: XD <-> X carry the same flux (λ = 1), and so do XT -> Xp and XpY -> X + Yp. On the open network,
    the uptake and the export carry the same flux, and a -> b and b -> a are coupled to nothing."""
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    found = kimonu.mass_balance_flux_coupling(crn)
    assert [sorted(g) for g in found.groups if len(g) > 1] == [[0, 1], [4, 7]] and found.blocked == []
    assert found.ratio[1] == found.ratio[7] == 1 and found.n_groups == 2 and found.largest_group == 2
    assert [sorted(g) for g in found.groups] == sorted(sorted(g) for g in _flux_balance(crn))
    crn, cx = network(OPEN, ["a", "b"])
    found = kimonu.mass_balance_flux_coupling(crn)
    assert found.groups == [{0, 3}, {1}, {2}] and list(found.group_of) == [1, 2, 3, 1]
    flux = kimonu.interior_flux(crn)
    for g in found.groups:
        i = min(g)
        assert all(flux[j] == found.ratio[j] * flux[i] for j in g)


def test_mass_balance_coupling_reports_blocked_reactions_and_ratios_other_than_one():
    """a -> 2b with b -> 0 and a -> c without an exit: the uptake and a -> 2b are coupled, the export of b runs at
    twice their rate, and a -> c is blocked, a zero row of K."""
    crn, cx = network([("0", "a"), ("a", "b+b"), ("b", "0"), ("a", "c")], ["a", "b", "c"])
    found = kimonu.mass_balance_flux_coupling(crn)
    assert found.blocked == [3] and found.groups == [{0, 1, 2}] and list(found.group_of) == [1, 1, 1, 0]
    assert found.ratio == {0: 1, 1: 1, 2: 2}


def test_unknown_route():
    crn, cx = network(OPEN, ["a", "b"])
    with pytest.raises(ValueError):
        kimonu.kinetic_flux_coupling(crn, route="sampling")
