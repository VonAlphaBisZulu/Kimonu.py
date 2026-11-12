"""build_crn on cobra models: the zero complex has its own row, so uptake and export reactions have
a substrate and a product complex. Without that row, mass action loses the uptake rates (W u = 0 is
then not the steady-state condition) and kinetic modules couple complexes that are not coupled."""
import os
from fractions import Fraction

import numpy as np
import pytest
import cobra

import kimonu
from kimonu.io import build_crn, exact_matrix, irreversible, substrate_complexes, zero_complexes
from kimonu.linalg import entries

SHARED = os.environ.get(
    "KIMONU_SHARED_MODELS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "shared_models"),
)


def _toy():
    """uptake of a, a <-> b, b -> export"""
    m = cobra.Model("toy")
    a, b = cobra.Metabolite("a"), cobra.Metabolite("b")
    for rid, st in [("up", {a: 1}), ("ab", {a: -1, b: 1}), ("ba", {b: -1, a: 1}), ("out", {b: -1})]:
        r = cobra.Reaction(rid, lower_bound=0, upper_bound=1000)
        m.add_reactions([r]); r.add_metabolites(st)
    return m


def test_zero_complex_has_its_own_row():
    crn = build_crn(_toy())
    A = crn.A.toarray()
    assert len(zero_complexes(crn.Y)) == 1
    assert all(sorted(A[:, j][A[:, j] != 0]) == [-1, 1] for j in range(A.shape[1]))
    assert np.array_equal((crn.Y @ crn.A).toarray(), crn.N_float.toarray())
    assert np.array_equal(crn.N.to_numpy().astype(float), crn.N_float.toarray())
    assert isinstance(crn, kimonu.CRN) and crn.pins == {}


def test_attached_uptake_fixes_its_product():
    """psi_0 = 1: with a constant inflow, a's steady state is fixed, as in the matrix-built test."""
    crn = build_crn(_toy())
    it = kimonu.iterated_route(crn.A, crn.Y, crn.N, detach_uptakes=False, metabolite_ids=crn.metabolite_ids)
    assert "a" in it.acr


def _e_coli_core():
    """e_coli_core with reversible reactions split and blocked reactions removed, as for the benchmarks"""
    import straindesign as sd
    path = os.path.join(SHARED, "e_coli_core.mat")
    if not os.path.exists(path):
        pytest.skip("shared model e_coli_core not available")
    cobra.Configuration().solver = "glpk"
    m = cobra.io.load_matlab_model(path)
    for r in list(m.reactions):
        if r.lower_bound < 0:
            rv = cobra.Reaction(r.id + "_rev", lower_bound=0, upper_bound=-r.lower_bound)
            m.add_reactions([rv]); rv.add_metabolites({k: -v for k, v in r.metabolites.items()})
            r.lower_bound = 0
    ranges = sd.fva(m, solver="glpk")
    m.remove_reactions([r for r in ranges.index if abs(ranges.at[r, "maximum"]) <= 1e-9], remove_orphans=True)
    return m


def test_e_coli_core_has_no_coupling_from_dropped_uptakes():
    """With the uptake rates missing from W, 2pg and 3pg came out as an ACRR pair; an exact second
    steady state with the same rate constants refutes it."""
    crn = build_crn(_e_coli_core())
    it = kimonu.iterated_route(crn.A, crn.Y, crn.N, detach_uptakes=True, metabolite_ids=crn.metabolite_ids)
    assert {"2pg_c", "3pg_c"} not in [set(p) for p in it.acrr]


def _bounded(**bounds):
    """x -> y with the given lower and upper bound, between an uptake of x and an export of y"""
    m = cobra.Model("bounded")
    x, y = cobra.Metabolite("x"), cobra.Metabolite("y")
    for rid, st, lb, ub in [("up", {x: 1}, 0, 1000), ("r", {x: -1, y: 1}, bounds["lb"], bounds["ub"]), ("out", {y: -1}, 0, 1000)]:
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        m.add_reactions([r])
        r.add_metabolites(st)
    return m


def test_a_network_without_reversible_reactions_is_left_as_it_is():
    m = _bounded(lb=0, ub=1000)
    assert irreversible(m) is m


def test_a_reversible_reaction_becomes_two():
    m = _bounded(lb=-5, ub=7)
    split = irreversible(m)
    assert m.reactions.r.lower_bound == -5 and len(m.reactions) == 3        # the input is not changed
    forward, reverse = split.reactions.r, split.reactions.r_rev
    assert (forward.lower_bound, forward.upper_bound, reverse.lower_bound, reverse.upper_bound) == (0, 7, 0, 5)
    assert {k.id: v for k, v in reverse.metabolites.items()} == {"x": 1, "y": -1}


def test_a_reaction_that_only_runs_backwards_is_replaced():
    split = irreversible(_bounded(lb=-5, ub=-2))
    assert "r" not in split.reactions and (split.reactions.r_rev.lower_bound, split.reactions.r_rev.upper_bound) == (2, 5)


def test_build_crn_gives_a_reversible_reaction_two_directions_between_the_same_complexes():
    crn = build_crn(_bounded(lb=-5, ub=7))
    assert crn.reaction_ids == ["up", "r", "out", "r_rev"]
    A = crn.A.toarray()
    assert (A[:, 1] == -A[:, 3]).all()
    assert crn.A.shape[0] == 3                                               # zero complex, x, y


def test_the_iterated_route_runs_on_a_model_with_reversible_reactions():
    """Written as one reaction, r would have to run forwards at every steady state of mass action."""
    m = _bounded(lb=-1000, ub=1000)
    crn = build_crn(m)
    it = kimonu.iterated_route(crn.A, crn.Y, crn.N, detach_uptakes=True)
    assert it.flux and all(x > 0 for x in it.flux)


def test_the_split_of_build_crn_is_the_split_by_hand():
    path = os.path.join(SHARED, "e_coli_core.mat")
    if not os.path.exists(path):
        pytest.skip("shared model e_coli_core not available")
    m = cobra.io.load_matlab_model(path)
    by_hand = _e_coli_core()
    crn = build_crn(m)
    assert set(by_hand.reactions.list_attr("id")) <= set(crn.reaction_ids)


def test_concordance_does_not_change_with_the_split():
    """Both directions together only add the kernel vector (1, 1) on them, whose row of A K is zero."""
    def sizes(r):
        return r.n_balanced, sorted(len(g) for g in r.modules)
    m = _bounded(lb=-1000, ub=1000)
    assert sizes(kimonu.concordance(m)) == sizes(kimonu.concordance(irreversible(m)))


def test_exact_matrix_keeps_the_declared_decimals():
    """-7.73333333333333e-08 has a denominator of 10**22, which int64 cannot clear; Fractions stay as they are."""
    M = exact_matrix([(0, 0, -7.73333333333333e-08), (0, 1, Fraction(1, 3))], (1, 2))
    assert sorted(entries(M)[0]) == [(0, 0, Fraction("-7.73333333333333e-08")), (0, 1, Fraction(1, 3))]


def test_substrate_complexes_of_the_toy():
    """up: 0 -> a, ab: a -> b, ba: b -> a, out: b -> 0, with the complexes 0, a, b in order of first use"""
    crn = build_crn(_toy())
    assert crn.complexes[0] == () and substrate_complexes(crn.A) == {0: 0, 1: 1, 2: 2, 3: 2}


def test_prepare_removes_blocked_reactions_and_leaves_the_model_alone():
    """A dead end y -> z (nothing feeds y) is blocked, and its metabolites go with it."""
    m = _toy()
    y, z = cobra.Metabolite("y"), cobra.Metabolite("z")
    dead = cobra.Reaction("dead", lower_bound=0, upper_bound=1000)
    m.add_reactions([dead]); dead.add_metabolites({y: -1, z: 1})
    crn = kimonu.prepare(m)
    assert crn.reaction_ids == ["up", "ab", "ba", "out"] and crn.metabolite_ids == ["a", "b"] and crn.pins == {}
    assert "dead" in m.reactions


def test_prepare_splits_reversible_reactions():
    crn = kimonu.prepare(_bounded(lb=-1000, ub=1000))
    assert crn.reaction_ids == ["up", "r", "out", "r_rev"]
