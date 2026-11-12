"""Networks and preprocessing: build_crn on cobra models, preprocess, as_crn, and the blocked-reaction LP.

The zero complex has its own row, so uptake and export reactions have a substrate and a product complex. Without that
row, mass action loses the uptake rates (W u = 0 is then not the steady-state condition) and kinetic modules couple
complexes that are not coupled."""
import os
from fractions import Fraction

import numpy as np
import pytest
import cobra
from scipy import sparse

import kimonu
from kimonu.linalg import matrix_product
from kimonu.network import _blocked_reactions, build_crn, make_irreversible, substrate_complexes, zero_complexes

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


def _dense_Y(crn):
    Y = np.zeros((len(crn.metabolite_ids), len(crn.Yt)))
    for c, y in enumerate(crn.Yt):
        for i, v in y:
            Y[i, c] = v
    return Y


def test_zero_complex_has_its_own_row():
    crn = build_crn(_toy())
    A = crn.A.toarray()
    assert zero_complexes(crn) == {0} and crn.Yt[0] == ()
    assert all(sorted(A[:, j][A[:, j] != 0]) == [-1, 1] for j in range(A.shape[1]))
    assert np.array_equal(_dense_Y(crn) @ A, crn.N.to_numpy().astype(float))
    assert crn.Yt == [(), ((0, 1),), ((1, 1),)] and crn.complexes == [(), (("a", 1),), (("b", 1),)]
    assert isinstance(crn, kimonu.CRN) and crn.pins == {}


def test_crn_from_matrices_is_the_crn_of_the_model():
    """The same network given as A and Y gives the CRN that build_crn gives for the cobra model."""
    crn = build_crn(_toy())
    again = kimonu.crn_from_matrices(crn.A, _dense_Y(crn), metabolite_ids=crn.metabolite_ids,
                                     reaction_ids=crn.reaction_ids)
    assert (again.A != crn.A).nnz == 0 and again.Yt == crn.Yt and again.complexes == crn.complexes
    assert np.array_equal(again.N.to_numpy(), crn.N.to_numpy()) and again.pins == {}
    unnamed = kimonu.crn_from_matrices(crn.A.toarray(), _dense_Y(crn))
    assert unnamed.metabolite_ids == [0, 1] and unnamed.reaction_ids == [0, 1, 2, 3]
    assert unnamed.complexes == [(), ((0, 1),), ((1, 1),)]


def test_compositions_are_exact_from_floats():
    """A float coefficient enters Yt as its decimal, and N = Y·A keeps it."""
    crn = kimonu.crn_from_matrices(np.array([[-1], [1]]), np.array([[0.1, 0.0], [0.0, 1.496]]))
    assert crn.Yt == [((0, Fraction(1, 10)),), ((1, Fraction(1496, 1000)),)]
    assert matrix_product(crn.A.T, crn.Yt) == [((0, Fraction(-1, 10)), (1, Fraction(1496, 1000)))]


def test_attached_uptake_fixes_its_product():
    """psi_0 = 1: with a constant inflow, a's steady state is fixed, as in the matrix-built test."""
    crn = build_crn(_toy())
    it = kimonu.iterated_route(crn, detach_uptakes=False)
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
    it = kimonu.iterated_route(crn, detach_uptakes=True)
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
    assert make_irreversible(m) is m


def test_a_reversible_reaction_becomes_two():
    m = _bounded(lb=-5, ub=7)
    split = make_irreversible(m)
    assert m.reactions.r.lower_bound == -5 and len(m.reactions) == 3        # the input is not changed
    forward, reverse = split.reactions.r, split.reactions.r_rev
    assert (forward.lower_bound, forward.upper_bound, reverse.lower_bound, reverse.upper_bound) == (0, 7, 0, 5)
    assert {k.id: v for k, v in reverse.metabolites.items()} == {"x": 1, "y": -1}


def test_a_reaction_that_only_runs_backwards_is_replaced():
    split = make_irreversible(_bounded(lb=-5, ub=-2))
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
    it = kimonu.iterated_route(crn, detach_uptakes=True)
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
    assert sizes(kimonu.concordance(m)) == sizes(kimonu.concordance(make_irreversible(m)))


def test_substrate_complexes_of_the_toy():
    """up: 0 -> a, ab: a -> b, ba: b -> a, out: b -> 0, with the complexes 0, a, b in order of first use"""
    crn = build_crn(_toy())
    assert crn.complexes[0] == () and substrate_complexes(crn) == {0: 0, 1: 1, 2: 2, 3: 2}


def test_preprocess_removes_blocked_reactions_and_leaves_the_model_alone():
    """A dead end y -> z (nothing feeds y) is blocked, and its metabolites go with it."""
    m = _toy()
    y, z = cobra.Metabolite("y"), cobra.Metabolite("z")
    dead = cobra.Reaction("dead", lower_bound=0, upper_bound=1000)
    m.add_reactions([dead]); dead.add_metabolites({y: -1, z: 1})
    crn = kimonu.preprocess(m)
    assert crn.reaction_ids == ["up", "ab", "ba", "out"] and crn.metabolite_ids == ["a", "b"] and crn.pins == {}
    assert "dead" in m.reactions


def test_preprocess_splits_reversible_reactions():
    crn = kimonu.preprocess(_bounded(lb=-1000, ub=1000))
    assert crn.reaction_ids == ["up", "r", "out", "r_rev"]


def test_as_crn_passes_a_crn_and_preprocesses_a_model():
    """A CRN comes back as it is. A cobra model is preprocessed, so its dead end goes, or only built with
    preprocess=False, so the dead end stays."""
    m = _toy()
    y, z = cobra.Metabolite("y"), cobra.Metabolite("z")
    dead = cobra.Reaction("dead", lower_bound=0, upper_bound=1000)
    m.add_reactions([dead]); dead.add_metabolites({y: -1, z: 1})
    crn = kimonu.preprocess(m)
    assert kimonu.as_crn(crn) is crn
    assert kimonu.as_crn(m).reaction_ids == crn.reaction_ids == ["up", "ab", "ba", "out"]
    assert kimonu.as_crn(m, preprocess=False).reaction_ids == ["up", "ab", "ba", "out", "dead"]


def test_every_analysis_takes_a_cobra_model():
    """Each entry point preprocesses a cobra model as as_crn does, and gives what it gives on the CRN."""
    m = _bounded(lb=-1000, ub=1000)
    crn = kimonu.preprocess(m)
    flux = kimonu.interior_flux(crn)
    assert kimonu.interior_flux(m) == flux
    assert kimonu.concordance(m).modules == kimonu.concordance(crn).modules
    for route in (kimonu.iterated_route, kimonu.tangent_route):
        assert route(m).modules == route(crn).modules and route(m).acrr == route(crn).acrr
    assert kimonu.reference_route(m).modules == kimonu.reference_route(crn).modules
    assert kimonu.certify(m).certified == kimonu.certify(crn).certified
    assert kimonu.kinetic_flux_coupling(m).groups == kimonu.kinetic_flux_coupling(crn).groups
    assert kimonu.mass_balance_flux_coupling(m).groups == kimonu.mass_balance_flux_coupling(crn).groups
    assert kimonu.rate_constants(m, flux, [2, 3]) == kimonu.rate_constants(crn, flux, [2, 3])
    assert kimonu.relation_values(m, flux).acrr == kimonu.relation_values(crn, flux).acrr
    assert kimonu.binomial_formulas(m).constraints == kimonu.binomial_formulas(crn).constraints
    assert kimonu.steady_state_map(m, flux) == kimonu.steady_state_map(crn, flux)
    assert kimonu.tangent_space(m, flux) == kimonu.tangent_space(crn, flux)
    B = kimonu.free_directions(crn, kimonu.iterated_route(crn).modules)
    assert kimonu.free_directions(m, kimonu.iterated_route(crn).modules) == B
    assert kimonu.readout(m, B) == kimonu.readout(crn, B)


def test_blocked_reactions_by_one_lp():
    """A chain runs only while its export is open, and a product that nothing consumes blocks its producers."""
    S = sparse.csr_matrix(np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]]))    # up a, a -> b, b -> c, c out
    assert _blocked_reactions(S).tolist() == [False] * 4
    assert _blocked_reactions(S, closed=[3]).tolist() == [True] * 4
    dead_end = sparse.csr_matrix(np.array([[1, -1, 0], [0, 1, 0], [0, 0, 1]]))           # b and the third are dead
    assert _blocked_reactions(dead_end).tolist() == [True, True, True]
