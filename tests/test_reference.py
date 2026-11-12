"""Reference-equivalent mode: upstream sets, the overlap merge and the Shinar-Feinberg rule.

"The reference" is the implementation of Langary et al. (2025). The two phases of the upstream set are small but easy
to get subtly wrong, and the failure mode is quiet. A module comes out one complex short, and the loss only shows up
as missing ACR or ACRR several steps later. The networks here are hand-built so the right answer is readable off the
graph, plus EnvZ-OmpR and the network of Figure 1, on which the mode stops at the coupled sets of the first round
of the iterated route and misses the ACR that the feedback finds.
"""
import numpy as np

import kimonu
from kimonu.reference import (complex_graph, linkage_classes, shinar_feinberg, strong_components, upstream_modules,
                              upstream_set)
from networks import ENVZ, ENVZ_SPECIES, FIG1, TOYS, network


def _envz():
    return network(ENVZ, ENVZ_SPECIES)


def _named(modules, cx):
    return sorted(sorted(cx[i] for i in m) for m in modules if len(m) > 1)


def test_terminality_is_judged_on_the_full_graph():
    """A component that feeds something outside the candidate set is still feeding.

    Chain a -> b -> c. Ask for the upstream set of {a, b} only. The component {b} has no successor within {a, b}, so
    judging terminality against the candidates would delete it and strand a. Against the full graph b still feeds c,
    so it survives. This is the distinction that decides whether whole modules keep their ACR and ACRR relations or
    fall apart into singletons.
    """
    crn, cx = network([("a", "b"), ("b", "c")], ["a", "b", "c"])
    out = upstream_set({cx.index("a"), cx.index("b")}, complex_graph(crn))
    assert cx.index("b") in out, "b feeds c in the full graph, so it is not terminal"


def test_phase_one_removes_entry_complexes_to_a_fixed_point():
    """Removal cascades: dropping one complex can turn its consumer into an entry complex.

    x -> a -> b -> c -> d with x outside the candidates. a is an entry complex (fed from x), and once it goes b is fed
    from outside too, and so on. A single pass would keep b and c.
    """
    crn, cx = network([("x", "a"), ("a", "b"), ("b", "c"), ("c", "d")], ["x", "a", "b", "c", "d"])
    out = upstream_set({cx.index(n) for n in ("a", "b", "c", "d")}, complex_graph(crn))
    for name in ("a", "b", "c"):
        assert cx.index(name) not in out, f"{name} should be removed by the phase 1 cascade"


def test_upstream_set_is_idempotent():
    """The upstream set of an upstream set is itself, so it is a closure."""
    crn, cx = _envz()
    graph = complex_graph(crn)
    res = kimonu.concordance(crn)
    for module in res.modules:
        once = upstream_set(set(module) | set(res.balanced), graph)
        if once:
            assert upstream_set(once | set(res.balanced), graph) - set(res.balanced) == once - set(res.balanced)


def test_complex_graph_of_envz():
    crn, cx = _envz()
    graph = complex_graph(crn)
    assert graph.succ[cx.index("XD")] == {cx.index("X")} and graph.pred[cx.index("XD")] == {cx.index("X")}
    assert graph.substrates[4] == [cx.index("XT")] and graph.producers[cx.index("Xp")] == [4]


def test_strong_components_partition_the_scope():
    """Every complex in scope lands in exactly one strong component, and XD <-> X <-> XT form one."""
    crn, cx = _envz()
    scope = set(range(len(cx)))
    comps = strong_components(complex_graph(crn).succ, scope)
    assert set().union(*comps) == scope
    assert sum(len(c) for c in comps) == len(scope), "components must be disjoint"
    assert {cx.index("XD"), cx.index("X"), cx.index("XT")} in comps


def test_linkage_classes_of_envz():
    """Four linkage classes, the EnvZ cycle and one for each of the three complexes with Y or Yp."""
    crn, cx = _envz()
    graph = complex_graph(crn)
    classes = linkage_classes(graph.succ, graph.pred, set(range(len(cx))))
    assert sorted(sorted(cx[c] for c in g) for g in classes) == [
        ["X", "XD", "XT", "Xp"], ["X+Yp", "Xp+Y", "XpY"], ["XD+Y", "XD+Yp", "XDYp"], ["XT+Y", "XT+Yp", "XTYp"]]


def test_envz_ompr_end_to_end():
    """Five balanced complexes, a seven-dimensional kernel, and three kinetic modules from upstream sets merged by
    overlap. No two coupled complexes differ in one metabolite, so the Shinar-Feinberg rule finds no ACR, where the
    iterated route proves the ACR of Yp."""
    crn, cx = _envz()
    ref = kimonu.reference_route(crn)
    assert ref.route == "reference" and ref.B is None and ref.flux is None and ref.dimension is None
    assert ref.concordance.n_balanced == 5 and ref.concordance.nullity == 7
    assert _named(ref.modules, cx) == [["X", "XD", "XT", "Xp+Y", "XpY"], ["XD+Yp", "XDYp"], ["XT+Yp", "XTYp"]]
    assert ref.acr == [] and ref.acrr == [] and ref.rounds == 2


def test_figure1_stops_at_the_coupled_sets_of_ker_w():
    """On the network of Figure 1 the mode gives the three coupled sets of ker(W) of the first round of the iterated
    route."""
    crn, cx = network(FIG1, list("ABCDEFGHI"))
    ref = kimonu.reference_route(crn)
    assert _named(ref.modules, cx) == [["A", "B", "C", "D+E", "F"], ["A+G", "I"], ["C+G", "H"]]


def test_upstream_modules_seed_with_the_given_modules_only():
    """With no seed, every complex is a singleton, except linkage classes made up entirely of balanced complexes."""
    crn, cx = _envz()
    res = kimonu.concordance(crn)
    modules, rounds = upstream_modules(crn, [], res.balanced)
    assert rounds == 1 and sorted(len(m) for m in modules) == [1] * len(cx)
    seeded, _ = upstream_modules(crn, [g for g in res.modules if len(g) > 1], res.balanced)
    assert seeded == kimonu.reference_route(crn).modules


def test_shinar_feinberg_rule():
    """A and A+G differ in G alone, so G has ACR. A+X and A+Y differ in X and Y with A in common, so X, Y have ACRR.
    A and B differ in two metabolites with nothing in common, which gives nothing."""
    species = ["A", "B", "G", "X", "Y"]
    Y = np.zeros((5, 5), dtype=np.int64)
    for c, comp in enumerate([{"A": 1}, {"A": 1, "G": 1}, {"A": 1, "X": 1}, {"A": 1, "Y": 1}, {"B": 1}]):
        for s, v in comp.items():
            Y[species.index(s), c] = v
    crn = kimonu.crn_from_matrices(np.zeros((5, 0), dtype=np.int64), Y)
    acr, acrr = shinar_feinberg(crn, [{0, 1}, {2, 3}, {0, 4}])
    assert acr == {2} and acrr == {(3, 4)}


def test_modules_are_a_partition_and_deterministic():
    """Kinetic modules cover every complex exactly once, singletons included, the same on every run."""
    for name, (R, S) in TOYS.items():
        crn, cx = network(R, S)
        runs = [kimonu.reference_route(crn).modules for _ in range(3)]
        assert sorted(c for m in runs[0] for c in m) == list(range(len(cx))), name
        assert runs[0] == runs[1] == runs[2], name


def test_iterated_route_contains_the_coupled_pairs_of_the_reference():
    """Every coupled complex pair of the reference-equivalent mode is also proven by the iterated route, in both uptake
    modes, with the zero complex left out."""
    for name, (R, S) in TOYS.items():
        crn, cx = network(R, S)
        zero = kimonu.network.zero_complexes(crn)
        pairs = lambda mods: {(a, b) for m in mods for a in m for b in m if a < b and a not in zero and b not in zero}
        ref = pairs(kimonu.reference_route(crn).modules)
        for detach in (True, False):
            assert ref <= pairs(kimonu.iterated_route(crn, detach_uptakes=detach).modules), (name, detach)
