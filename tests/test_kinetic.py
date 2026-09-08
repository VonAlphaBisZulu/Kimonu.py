"""Kinetic-module tests: the upstream algorithm and what it feeds.

The two phases of Remark S2-1 are small but easy to get subtly wrong, and the failure mode is
quiet — a module comes out one complex short, and the loss only shows up as missing ACR/ACRR
several steps later. The networks here are hand-built so the right answer is readable off the
graph, plus EnvZ-OmpR as a literature-anchored end-to-end check.
"""
import numpy as np
import pytest
from scipy import sparse

import kimonu
from kimonu.kinetic import _Topology, _sccs, upstream_set, classical_deficiency


def _crn(reactions, species):
    """Build (A, Y, N) from a list of (substrate, product) complexes written as species multisets."""
    complexes = []
    for s, p in reactions:
        for c in (s, p):
            if c not in complexes:
                complexes.append(c)
    A = np.zeros((len(complexes), len(reactions)), dtype=np.int64)
    for j, (s, p) in enumerate(reactions):
        A[complexes.index(s), j] = -1
        A[complexes.index(p), j] = 1
    Y = np.zeros((len(species), len(complexes)), dtype=np.int64)
    for k, c in enumerate(complexes):
        for tok in c.split("+"):
            Y[species.index(tok), k] += 1
    return sparse.csr_matrix(A), sparse.csr_matrix(Y), sparse.csr_matrix(Y @ A), complexes


ENVZ = [("XD", "X"), ("X", "XD"), ("X", "XT"), ("XT", "X"), ("XT", "Xp"),
        ("Xp+Y", "XpY"), ("XpY", "Xp+Y"), ("XpY", "X+Yp"),
        ("XT+Yp", "XTYp"), ("XTYp", "XT+Yp"), ("XTYp", "XT+Y"),
        ("XD+Yp", "XDYp"), ("XDYp", "XD+Yp"), ("XDYp", "XD+Y")]
ENVZ_SPECIES = ["XD", "X", "XT", "Xp", "Y", "XpY", "Yp", "XTYp", "XDYp"]


def _envz():
    return _crn(ENVZ, ENVZ_SPECIES)


def test_terminality_is_judged_on_the_full_graph():
    """A component that feeds something outside the candidate set is still feeding.

    Chain a -> b -> c. Ask for the upstream set of {a, b} only. The component {b} has no successor
    *within* {a, b}, so judging terminality against the candidate would delete it and strand a.
    Against the full graph b still feeds c, so it survives. This is the distinction that decides
    whether whole modules keep their ACR/ACRR relations or fall apart into singletons.
    """
    A, Y, N, cx = _crn([("a", "b"), ("b", "c")], ["a", "b", "c"])
    topo = _Topology(A)
    out = upstream_set({cx.index("a"), cx.index("b")}, topo)
    assert cx.index("b") in out, "b feeds c in the full graph, so it is not terminal"


def test_phase_one_removes_entry_complexes_to_a_fixed_point():
    """Removal cascades: dropping one complex can turn its consumer into an entry complex.

    x -> a -> b -> c with x outside the candidate. `a` is an entry complex (fed from x), and once
    it goes `b` is fed from outside too, and so on. A single pass would keep b and c.
    """
    A, Y, N, cx = _crn([("x", "a"), ("a", "b"), ("b", "c"), ("c", "d")], ["x", "a", "b", "c", "d"])
    topo = _Topology(A)
    inside = {cx.index(n) for n in ("a", "b", "c", "d")}
    out = upstream_set(inside, topo)
    for name in ("a", "b", "c"):
        assert cx.index(name) not in out, f"{name} should be removed by the phase I cascade"


def test_upstream_set_is_idempotent():
    """The upstream set of an upstream set is itself — it is a closure."""
    A, Y, N, cx = _envz()
    topo = _Topology(A)
    res = kimonu.concordance_from_matrices(A, N)
    seed = set(res.modules[0]) | set(res.balanced)
    once = upstream_set(seed, topo)
    if once:
        assert upstream_set(once | set(res.balanced), topo) - set(res.balanced) == once - set(res.balanced)


def test_envz_ompr_end_to_end():
    """EnvZ-OmpR is the Shinar-Feinberg ACR network; the answer is known independently.

    Pinning the whole chain here means a regression anywhere between the kernel and the span test
    shows up as a failure on a network whose answer is in the literature.
    """
    A, Y, N, cx = _envz()
    res = kimonu.concordance_from_matrices(A, N)
    assert res.n_balanced == 5
    assert res.nullity == 7
    kin = kimonu.kinetic_modules(A, res.modules, res.balanced, Y=Y)
    assert kin.giant_module == 9, "the four terminal complexes are exactly what phase IV removes"
    acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=ENVZ_SPECIES)
    assert sorted(acr.acr) == ["Yp"], "phosphorylated OmpR is the ACR species"
    assert len(acr.acrr) == 15, "one 6-clique among X, XD, XT, XpY, XTYp, XDYp"


def test_deficiency_matches_the_definition():
    """delta = n - l - rank(N), computed independently of the implementation."""
    A, Y, N, cx = _envz()
    topo = _Topology(A)
    delta = classical_deficiency(Y, A, topo.succ, topo.pred)
    n = A.shape[0]
    rank = np.linalg.matrix_rank(np.asarray(N.todense(), dtype=float))
    linkage = len({frozenset(c) for c in _weakly_connected(topo, n)})
    assert delta == n - linkage - rank


def _weakly_connected(topo, n):
    seen, out = set(), []
    for start in range(n):
        if start in seen:
            continue
        comp, stack = set(), [start]
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u)
            stack.extend((topo.succ[u] | topo.pred[u]) - comp)
        seen |= comp
        out.append(comp)
    return out


def test_sccs_partition_the_scope():
    """Every complex in scope lands in exactly one strong component."""
    A, Y, N, cx = _envz()
    topo = _Topology(A)
    scope = set(range(A.shape[0]))
    comps = _sccs(topo.succ, scope)
    assert set().union(*comps) == scope
    assert sum(len(c) for c in comps) == len(scope), "components must be disjoint"


def test_modules_are_a_partition_of_the_complexes():
    """Kinetic modules cover every complex exactly once, singletons included."""
    A, Y, N, cx = _envz()
    res = kimonu.concordance_from_matrices(A, N)
    kin = kimonu.kinetic_modules(A, res.modules, res.balanced, Y=Y)
    union = set()
    for mod in kin.modules:
        assert not (union & mod), "modules must be disjoint"
        union |= mod
    assert union == set(range(A.shape[0]))


def test_acrr_is_transitively_closed():
    """ACRR is an equivalence relation, so the reported pairs must be a union of cliques.

    Kimonu gets this by construction — it groups identical rows rather than testing pairs — so a
    violation here means the grouping itself broke.
    """
    A, Y, N, cx = _envz()
    res = kimonu.concordance_from_matrices(A, N)
    kin = kimonu.kinetic_modules(A, res.modules, res.balanced, Y=Y)
    acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=ENVZ_SPECIES)
    adj = {}
    for a, b in acr.acrr:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    for v, nbrs in adj.items():
        for u in nbrs:
            assert adj[u] - {v} >= nbrs - {u}, f"{v} and {u} are ACRR but their neighbourhoods differ"


def test_acr_pairs_are_excluded_from_acrr():
    """A pair whose members both have ACR has a constant ratio trivially; it is not reported."""
    A, Y, N, cx = _envz()
    res = kimonu.concordance_from_matrices(A, N)
    kin = kimonu.kinetic_modules(A, res.modules, res.balanced, Y=Y)
    acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=ENVZ_SPECIES)
    have = set(acr.acr)
    assert not [p for p in acr.acrr if p[0] in have and p[1] in have]


def test_kinetic_modules_are_deterministic():
    """Same input, same output — no set-iteration order leaking into the result."""
    A, Y, N, cx = _envz()
    res = kimonu.concordance_from_matrices(A, N)
    runs = [kimonu.kinetic_modules(A, res.modules, res.balanced, Y=Y) for _ in range(3)]
    first = frozenset(frozenset(m) for m in runs[0].modules)
    for r in runs[1:]:
        assert frozenset(frozenset(m) for m in r.modules) == first


def test_without_Y_only_overlapping_sets_merge():
    """Y gates the coupling and deficiency merges, so omitting it can only give a finer partition."""
    A, Y, N, cx = _envz()
    res = kimonu.concordance_from_matrices(A, N)
    with_y = kimonu.kinetic_modules(A, res.modules, res.balanced, Y=Y)
    without = kimonu.kinetic_modules(A, res.modules, res.balanced, Y=None)
    assert without.giant_module <= with_y.giant_module
