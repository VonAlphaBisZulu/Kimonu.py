"""Iterated-route tests: the coupled sets of ker(W), the readout of ker((YΔ)ᵀ), and the feedback between them.

EnvZ-OmpR is the literature anchor. Yp has ACR for every choice of rate constants, and reaching it needs the feedback
step, since the coupled sets of ker(W) alone leave the complexes that differ by Yp in separate sets. The small open
network pins the attached and detached uptake modes.
"""
import kimonu
from kimonu.io import zero_complexes
from kimonu.iterated import coupled_columns
from networks import ENVZ, ENVZ_SPECIES, OPEN, matrices


def _partition(result):
    return frozenset(frozenset(m) for m in result.modules if len(m) > 1)


def test_envz_ompr_end_to_end():
    """Nine non-terminal complexes and X+Yp in one module, Yp the only ACR species, one 6-clique of ACRR."""
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    it = kimonu.iterated_route(A, Y, N, metabolite_ids=ENVZ_SPECIES)
    assert it.route == "iterated" and it.n_modules == 2 and it.giant_module == 10
    giant = next(m for m in it.modules if len(m) == 10)
    assert {cx[i] for i in giant} == {"XD", "X", "XT", "Xp+Y", "XpY", "XT+Yp", "XTYp", "XD+Yp", "XDYp", "X+Yp"}
    assert it.acr == ["Yp"]
    assert len(it.acrr) == 15
    assert it.rounds == 2 and it.dimension == 2 and len(it.B) == len(ENVZ_SPECIES)


def test_feedback_is_what_couples_the_yp_complexes():
    """Without the feedback the coupled sets of ker(W) leave XD and XD+Yp apart; their ratio is x_Yp."""
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    sets = coupled_columns(A, Y, N, kimonu.interior_flux(N), detach_uptakes=False)
    xd, xdyp = cx.index("XD"), cx.index("XD+Yp")
    assert not any(xd in s and xdyp in s for s in sets)


def test_coupled_columns_partition_the_columns_and_take_the_feedback_rows():
    """Every column of W lands in one class, and a row u_c = u_d puts c and d into one class."""
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    flux = kimonu.interior_flux(N)
    classes = coupled_columns(A, Y, N, flux)
    assert sorted(c for g in classes for c in g) == list(range(A.shape[0]))
    xd, xdyp = cx.index("XD"), cx.index("XD+Yp")
    joined = coupled_columns(A, Y, N, flux, equal=[(xd, xdyp)])
    assert any(xd in g and xdyp in g for g in joined)


def test_terminal_complexes_join_only_where_coupled():
    """X+Yp is coupled to the module through Yp, which has ACR (psi_{X+Yp}/psi_X = x_Yp). XD+Y and XT+Y are
    coupled to each other (ratio x_XD/x_XT), and Xp stays alone. Joining changes no ACR or ACRR."""
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    it = kimonu.iterated_route(A, Y, N, metabolite_ids=ENVZ_SPECIES)
    assert {frozenset(cx[i] for i in m) for m in it.modules if len(m) == 2} == {frozenset({"XD+Y", "XT+Y"})}
    assert it.module_of[cx.index("Xp")] not in {it.module_of[i] for m in it.modules if len(m) > 1 for i in m}
    assert it.acr == ["Yp"] and len(it.acrr) == 15


def test_modules_are_a_partition_of_the_complexes():
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    it = kimonu.iterated_route(A, Y, N)
    union = set()
    for mod in it.modules:
        assert not (union & mod), "modules must be disjoint"
        union |= mod
    assert union == set(range(A.shape[0]))


def test_result_does_not_depend_on_the_seed_or_the_draw():
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    first = kimonu.iterated_route(A, Y, N, seed=0)
    for seed in (1, 2, 3):
        for draw in ("uniform", "lp"):
            other = kimonu.iterated_route(A, Y, N, seed=seed, draw=draw)
            assert _partition(other) == _partition(first) and (other.acr, other.acrr) == (first.acr, first.acrr)


def test_a_given_flux_is_used():
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    flux = kimonu.interior_flux(N, seed=5, draw="lp")
    assert kimonu.iterated_route(A, Y, N, flux=flux).flux is flux


def test_attached_uptake_couples_through_the_zero_complex():
    """psi_0 = 1: with a constant inflow, a's steady state is fixed, so a has ACR."""
    A, Y, N, cx = matrices(OPEN, ["a", "b"])
    it = kimonu.iterated_route(A, Y, N, detach_uptakes=False, metabolite_ids=["a", "b"])
    zero, a = cx.index("0"), cx.index("a")
    assert it.module_of[zero] == it.module_of[a]
    assert "a" in it.acr


def test_detached_uptake_drops_the_zero_complex_coupling():
    """Each uptake with its own source: nothing is coupled to the zero complex any more, which is a module of its own,
    listed last."""
    A, Y, N, cx = matrices(OPEN, ["a", "b"])
    it = kimonu.iterated_route(A, Y, N, detach_uptakes=True, metabolite_ids=["a", "b"])
    zero = cx.index("0")
    assert not any(zero in m for m in it.modules if len(m) > 1)
    assert it.modules[-1] == {zero} and zero_complexes(Y) == {zero}
    assert "a" not in it.acr


def test_acrr_is_transitively_closed():
    """ACRR is an equivalence relation, so the reported pairs must be a union of cliques."""
    A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
    it = kimonu.iterated_route(A, Y, N, metabolite_ids=ENVZ_SPECIES)
    adj = {}
    for a, b in it.acrr:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    for v, nbrs in adj.items():
        for u in nbrs:
            assert adj[u] - {v} >= nbrs - {u}, f"{v} and {u} are ACRR but their neighbourhoods differ"
