"""Kinetic-module tests: the psi null space, the log-space closure, and the feedback between them.

EnvZ-OmpR is the literature anchor: Yp has ACR for every choice of rate constants, and reaching it
needs the feedback step, since the proportional-row readout alone leaves the complexes that differ
by Yp in separate sets. The small open networks pin the attached and detached uptake modes.
"""
from fractions import Fraction

import numpy as np
from scipy import sparse

import kimonu
from kimonu.kinetic import positive_flux, _psi_coupled


def _crn(reactions, species):
    """Build (A, Y, N) from a list of (substrate, product) complexes written as species multisets;
    "0" is the zero complex."""
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
        if c == "0":
            continue
        for tok in c.split("+"):
            Y[species.index(tok), k] += 1
    return sparse.csr_matrix(A), sparse.csr_matrix(Y), sparse.csr_matrix(Y @ A), complexes


ENVZ = [("XD", "X"), ("X", "XD"), ("X", "XT"), ("XT", "X"), ("XT", "Xp"),
        ("Xp+Y", "XpY"), ("XpY", "Xp+Y"), ("XpY", "X+Yp"),
        ("XT+Yp", "XTYp"), ("XTYp", "XT+Yp"), ("XTYp", "XT+Y"),
        ("XD+Yp", "XDYp"), ("XDYp", "XD+Yp"), ("XDYp", "XD+Y")]
ENVZ_SPECIES = ["XD", "X", "XT", "Xp", "Y", "XpY", "Yp", "XTYp", "XDYp"]

# uptake of a, a <-> b, b -> 0; a is fed from the zero complex at constant rate
OPEN = [("0", "a"), ("a", "b"), ("b", "a"), ("b", "0")]


def _partition(result):
    return frozenset(frozenset(m) for m in result.modules if len(m) > 1)


def test_positive_flux_is_an_exact_interior_steady_state():
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    flux = positive_flux(N)
    assert all(isinstance(x, Fraction) and x > 0 for x in flux)
    Nd = N.toarray()
    for i in range(Nd.shape[0]):
        assert sum(Fraction(int(Nd[i, j])) * flux[j] for j in range(Nd.shape[1])) == 0


def test_envz_ompr_end_to_end():
    """Nine non-terminal complexes and X+Yp in one module, Yp the only ACR species, one 6-clique of ACRR."""
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    kin = kimonu.kinetic_modules(A, Y, N)
    assert kin.n_modules == 2
    assert kin.giant_module == 10
    giant = next(m for m in kin.modules if len(m) == 10)
    assert {cx[i] for i in giant} == {"XD", "X", "XT", "Xp+Y", "XpY", "XT+Yp", "XTYp", "XD+Yp", "XDYp", "X+Yp"}
    acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=ENVZ_SPECIES, known_acr=kin.known_acr)
    assert sorted(acr.acr) == ["Yp"]
    assert len(acr.acrr) == 15


def test_feedback_is_what_couples_the_yp_complexes():
    """Without the feedback the psi readout leaves XD and XD+Yp apart; their ratio is x_Yp."""
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    sets, _ = _psi_coupled(A, Y, N, positive_flux(N), detach_uptakes=False)
    xd, xdyp = cx.index("XD"), cx.index("XD+Yp")
    assert not any(xd in s and xdyp in s for s in sets)


def test_terminal_complexes_join_only_where_coupled():
    """X+Yp is coupled to the module through Yp, which has ACR (psi_{X+Yp}/psi_X = x_Yp); XD+Y and XT+Y are
    coupled to each other (ratio x_XD/x_XT); Xp stays alone. Joining changes no ACR or ACRR."""
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    kin = kimonu.kinetic_modules(A, Y, N)
    giant = next(m for m in kin.modules if len(m) == 10)
    assert {cx[i] for i in giant} == {"XD", "X", "XT", "Xp+Y", "XpY", "XT+Yp", "XTYp", "XD+Yp", "XDYp", "X+Yp"}
    assert {frozenset(cx[i] for i in m) for m in kin.modules if len(m) == 2} == {frozenset({"XD+Y", "XT+Y"})}
    assert kin.module_of[cx.index("Xp")] not in {kin.module_of[i] for m in kin.modules if len(m) > 1 for i in m}
    acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=ENVZ_SPECIES, known_acr=kin.known_acr)
    assert sorted(acr.acr) == ["Yp"] and len(acr.acrr) == 15


def test_modules_are_a_partition_of_the_complexes():
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    kin = kimonu.kinetic_modules(A, Y, N)
    union = set()
    for mod in kin.modules:
        assert not (union & mod), "modules must be disjoint"
        union |= mod
    assert union == set(range(A.shape[0]))


def test_result_does_not_depend_on_the_seed():
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    first = _partition(kimonu.kinetic_modules(A, Y, N, seed=0))
    for seed in (1, 2, 3):
        assert _partition(kimonu.kinetic_modules(A, Y, N, seed=seed)) == first


def test_attached_uptake_couples_through_the_zero_complex():
    """psi_0 = 1: with a constant inflow, a's steady state is fixed, so a has ACR."""
    A, Y, N, cx = _crn(OPEN, ["a", "b"])
    kin = kimonu.kinetic_modules(A, Y, N)
    zero, a = cx.index("0"), cx.index("a")
    assert kin.module_of[zero] == kin.module_of[a]
    acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=["a", "b"], known_acr=kin.known_acr)
    assert "a" in acr.acr


def test_detached_uptake_drops_the_zero_complex_coupling():
    """Each uptake with its own source: nothing is coupled to the zero complex any more."""
    A, Y, N, cx = _crn(OPEN, ["a", "b"])
    kin = kimonu.kinetic_modules(A, Y, N, detach_uptakes=True)
    zero = cx.index("0")
    assert not any(zero in m for m in kin.modules if len(m) > 1)
    acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=["a", "b"], known_acr=kin.known_acr)
    assert "a" not in acr.acr


def test_acrr_is_transitively_closed():
    """ACRR is an equivalence relation, so the reported pairs must be a union of cliques."""
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    kin = kimonu.kinetic_modules(A, Y, N)
    acr = kimonu.acr_acrr(Y, kin.modules, metabolite_ids=ENVZ_SPECIES, known_acr=kin.known_acr)
    adj = {}
    for a, b in acr.acrr:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    for v, nbrs in adj.items():
        for u in nbrs:
            assert adj[u] - {v} >= nbrs - {u}, f"{v} and {u} are ACRR but their neighbourhoods differ"


def test_blocked_network_is_rejected():
    """a -> b with nothing consuming b has no positive steady state."""
    A, Y, N, cx = _crn([("a", "b")], ["a", "b"])
    try:
        positive_flux(N)
    except ValueError:
        return
    raise AssertionError("a network with blocked reactions must be rejected")
