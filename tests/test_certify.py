"""Both routes at one reference flux, and the routes on a cobra model.

EnvZ-OmpR has one steady state per rate-constant set up to its conservation laws, so the routes agree and the list is
certified. On the attached kinase of Figure 2 they differ, and X, C and P stay open, of which X is true.
"""
import warnings

import cobra
import pytest

import kimonu
from networks import ENVZ, ENVZ_SPECIES, KINASE, KINASE_SPECIES, network


def test_envz_is_certified_without_a_warning():
    crn, cx = network(ENVZ, ENVZ_SPECIES)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cert = kimonu.certify(crn)
    assert cert.certified and cert.open_pairs == [] and cert.open_acr == [] and cert.open_acrr == []
    assert cert.acr == ["Yp"] and len(cert.acrr) == 15 and cert.modules == cert.iterated.modules
    assert cert.iterated.flux is cert.tangent.flux is cert.flux


def test_kinase_gap_gives_a_warning_and_the_open_relations():
    """Attached, the iterated route proves the module {C+2P, P} and no ACR. The tangent route reports one module of
    all seven complexes and ACR for X, C and P."""
    crn, cx = network(KINASE, KINASE_SPECIES)
    with pytest.warns(UserWarning, match="stay open"):
        cert = kimonu.certify(crn, detach_uptakes=False)
    assert not cert.certified and cert.acr == [] and cert.open_acr == ["C", "P", "X"]
    n = len(cx)
    assert len(cert.open_pairs) == n * (n - 1) // 2 - 1
    assert (cx.index("C+P+P"), cx.index("P")) not in cert.open_pairs


def _toy():
    """0 -> a, a <-> b, b -> 0 and a dead end y -> z as a cobra model"""
    m = cobra.Model("toy")
    a, b, y, z = (cobra.Metabolite(s) for s in "abyz")
    for rid, st, lb in [("up", {a: 1}, 0), ("ab", {a: -1, b: 1}, -1000), ("out", {b: -1}, 0),
                        ("dead", {y: -1, z: 1}, 0)]:
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=1000)
        m.add_reactions([r]); r.add_metabolites(st)
    return m


def test_every_route_runs_on_a_cobra_model():
    """Detached, the supply is free and a, b only keep their ratio. Attached, both have ACR."""
    crn = kimonu.preprocess(_toy())
    assert crn.reaction_ids == ["up", "ab", "out", "ab_rev"]
    it = kimonu.iterated_route(_toy())
    assert it.route == "iterated" and it.acr == [] and it.acrr == [("a", "b")] and it.acrr_groups == [["a", "b"]]
    assert kimonu.iterated_route(crn, detach_uptakes=False).acr == ["a", "b"]
    assert kimonu.tangent_route(_toy()).acrr == [("a", "b")]
    assert kimonu.reference_route(_toy()).route == "reference"
    both = kimonu.certify(_toy())
    assert isinstance(both, kimonu.Certificate) and both.certified
