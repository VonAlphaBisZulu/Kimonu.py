"""Phase-1 concordance core tests.

Validates the *method* (N = Y·A, K = ker(N), M = A·K, proportional-row grouping) structurally and
for internal consistency. Exact agreement with the Julia PoC counts (EnvZ 13/5/3/8, e_coli
895/209/110/2436) additionally requires the elementary-step splitting preprocessing the PoC used —
see ROADMAP (Phase 1b); the shared SBML here are the un-split models.
"""
import os
from fractions import Fraction

import numpy as np
import pytest
from cobra.io import read_sbml_model

import kimonu
from kimonu.io import build_crn, irreversible

SHARED = os.environ.get(
    "KIMONU_SHARED_MODELS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "shared_models"),
)


def _model(name):
    path = os.path.join(SHARED, name + ".xml")
    if not os.path.exists(path):
        pytest.skip(f"shared model {name} not available")
    return read_sbml_model(path)


def test_N_equals_Y_times_A():
    """The complex factorization must reproduce the stoichiometric matrix: N == Y·A."""
    from cobra.util.array import create_stoichiometric_matrix
    m = _model("EnvZ_OmpR")
    crn = build_crn(m)
    N_ref = create_stoichiometric_matrix(irreversible(m))          # metabolites x reactions
    N_ya = np.asarray((crn.Y @ crn.A).todense(), dtype=float)
    assert np.array_equal(N_ref, N_ya)


def test_fractional_coefficients_are_exact():
    """Fractional stoichiometry must survive into N exactly.

    e_coli_core has 0.5 O2 in CYTBD and a biomass reaction full of decimals. Rounding those to
    integers still produces a plausible-looking network with the same complexes, so the damage is
    invisible in the counts and only shows up as a wrong nullspace.
    """
    from cobra.util.array import create_stoichiometric_matrix
    m = _model("e_coli_core")
    crn = build_crn(m)
    S = create_stoichiometric_matrix(irreversible(m))

    assert np.array_equal(S, np.asarray((crn.Y @ crn.A).todense(), dtype=float))
    assert any(float(c) != int(c) for comp in crn.complexes for _, c in comp), \
        "test model has no fractional coefficients — it cannot catch the rounding it is here to catch"

    # N_exact holds the coefficients themselves, so it must reproduce S entry for entry.
    assert np.allclose(np.asarray(crn.N_exact.to_numpy(), dtype=float), S)


def test_partition_is_internally_consistent():
    """Balanced complexes have zero activity; complexes in a module have proportional M-rows."""
    import straindesign as sd
    m = _model("EnvZ_OmpR")
    crn = build_crn(m)
    K = sd.sparse_nullspace(crn.N)
    assert not isinstance(K, sd.ExactCOO)         # small model: fast int64 path
    M = np.asarray((crn.A @ K).todense(), dtype=object)

    r = kimonu.concordance(m, preprocess=False)      # partition of this exact crn
    for i in r.balanced:
        assert all(v == 0 for v in M[i]), "balanced complex has nonzero activity"
    for grp in r.modules:
        if len(grp) < 2:
            continue
        base = M[grp[0]]
        b0 = next(x for x in base if x != 0)
        ref = [Fraction(int(x), int(b0)) for x in base]
        for i in grp[1:]:
            row = M[i]
            r0 = next(x for x in row if x != 0)
            assert [Fraction(int(x), int(r0)) for x in row] == ref, "module rows not proportional"


def test_deterministic_and_runs_on_ecoli():
    m = _model("e_coli_core")
    r1 = kimonu.concordance(m)
    r2 = kimonu.concordance(m)
    assert (r1.n_complexes, r1.n_balanced, r1.n_modules, r1.n_concordant_pairs) == \
           (r2.n_complexes, r2.n_balanced, r2.n_modules, r2.n_concordant_pairs)
    # balanced + partitioned complexes account for every complex
    assert r1.n_balanced + sum(len(g) for g in r1.modules) == r1.n_complexes


def _toy(reactions, pins=None):
    """A cobra model from {id: {metabolite: coefficient}}; reactions in ``pins`` get lb = ub = value."""
    import cobra
    m = cobra.Model("toy")
    m.add_metabolites([cobra.Metabolite(s) for s in sorted({s for st in reactions.values() for s in st})])
    for rid, st in reactions.items():
        v = (pins or {}).get(rid)
        r = cobra.Reaction(rid, lower_bound=v if v else 0, upper_bound=v if v else 1000)
        m.add_reactions([r])
        r.add_metabolites({m.metabolites.get_by_id(s): c for s, c in st.items()})
    return m


# r_i = r_k + r_x and r_j = r_l + r_x: rows i and j of ker(N) are not proportional
TWO_FEEDS = {"k": {"U": 1}, "l": {"W": 1}, "x": {"U": 1, "W": 1}, "i": {"U": -1}, "j": {"W": -1}}


def _coupled(model):
    """Reaction ids grouped by proportional rows of the kernel that concordance() computes."""
    import straindesign as sd
    from kimonu.concordance import _pins, _with_pins
    from kimonu.linalg import null_space_rows, proportional_classes
    crn = build_crn(model)
    rows = null_space_rows(_with_pins(crn.N_exact, crn.reaction_ids, _pins(model, sd.fva(model))))
    return {frozenset(crn.reaction_ids[i] for i in c) for c in proportional_classes(rows)[1] if len(c) > 1}


def test_two_equal_pins_couple_reactions_that_still_vary():
    """With r_k = r_l = 5, r_i = r_j = 5 + r_x for every r_x: a coupling ker(N) alone does not show."""
    assert _coupled(_toy(TWO_FEEDS)) == set()
    assert _coupled(_toy(TWO_FEEDS, {"k": 5, "l": 5})) == {frozenset("kl"), frozenset("ij")}


def test_one_pin_or_mismatched_pins_add_no_coupling():
    """One pin leaves the span of the steady states at ker(N); pins at 5 and 3 only fix k against l."""
    assert _coupled(_toy(TWO_FEEDS, {"k": 5})) == set()
    assert _coupled(_toy(TWO_FEEDS, {"k": 5, "l": 3})) == {frozenset("kl")}


def test_pins_reach_the_concordance_partition():
    """Two feed cycles through the zero complex: the zero complex is a module of its own, and pinning
    both feeds makes every complex concordant."""
    net = {"k": {"X": 1}, "y": {"Y": 1}, "m": {"X": -1, "Y": -1},
           "l": {"W": 1}, "v": {"V": 1}, "n": {"W": -1, "V": -1}}
    sizes = lambda r: sorted(len(g) for g in r.modules)
    assert sizes(kimonu.concordance(_toy(net))) == [1, 3, 3]
    assert sizes(kimonu.concordance(_toy(net, {"k": 2}))) == [1, 3, 3]
    assert sizes(kimonu.concordance(_toy(net, {"k": 2, "l": 7}))) == [7]


def _random_irreversible(rng):
    """A small random network with bounds 0 to 1000; dead ends and one-way loops make some reactions blocked."""
    import cobra
    m = cobra.Model("random")
    n_met, n_rxn = rng.integers(3, 12), rng.integers(4, 30)
    mets = [cobra.Metabolite(f"M{i}") for i in range(n_met)]
    m.add_metabolites(mets)
    for j in range(n_rxn):
        r = cobra.Reaction(f"R{j}", lower_bound=0, upper_bound=1000)
        m.add_reactions([r])
        picked = rng.choice(n_met, size=rng.integers(1, 4), replace=False)
        r.add_metabolites({mets[i]: float(rng.choice([-2, -1, 1, 2])) for i in picked})
    return m


def test_one_lp_finds_the_blocked_reactions_of_fva():
    import cobra
    import straindesign as sd
    from kimonu.concordance import _preprocess
    cobra.Configuration().solver = "glpk"
    rng = np.random.default_rng(0)
    seen = 0
    for _ in range(60):
        m = _random_irreversible(rng)
        ranges = sd.fva(m, solver="glpk")
        fva = {r for r in ranges.index if abs(ranges.at[r, "maximum"]) <= 1e-9}
        blocked, pins = _preprocess(m)
        assert set(blocked) == fva and pins == {}
        seen += 0 < len(fva) < len(m.reactions)
    assert seen >= 10


def test_one_lp_keeps_both_branches_that_compete_for_flux():
    """r_u = r_a + r_b: maximizing the sum of the rates under bounds can stop with r_b = 0."""
    import cobra
    from kimonu.concordance import _preprocess
    m = _toy({"u": {"X": 1}, "a": {"X": -1}, "b": {"X": -1}, "dead": {"Y": 1}})
    assert _preprocess(m)[0] == ["dead"]


def test_one_lp_respects_closed_reactions():
    from kimonu.concordance import _preprocess
    m = _toy({"u": {"X": 1}, "a": {"X": -1}})
    m.reactions.get_by_id("a").upper_bound = 0
    assert set(_preprocess(m)[0]) == {"u", "a"}


def test_forced_flux_goes_through_fva_and_finds_pins():
    from kimonu.concordance import _preprocess
    blocked, pins = _preprocess(_toy(TWO_FEEDS, {"k": 5, "l": 5}))
    assert blocked == [] and pins == {"k": Fraction(5), "l": Fraction(5)}
