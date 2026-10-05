"""build_crn on cobra models: the zero complex has its own row, so uptake and export reactions have
a substrate and a product complex. Without that row, mass action loses the uptake rates (W u = 0 is
then not the steady-state condition) and kinetic modules couple complexes that are not coupled."""
import os

import numpy as np
import pytest
import cobra

import kimonu
from kimonu.io import build_crn, zero_complexes

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
    assert np.array_equal((crn.Y @ crn.A).toarray(), crn.N.toarray())


def test_attached_uptake_fixes_its_product():
    """psi_0 = 1: with a constant inflow, a's steady state is fixed, as in the matrix-built test."""
    crn = build_crn(_toy())
    kin = kimonu.kinetic_modules(crn.A, crn.Y, crn.N_exact)
    acr = kimonu.acr_acrr(crn.Y, kin.modules, metabolite_ids=crn.metabolite_ids, known_acr=kin.known_acr)
    assert "a" in acr.acr


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
    kin = kimonu.kinetic_modules(crn.A, crn.Y, crn.N_exact, detach_uptakes=True)
    acr = kimonu.acr_acrr(crn.Y, kin.modules, metabolite_ids=crn.metabolite_ids, known_acr=kin.known_acr)
    assert {"2pg_c", "3pg_c"} not in [set(p) for p in acr.acrr]
