"""Jacobian route against the iterated route on networks whose answer is known.

The one-species network 0 <-> A, 2A <-> 3A has three positive steady states for suitable rate
constants (x_A = 1, 2, 3 when dA/dt = -(x-1)(x-2)(x-3)). The tangent space at each is {0}, so the
Jacobian route calls A fixed; it is, along each branch, but at three different values. The iterated
route makes only global statements and reports no ACR. This is the case a multistationarity
certificate has to catch.
"""
import os, sys
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src")); sys.path.insert(0, os.path.join(os.path.dirname(B), "tests"))
import kimonu
from test_kinetic import _crn, ENVZ, ENVZ_SPECIES
from kinetic_variants.tangent_one_shot import one_shot_exact

CUBIC = [("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")]

def test_cubic_is_local_only():
    A, Y, N, cx = _crn(CUBIC, ["A"])
    _, acr, _, dim = one_shot_exact(A, Y, N, detach=False)
    assert dim == 0 and acr == [0], "each branch is an isolated point, so A is locally fixed"
    it = kimonu.kinetic_modules(A, Y, N)
    assert kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr).acr == [], "A takes three values"

def test_cubic_really_has_three_steady_states():
    roots = np.sort(np.roots([-1, 6, -11, 6]).real)        # 6 - 11x + 6x^2 - x^3
    assert np.allclose(roots, [1, 2, 3])

def test_envz_both_routes_agree():
    A, Y, N, cx = _crn(ENVZ, ENVZ_SPECIES)
    mods, acr, groups, dim = one_shot_exact(A, Y, N, detach=False, want_groups=True)
    it = kimonu.kinetic_modules(A, Y, N)
    assert {frozenset(m) for m in mods if len(m) > 1} == {frozenset(m) for m in it.modules if len(m) > 1}
    assert [ENVZ_SPECIES[i] for i in acr] == ["Yp"] and dim == 2
    assert sorted(len(g) for g in groups) == [6]


def test_descartes_step_on_known_networks():
    """silent on the cubic (three sign changes), fires on constant supply with first- and second-order loss,
    and treats detached uptakes as free, not constant"""
    from kinetic_variants.descartes import kinetic_modules_descartes
    A, Y, N, cx = _crn(CUBIC, ["A"])
    assert kinetic_modules_descartes(A, Y, N)[1] == []
    dimer = [("0", "X"), ("X", "0"), ("X+X", "0")]
    A, Y, N, cx = _crn(dimer, ["X"])
    it = kimonu.kinetic_modules(A, Y, N)
    assert kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr).acr == [], "the iterated route alone misses it"
    assert kinetic_modules_descartes(A, Y, N)[1] == [0], "one sign change: exactly one positive x"
    assert kinetic_modules_descartes(A, Y, N, detach_uptakes=True)[1] == [], "a free source leaves x free"
