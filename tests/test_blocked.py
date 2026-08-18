"""LP-free blocked-reaction detection must never over-report.

Soundness is the property that matters: everything flagged has to be genuinely blocked, because
callers use the result to prune. Completeness is not claimed — deciding the irreversibility half
exactly is an implicit-equality question for a polyhedral cone.
"""
import numpy as np
from scipy import sparse

from kimonu.blocked import blocked_reactions, stoichiometrically_blocked


def test_zero_row_of_the_kernel_is_blocked():
    """A reaction consuming a metabolite nothing produces cannot carry flux."""
    #  r0: -> A      r1: A ->      r2: B ->      (B has no producer)
    N = sparse.csr_matrix(np.array([[1, -1, 0], [0, 0, -1]], dtype=np.int64))
    result = blocked_reactions(N, [0, 0, 0], [1000, 1000, 1000])
    assert 2 in result.blocked
    assert 2 in stoichiometrically_blocked(N)
    assert 0 not in result.blocked and 1 not in result.blocked


def test_reversible_reaction_pinned_from_two_sides():
    """Two balances can force a reversible reaction to zero even though neither does alone.

    e_e leaves only through an export-only exchange, so the reversible transporter must run
    backwards to supply it; the e_c balance then requires it to run forwards.
    """
    #            EX(e_e ->)   T(e_e <=> e_c)   C(e_c ->)
    N = sparse.csr_matrix(np.array([[-1, -1, 0],
                                    [0, 1, -1]], dtype=np.int64))
    result = blocked_reactions(N, [0, -1000, 0], [1000, 1000, 1000])
    assert set(result.blocked) == {0, 1, 2}


def test_backward_only_reactions_are_oriented_correctly():
    """A reaction confined to negative flux runs backwards; its coefficients read inverted.

    Ignoring that inverts the direction of the dead-end argument and blocks reactions that are
    perfectly able to carry flux.
    """
    #  r0: -> A (forward)     r1: A -> B written backwards, bounds (-1000, 0), so it makes A
    N = sparse.csr_matrix(np.array([[1, -1], [0, 1]], dtype=np.int64))
    result = blocked_reactions(N, [0, -1000], [1000, 0])
    assert result.blocked == [] or 0 not in result.blocked
