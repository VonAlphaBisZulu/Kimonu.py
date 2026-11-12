"""The exact linear algebra helpers: the row form, reading matrices, nullspaces, proportional rows, blocked reactions."""
from fractions import Fraction

import numpy as np
import straindesign as sd
from scipy import sparse

from kimonu.linalg import (blocked_reactions, entries, exact, fraction, matrix_product, null_space_rows,
                           proportional_classes, row_form, to_float)

N = sparse.csr_matrix(np.array([[1, -1, 0, 0], [0, 1, -1, -1]]))      # uptake of a, a -> b, two exports of b


def test_fraction_reads_the_declared_decimal():
    assert fraction(0.1) == Fraction(1, 10) and fraction(np.float64(1.496)) == Fraction(1496, 1000)
    assert fraction(np.int64(3)) == 3 and fraction(Fraction(2, 3)) == Fraction(2, 3) and fraction(10**30) == 10**30


def test_entries_read_every_exact_type_alike():
    exact_N = exact([(0, 0, 1), (0, 1, -1), (1, 1, 1), (1, 2, -1), (1, 3, -1)], (2, 4))
    assert sorted(entries(N)[0]) == sorted(entries(exact_N)[0]) and entries(exact_N)[1] == (2, 4)
    big = sd.sparse_nullspace(exact([(0, 0, Fraction(1, 10**30)), (0, 1, Fraction(1))], (1, 2)))
    (_, _, v), = [t for t in entries(big)[0] if t[0] == 1]
    assert isinstance(v, Fraction)


def test_exact_keeps_the_declared_decimals_and_drops_zeros():
    """-7.73333333333333e-08 has a denominator of 10**22, which int64 cannot clear; Fractions stay as they are."""
    M = exact([(0, 0, -7.73333333333333e-08), (0, 1, Fraction(1, 3)), (0, 2, 0)], (1, 3))
    assert sorted(entries(M)[0]) == [(0, 0, Fraction("-7.73333333333333e-08")), (0, 1, Fraction(1, 3))]
    assert np.allclose(to_float(M).toarray(), [[-7.73333333333333e-08, 1 / 3, 0]])


def test_row_form_is_sorted_sparse_rows():
    M = np.array([[1, 0, 0.5], [0, 0, 0], [0, 3, 0]])
    assert row_form(M) == [((0, 1), (2, Fraction(1, 2))), (), ((1, 3),)]
    assert row_form(exact(entries(M)[0], (3, 3))) == row_form(sparse.csr_matrix(M))


def test_null_space_rows_and_proportional_classes():
    """ker(N) is two-dimensional. The uptake and a -> b both carry the sum of the two exports, so rows 0 and 1 of K
    are equal, and rows 2 and 3 are proportional to neither them nor each other."""
    rows = null_space_rows(N)
    assert len(rows) == 4 and len({j for r in rows for j, _ in r}) == 2
    assert all(r == tuple(sorted(r)) for r in rows)
    Nd = N.toarray()
    for col in range(2):
        k = [dict(r).get(col, 0) for r in rows]
        assert all(sum(int(Nd[i, j]) * k[j] for j in range(4)) == 0 for i in range(2))
    zero, classes = proportional_classes(rows)
    assert zero == [] and [0, 1] in classes and len(classes) == 3
    assert null_space_rows(row_form(N), 4) == rows


def test_null_space_of_a_matrix_without_rows_is_the_identity():
    assert null_space_rows([], 3) == [((0, 1),), ((1, 1),), ((2, 1),)]


def test_matrix_product_takes_row_forms_or_matrices():
    A = sparse.csr_matrix(np.array([[1, 0], [1, 1], [0, 0]]))
    K = [((0, Fraction(2)),), ((0, Fraction(1)), (1, Fraction(-2)))]
    assert matrix_product(A, K) == [((0, 2),), ((0, 3), (1, -2)), ()]
    assert matrix_product(row_form(A), exact([(0, 0, 2), (1, 0, 1), (1, 1, -2)], (2, 2))) == \
        [((0, 2),), ((0, 3), (1, -2)), ()]


def test_blocked_reactions_by_one_lp():
    """A chain runs only while its export is open, and a product that nothing consumes blocks its producers."""
    S = sparse.csr_matrix(np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]]))    # up a, a -> b, b -> c, c out
    assert blocked_reactions(S).tolist() == [False] * 4
    assert blocked_reactions(S, closed=[3]).tolist() == [True] * 4
    dead_end = sparse.csr_matrix(np.array([[1, -1, 0], [0, 1, 0], [0, 0, 1]]))           # b and the third are dead
    assert blocked_reactions(dead_end).tolist() == [True, True, True]
