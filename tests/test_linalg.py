"""The exact linear algebra helpers: the row form, reading matrices, products, nullspaces, and the row readout."""
from fractions import Fraction

import numpy as np
import straindesign as sd
from scipy import sparse

from kimonu.linalg import entries, exact, fraction, matrix_product, nullspace, row_form, to_float
from kimonu.readout import row_readout

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


def test_nullspace_and_row_readout():
    """ker(N) is two-dimensional. The uptake and a -> b both carry the sum of the two exports, so rows 0 and 1 of K
    are equal, and rows 2 and 3 are proportional to neither them nor each other."""
    rows = nullspace(N)
    assert len(rows) == 4 and len({j for r in rows for j, _ in r}) == 2
    assert all(r == tuple(sorted(r)) for r in rows)
    Nd = N.toarray()
    for col in range(2):
        k = [dict(r).get(col, 0) for r in rows]
        assert all(sum(int(Nd[i, j]) * k[j] for j in range(4)) == 0 for i in range(2))
    zero, groups, ratio = row_readout(rows)
    assert zero == [] and [0, 1] in groups and len(groups) == 3 and ratio[1] == 1
    assert nullspace(row_form(N), 4) == rows


def test_row_readout_gives_zero_rows_groups_and_ratios():
    """Rows 0, 2 and 4 are proportional, row 2 = -2 · row 0 and row 4 = 1/2 · row 0, row 1 is zero, and row 3 is
    proportional to no other row. The groups are in the order of their first row."""
    F = Fraction
    rows = [((0, F(2)), (3, F(-1))), (), ((0, F(-4)), (3, F(2))), ((1, F(1)),), ((0, F(1)), (3, F(-1, 2)))]
    found = row_readout(rows)
    assert found.zero == [1] and found.groups == [[0, 2, 4], [3]]
    assert found.ratio == {0: 1, 2: -2, 4: F(1, 2), 3: 1}
    for g in found.groups:
        for i in g:
            assert [(j, found.ratio[i] * v) for j, v in rows[g[0]]] == list(rows[i])


def test_nullspace_of_a_matrix_without_rows_is_the_identity():
    assert nullspace([], 3) == [((0, 1),), ((1, 1),), ((2, 1),)]


def test_matrix_product_takes_row_forms_or_matrices():
    A = sparse.csr_matrix(np.array([[1, 0], [1, 1], [0, 0]]))
    K = [((0, Fraction(2)),), ((0, Fraction(1)), (1, Fraction(-2)))]
    assert matrix_product(A, K) == [((0, 2),), ((0, 3), (1, -2)), ()]
    assert matrix_product(row_form(A), exact([(0, 0, 2), (1, 0, 1), (1, 1, -2)], (2, 2))) == \
        [((0, 2),), ((0, 3), (1, -2)), ()]
