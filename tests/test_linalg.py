"""The exact linear algebra helpers: reading matrices, nullspaces in row form, proportional rows, blocked reactions."""
from fractions import Fraction

import numpy as np
import straindesign as sd
from scipy import sparse

from kimonu.linalg import (blocked_reactions, entries, exact, left_null_rows, matrix_product, null_space_rows,
                           proportional_classes, to_float)

N = sparse.csr_matrix(np.array([[1, -1, 0, 0], [0, 1, -1, -1]]))      # uptake of a, a -> b, two exports of b


def test_entries_read_every_exact_type_alike():
    exact_N = exact([(0, 0, Fraction(1)), (0, 1, Fraction(-1)), (1, 1, Fraction(1)), (1, 2, Fraction(-1)),
                     (1, 3, Fraction(-1))], (2, 4))
    assert sorted(entries(N)[0]) == sorted(entries(exact_N)[0]) and entries(exact_N)[1] == (2, 4)
    big = sd.sparse_nullspace(exact([(0, 0, Fraction(1, 10**30)), (0, 1, Fraction(1))], (1, 2)))
    (_, _, v), = [t for t in entries(big)[0] if t[0] == 1]
    assert isinstance(v, Fraction)


def test_exact_drops_zero_entries_and_to_float_reads_back():
    M = exact([(0, 0, Fraction(1, 3)), (0, 1, Fraction(0))], (1, 2))
    assert entries(M)[0] == [(0, 0, Fraction(1, 3))]
    assert np.allclose(to_float(M).toarray(), [[1 / 3, 0]])


def test_null_space_rows_and_proportional_classes():
    """ker(N) is two-dimensional. The uptake and a -> b both carry the sum of the two exports, so rows 0 and 1 of K
    are equal, and rows 2 and 3 are proportional to neither them nor each other."""
    rows = null_space_rows(N)
    assert len(rows) == 4 and len({j for r in rows for j in r}) == 2
    Nd = N.toarray()
    for col in range(2):
        k = [r.get(col, 0) for r in rows]
        assert all(sum(int(Nd[i, j]) * k[j] for j in range(4)) == 0 for i in range(2))
    zero, classes = proportional_classes(rows)
    assert zero == [] and [0, 1] in classes and len(classes) == 3


def test_matrix_product_takes_a_matrix_or_rows():
    A = sparse.csr_matrix(np.array([[1, 0], [1, 1], [0, 0]]))
    K = [{0: Fraction(2)}, {1: Fraction(-2), 0: Fraction(1)}]
    assert matrix_product(A, K) == [{0: 2}, {0: 3, 1: -2}, {}]
    assert matrix_product(A, exact([(0, 0, Fraction(2)), (1, 0, Fraction(1)), (1, 1, Fraction(-2))], (2, 2))) == \
        [{0: 2}, {0: 3, 1: -2}, {}]


def test_left_null_rows():
    """ker(Nᵀ) of a matrix with two equal rows: one vector, with opposite entries on the two rows."""
    M = sparse.csr_matrix(np.array([[1, 2], [1, 2]]))
    (a,), (b,) = left_null_rows(M)
    assert a[1] == -b[1] != 0
    assert left_null_rows(sparse.csr_matrix((2, 0))) is None


def test_blocked_reactions_by_one_lp():
    """A chain runs only while its export is open, and a product that nothing consumes blocks its producers."""
    S = sparse.csr_matrix(np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]]))    # up a, a -> b, b -> c, c out
    assert blocked_reactions(S).tolist() == [False] * 4
    assert blocked_reactions(S, closed=[3]).tolist() == [True] * 4
    dead_end = sparse.csr_matrix(np.array([[1, -1, 0], [0, 1, 0], [0, 0, 1]]))           # b and the third are dead
    assert blocked_reactions(dead_end).tolist() == [True, True, True]
