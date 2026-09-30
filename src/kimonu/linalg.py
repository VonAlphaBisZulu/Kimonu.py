"""Exact sparse matrices, whatever type straindesign hands back.

straindesign's exact matrices come in three types: a RationalMatrix (what Kimonu builds N as), and from
sparse_nullspace either a scipy integer matrix (entries fit in int64) or an ExactCOO (they do not).
The rest of Kimonu reads matrices only through these two functions and never checks the type.
"""
from fractions import Fraction

from scipy import sparse

import straindesign as sd


def entries(M):
    """(row, column, Fraction) triples of M, and its shape."""
    if isinstance(M, sd.RationalMatrix):
        M = M.to_coo_exact()
    if isinstance(M, sd.ExactCOO):
        return [(int(i), int(j), Fraction(int(v), int(M.denom)))
                for i, j, v in zip(M.rows, M.cols, M.data)], tuple(M.shape)
    M = sparse.coo_matrix(M)
    return [(int(i), int(j), Fraction(str(v))) for i, j, v in zip(M.row, M.col, M.data)], M.shape


def product_rows(A, K):
    """Rows of A·K as {column: value}, one per row of A; a zero row is an empty dict.

    A is an integer matrix (the complex incidence matrix), K any exact matrix.
    """
    k_rows = {}
    for i, j, v in entries(K)[0]:
        k_rows.setdefault(i, []).append((j, v))
    A = sparse.csr_matrix(A)
    rows = []
    for c in range(A.shape[0]):
        acc = {}
        for j, a in zip(A.indices[A.indptr[c]:A.indptr[c + 1]], A.data[A.indptr[c]:A.indptr[c + 1]]):
            for col, v in k_rows.get(int(j), ()):
                acc[col] = acc.get(col, 0) + int(a) * v
        rows.append({col: v for col, v in acc.items() if v != 0})
    return rows
