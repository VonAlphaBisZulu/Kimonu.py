"""Exact linear algebra behind the analyses: reading matrices, nullspaces in row form, blocked reactions.

straindesign's exact matrices come in three types: a RationalMatrix (what Kimonu builds N as), and from
sparse_nullspace either a scipy integer matrix (entries fit in int64) or an ExactCOO (they do not).
Everything here reads them through entries() and never checks the type elsewhere.

An ExactCOO holds every entry over one common denominator, the lcm of all denominators in the matrix. On a
kernel at a generic flux that lcm is huge, so every entry becomes a huge integer. Kernels are therefore read
as fractions from the RationalMatrix itself (null_space_rows).
"""
from fractions import Fraction

import numpy as np
import straindesign as sd
from scipy import sparse
from straindesign.compression import nullspace
from scipy.optimize import linprog


def entries(M):
    """(row, column, Fraction) triples of M, and its shape."""
    if isinstance(M, sd.RationalMatrix):
        shape = (M.get_row_count(), M.get_column_count())
        if 0 in shape:                                   # an empty kernel has no sparse storage to read
            return [], shape
        return [(int(i), int(j), v) for i, row in M.to_sparse_pattern()[1].items() for j, v in row.items()], shape
    if isinstance(M, sd.ExactCOO):
        return [(int(i), int(j), Fraction(int(v), int(M.denom)))
                for i, j, v in zip(M.rows, M.cols, M.data)], tuple(M.shape)
    M = sparse.coo_matrix(M)
    return [(int(i), int(j), Fraction(str(v))) for i, j, v in zip(M.row, M.col, M.data)], M.shape


def to_float(M):
    """M as a scipy float matrix."""
    triples, shape = entries(M)
    return sparse.csr_matrix(([float(v) for *_, v in triples], ([i for i, *_ in triples], [j for _, j, _ in triples])),
                             shape=shape)


def exact(triples, shape):
    """A straindesign RationalMatrix from (row, column, Fraction) triples; zero entries are dropped."""
    return sd.RationalMatrix.from_fractions(((i, j, v) for i, j, v in triples if v != 0), shape)


def matrix_product(A, K):
    """A·K in row form, one {column: value} dict per row of A, for an integer matrix A (the incidence matrix).

    K is any exact matrix, or a matrix in row form (a list with one {column: value} dict per row).
    """
    k_rows = {}
    if isinstance(K, list):
        for i, row in enumerate(K):
            if row:
                k_rows[i] = list(row.items())
    else:
        for i, j, v in entries(K)[0]:
            k_rows.setdefault(i, []).append((j, v))
    A = sparse.csr_matrix(A)
    # the result in row form: one {column: value} dict per row of A, an empty dict for a zero row
    rows = []
    for c in range(A.shape[0]):
        acc = {}
        for j, a in zip(A.indices[A.indptr[c]:A.indptr[c + 1]], A.data[A.indptr[c]:A.indptr[c + 1]]):
            for col, v in k_rows.get(int(j), ()):
                acc[col] = acc.get(col, 0) + int(a) * v
        rows.append({col: v for col, v in acc.items() if v != 0})
    return rows


def null_space_rows(M):
    """An exact basis of ker(M) in row form: one {basis vector: value} dict per column of M."""
    if not isinstance(M, sd.RationalMatrix):
        M = exact(*entries(M))
    rows = [{} for _ in range(M.get_column_count())]
    for i, j, v in entries(nullspace(M))[0]:
        rows[i][j] = v
    return rows


def proportional_classes(rows):
    """Indices of the zero rows, and the classes of mutually proportional nonzero rows.

    Normalizing each row so its first nonzero is 1 makes proportional rows identical, so grouping equal
    keys finds every class at once (and is transitive by construction).
    """
    zero, classes = [], {}
    for i, row in enumerate(rows):
        if not row:
            zero.append(i)
            continue
        nz = sorted(row.items())
        classes.setdefault(tuple((col, v / nz[0][1]) for col, v in nz), []).append(i)
    return zero, list(classes.values())


def left_null_rows(M):
    """Rows of an exact basis of ker(Mᵀ), one hashable tuple per row of M; None if M has no columns."""
    if M.shape[1] == 0:
        return None
    rows = null_space_rows(sparse.csr_matrix(M).T.tocsr())
    return [tuple(sorted(r.items())) for r in rows]


def blocked_reactions(N, closed=()):
    """Boolean array, True for the reactions of N that no steady state with r >= 0 can run.

    One LP: max sum(z) subject to N r = 0, r >= 0, 0 <= z <= 1 and z <= r. The steady states form a cone,
    and two steady states add up to one that runs the reactions of both, so a single r runs every reaction
    that can run at all. Scaling it brings each of them to at least 1, which makes z_j = 1 for those and
    z_j = 0 for the blocked ones. Maximizing sum(r) under upper bounds instead stops at a vertex and leaves
    reactions that can run at zero. ``closed`` lists the indices of reactions with an upper bound of zero.
    """
    Nf = to_float(N); m, q = Nf.shape
    closed = set(closed)
    bounds = [(0, 0 if j in closed else None) for j in range(q)] + [(0, 1)] * q
    A_eq = sparse.hstack([Nf, sparse.csr_matrix((m, q))], format="csr")
    A_ub = sparse.hstack([-sparse.identity(q), sparse.identity(q)], format="csr")
    lp = linprog(np.r_[np.zeros(q), -np.ones(q)], A_ub=A_ub, b_ub=np.zeros(q), A_eq=A_eq, b_eq=np.zeros(m),
                 bounds=bounds, method="highs")
    if lp.status != 0:
        raise ValueError(f"blocked-reaction LP failed: {lp.message}")
    return lp.x[q:] <= 0.5
