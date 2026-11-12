"""Exact linear algebra behind the analyses: the row form of a matrix, nullspaces, blocked reactions.

The row form. Every matrix that Kimonu computes, K, M = A·K, W, Yᵀ, (YΔ)ᵀ, B and Yᵀ B, is held in one form, a list
with one entry per row of the matrix. Each row is a tuple of (column, value) pairs for its nonzero entries, sorted
by column, with exact values (int or Fraction). The matrix

    [[1, 0, -1/2],
     [0, 0,    0],
     [0, 3,    0]]

is [((0, 1), (2, Fraction(-1, 2))), (), ((1, 3),)]. A zero row is the empty tuple, and two rows are equal exactly
when their tuples are, so zero rows, equal rows and proportional rows are found by grouping the rows in a
dictionary, without a tolerance (:func:`proportional_classes`, :func:`kimonu.readout.readout`). The number of
columns is not stored, and :func:`null_space_rows` takes it with a row form. The products of the paper, A·K, Yᵀ B
and W Yᵀ, are all :func:`matrix_product`.

The network's own matrices keep the containers they come in. A (complexes x reactions) is a scipy integer matrix,
and N is a straindesign RationalMatrix, the input of every nullspace of N. Y is held exactly, as Yᵀ in row form
(:class:`kimonu.CRN`).

straindesign's exact matrices come in three types, a RationalMatrix, and from sparse_nullspace either a scipy
integer matrix (entries fit in int64) or an ExactCOO (they do not). :func:`entries` reads all of them, and scipy
and numpy matrices too, and no other function checks the type. An ExactCOO holds every entry over one common
denominator, the lcm of all denominators in the matrix, which on a kernel at a generic flux is a huge integer, so
kernels are read as fractions from the RationalMatrix itself (:func:`null_space_rows`).
"""
from fractions import Fraction
from numbers import Integral

import numpy as np
import straindesign as sd
from scipy import sparse
from straindesign.compression import nullspace
from scipy.optimize import linprog


def fraction(x):
    """The exact rational value of a number as it is written, an int or Fraction as it is and a float as its decimal.

    ``Fraction(str(x))`` reproduces the decimal that a model declares (1.496 -> 1496/1000), whereas ``Fraction(x)``
    would recover the binary approximation of it. Coefficients are frequently fractional, as in biomass reactions and
    half-oxygen terms such as 0.5 O2 in CYTBD, so this is the difference between analysing the model and analysing a
    rounded version of it. Every number that enters Kimonu from outside is read here.
    """
    if isinstance(x, Fraction):
        return x
    return Fraction(int(x)) if isinstance(x, Integral) else Fraction(str(x))


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
    return [(int(i), int(j), fraction(v)) for i, j, v in zip(M.row, M.col, M.data)], M.shape


def exact(triples, shape):
    """A straindesign RationalMatrix from (row, column, value) triples, every value read by :func:`fraction`.

    Zero entries are dropped. Clearing denominators to int64 instead would overflow on models whose coefficients
    are written as 15-significant-digit decimals. iYS854's -7.73333333333333e-08 has a denominator of 10**22, and
    Chlorella_variabilis_iAJ526 mixes seven such terms in one row.
    """
    return sd.RationalMatrix.from_fractions(((int(i), int(j), fraction(v)) for i, j, v in triples), shape)


def row_form(M):
    """The row form of any matrix that :func:`entries` reads."""
    triples, (m, _) = entries(M)
    rows = [[] for _ in range(m)]
    for i, j, v in triples:
        if v:
            rows[i].append((j, v))
    return [tuple(sorted(r)) for r in rows]


def to_float(M):
    """M as a scipy float matrix."""
    triples, shape = entries(M)
    return sparse.csr_matrix(([float(v) for *_, v in triples], ([i for i, *_ in triples], [j for _, j, _ in triples])),
                             shape=shape)


def matrix_product(X, Z):
    """X·Z in row form. Each factor is a row form or any matrix that :func:`entries` reads.

    Row i of X·Z is the sum of the rows of Z weighted by the entries of row i of X, so the product touches only
    nonzero entries. A·K, Yᵀ B and W Yᵀ are all computed this way.
    """
    X = X if isinstance(X, list) else row_form(X)
    Z = Z if isinstance(Z, list) else row_form(Z)
    product = []
    for row in X:
        acc = {}
        for k, x in row:
            for col, z in Z[k]:
                acc[col] = acc.get(col, 0) + x * z
        product.append(tuple(sorted((col, v) for col, v in acc.items() if v)))
    return product


def null_space_rows(M, n_columns=None):
    """An exact basis of ker(M) in row form, one row per column of M and one column per basis vector.

    M is a row form with ``n_columns`` columns, or any matrix that :func:`entries` reads. The basis is
    straindesign's RREF basis, with 1 at each free coordinate and 0 at the other free coordinates, so the row of a
    free coordinate has a single entry, 1 at its own basis vector. A matrix without rows has the identity as its
    basis.
    """
    if isinstance(M, list):
        M = sd.RationalMatrix.from_fractions(((i, j, v) for i, row in enumerate(M) for j, v in row),
                                             (len(M), n_columns))
    elif not isinstance(M, sd.RationalMatrix):
        M = exact(*entries(M))
    K = nullspace(M)
    rows = [()] * M.get_column_count()
    if K.get_column_count():
        for i, row in K.to_sparse_pattern()[1].items():
            rows[int(i)] = tuple(sorted((int(j), v) for j, v in row.items()))
    return rows


def proportional_classes(rows):
    """Indices of the zero rows, and the classes of mutually proportional nonzero rows, of a row form.

    Dividing each row by its first nonzero entry makes proportional rows identical, so grouping equal keys finds
    every class at once (and is transitive by construction). The classes are in the order of their first row.
    """
    zero, classes = [], {}
    for i, row in enumerate(rows):
        if not row:
            zero.append(i)
            continue
        first = row[0][1]
        classes.setdefault(tuple((col, v / first) for col, v in row), []).append(i)
    return zero, list(classes.values())


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
