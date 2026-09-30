"""Exact linear algebra behind the analyses: reading matrices, null spaces in row form, interior points.

straindesign's exact matrices come in three types: a RationalMatrix (what Kimonu builds N as), and from
sparse_nullspace either a scipy integer matrix (entries fit in int64) or an ExactCOO (they do not).
Everything here reads them through entries() and never checks the type elsewhere.
"""
from fractions import Fraction

import numpy as np
import straindesign as sd
from scipy import sparse
from scipy.optimize import linprog


def entries(M):
    """(row, column, Fraction) triples of M, and its shape."""
    if isinstance(M, sd.RationalMatrix):
        M = M.to_coo_exact()
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
    """A·K for an integer matrix A (the complex incidence matrix) and any exact matrix K."""
    k_rows = {}
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
    rows = [{} for _ in range(M.get_column_count() if isinstance(M, sd.RationalMatrix) else M.shape[1])]
    for i, j, v in entries(sd.sparse_nullspace(M))[0]:
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


def positive_flux(N, seed=0):
    """An exact vector r in ker(N) with every entry strictly positive.

    An LP point with random lower bounds r >= 1 + U(0, 1) and a random objective is snapped onto the
    exact kernel through its free coordinates, the rows where the basis has a single entry. The
    steady-state set is a cone, so no upper bounds are needed. The LP is solved with HiGHS (bundled
    with scipy): forcing every rate to at least 1 spans rates over about eight orders of magnitude on
    genome-scale models, and HiGHS solved it on every model and seed tried, where SCIP, GLPK and CPLEX
    through straindesign often did not. Raises ValueError when no such vector exists, which means the
    network still has blocked reactions.
    """
    Nf = to_float(N); m, q = Nf.shape
    rng = np.random.default_rng(seed)
    lower = 1 + rng.random(q)
    lp = linprog(rng.random(q), A_eq=Nf, b_eq=np.zeros(m), bounds=[(lo, None) for lo in lower], method="highs")
    if lp.status != 0:
        raise ValueError("no strictly positive steady-state flux; remove blocked reactions first")
    rows = null_space_rows(N)
    free = {}
    for i, row in enumerate(rows):
        if len(row) == 1:
            (j, v), = row.items()
            free.setdefault(j, Fraction(float(lp.x[i])) / v)
    flux = [sum((v * free.get(j, 0) for j, v in row.items()), Fraction(0)) for row in rows]
    if not all(x > 0 for x in flux):
        raise ValueError("snapped flux is not strictly positive")
    return flux
