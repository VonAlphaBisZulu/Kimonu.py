"""Absolute concentration robustness, exactly.

A metabolite has **ACR** when its concentration is the same in every positive steady state; a pair
has **ACRR** when their ratio is. Following Langary et al. (2025), both are read off the composition
differences within kinetic modules: collect the columns Y(C_i) - Y(C_j) for complexes C_i, C_j in a
common kinetic module into a matrix YΔ, and then

    metabolite i has ACR         <=>  e_i        in im(YΔ)
    metabolites i, j have ACRR   <=>  e_i - e_j  in im(YΔ)

Span membership is decided here through the left nullspace rather than by projection. With
L = ker(YΔᵀ), a vector v lies in im(YΔ) exactly when Lᵀv = 0, so

    ACR   <=>  row i of L is zero
    ACRR  <=>  rows i and j of L are equal

which is the same primitive the concordance core already uses — one exact nullspace, then a grouping
of identical rows. It costs one kernel computation instead of a rank test per candidate, it makes
ACRR a grouping rather than a sweep over metabolite pairs, and it is exact: the alternative route
through a QR decomposition needs a numerical rank tolerance, and a tolerance is precisely what
decides borderline cases here.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np
from scipy import sparse

import straindesign as sd

AcrResult = namedtuple("AcrResult", ["acr", "acrr", "n_columns"])


def composition_differences(Y, modules, known_acr=()):
    """YΔ: the composition differences within each kinetic module.

    Each module contributes the differences of its first complex against every other member, which
    spans the same subspace as all pairwise differences would.

    Metabolites already known to have ACR contribute their unit vector as an extra column
    (Remark S3-6). This is what lets the analysis make progress: once a metabolite is known to be
    robust, differences involving it become decidable, which can license merges that were not
    decidable before.
    """
    Y = Y.tocsc()
    columns = []
    for module in modules:
        members = sorted(module)
        if len(members) < 2:
            continue
        reference = Y[:, members[0]]
        for other in members[1:]:
            columns.append(reference - Y[:, other])
    for met in sorted(known_acr):
        unit = sparse.csc_matrix(([1], ([met], [0])), shape=(Y.shape[0], 1), dtype=np.int64)
        columns.append(unit)
    if not columns:
        return sparse.csc_matrix((Y.shape[0], 0), dtype=np.int64)
    return sparse.hstack(columns).tocsc()


def _rows_of_left_nullspace(YD):
    """Rows of L = ker(YΔᵀ), as hashable exact tuples."""
    if YD.shape[1] == 0:
        # Nothing is in the span of an empty matrix except the zero vector, so no metabolite has
        # ACR; the left nullspace is the whole space and every row is distinct.
        return None

    L = sd.sparse_nullspace(YD.T.tocsr())
    if isinstance(L, sd.ExactCOO):
        rows = {}
        for r, c, v in zip(L.rows, L.cols, L.data):
            rows.setdefault(int(r), {})[int(c)] = Fraction(int(v))
        return [tuple(sorted(rows.get(i, {}).items())) for i in range(YD.shape[0])]

    L = L.tocsr()
    out = []
    for i in range(L.shape[0]):
        start, end = L.indptr[i], L.indptr[i + 1]
        order = np.argsort(L.indices[start:end])
        out.append(tuple((int(L.indices[start:end][k]), Fraction(int(L.data[start:end][k])))
                         for k in order))
    return out


def coupling_signatures(Y, modules, known_acr=()):
    """Signature L'Y[:, c] of each module's representative, for the Proposition S3-4 merge.

    Two modules may be merged when the composition difference of their representatives lies in
    im(YΔ), which holds exactly when their signatures agree — so merging is again a grouping of
    equal exact vectors rather than a rank test per pair.
    """
    YD = composition_differences(Y, modules, known_acr)
    rows = _rows_of_left_nullspace(YD)
    if rows is None:
        return None

    Ycsc = Y.tocsc()

    def signature(complex_index):
        start, end = Ycsc.indptr[complex_index], Ycsc.indptr[complex_index + 1]
        acc = {}
        for met, coeff in zip(Ycsc.indices[start:end], Ycsc.data[start:end]):
            for col, value in rows[int(met)]:
                acc[col] = acc.get(col, Fraction(0)) + Fraction(str(coeff)) * value
        return tuple(sorted((c, v) for c, v in acc.items() if v != 0))

    return [signature(min(module)) for module in modules]


def acr_acrr(Y, modules, metabolite_ids=None, known_acr=()):
    """ACR metabolites and ACRR metabolite pairs from the kinetic modules.

    Args:
        Y: metabolites x complexes composition matrix.
        modules: kinetic modules as iterables of complex indices.
        metabolite_ids: optional names; indices are returned when omitted.

    ACRR pairs where both metabolites already have ACR are excluded: their ratio is constant only
    because each concentration is, which is not an independent finding.
    """
    YD = composition_differences(Y, modules, known_acr)
    rows = _rows_of_left_nullspace(YD)
    n = Y.shape[0]
    names = metabolite_ids if metabolite_ids is not None else list(range(n))

    if rows is None:
        return AcrResult(acr=[], acrr=[], n_columns=0)

    acr = [i for i in range(n) if not rows[i]]
    acr_set = set(acr)

    groups = {}
    for i in range(n):
        if i in acr_set or not rows[i]:
            continue
        groups.setdefault(rows[i], []).append(i)

    acrr = []
    for members in groups.values():
        for a in range(len(members)):
            for b in range(a + 1, len(members)):
                acrr.append((names[members[a]], names[members[b]]))

    return AcrResult(acr=[names[i] for i in acr], acrr=sorted(acrr), n_columns=YD.shape[1])
