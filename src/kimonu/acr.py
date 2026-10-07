"""Fixed log-combinations: ACR, ACRR and kinetic modules, read off one exact basis.

Some combinations vᵀ ln x of the log concentrations may take the same value at every positive steady
state. A metabolite has **absolute concentration robustness** (ACR) when ln x_i is such a
combination, a pair has **ACRR** when ln x_i - ln x_j is, and two complexes are **kinetically coupled**
when ln ψ_c - ln ψ_d = (y_c - y_d)ᵀ ln x is, with y_c the composition of complex c. All three questions
are membership in one subspace C of fixed combinations:

    metabolite i has ACR         <=>  e_i        in C
    metabolites i, j have ACRR   <=>  e_i - e_j  in C
    complexes c, d are coupled   <=>  y_c - y_d  in C

Both routes to kinetic modules produce such a C. The iterated route (kinetic.py) collects the
composition differences within its modules, plus the unit vectors of metabolites already shown to have
ACR, as the columns of a matrix YΔ (Langary et al., 2025) and takes C = im(YΔ). The tangent route
(tangent.py) takes C as the orthogonal complement of the tangent space T of the steady states.

Membership is decided through a basis B of the complement C^⊥ rather than by projection: v is in C
exactly when Bᵀ v = 0. For the tangent route B = T. For the iterated route B = L = ker(YΔᵀ), by the
lemma

    v in im(YΔ)  <=>  Lᵀ v = 0,

which holds since im(YΔ) is the orthogonal complement of ker(YΔᵀ) (exactly, over the rationals). So

    ACR   <=>  row i of B is zero
    ACRR  <=>  rows i and j of B are equal
    coupled complexes  <=>  rows c and d of Yᵀ B are equal

which is the primitive the concordance code already uses: one exact null space, then a grouping of
identical rows. It costs one kernel computation instead of a span test per candidate, it makes ACRR and
the modules a grouping rather than a sweep over pairs, and it is exact; a QR-based rank test needs a
numerical tolerance, and a tolerance is precisely what decides borderline cases here.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np
from scipy import sparse

from kimonu.linalg import exact, null_space_rows

AcrResult = namedtuple("AcrResult", ["acr", "acrr", "n_columns"])
Readout = namedtuple("Readout", ["modules", "acr", "acrr"])


def composition_differences(Y, modules, known_acr=()):
    """YΔ: the composition differences within each kinetic module.

    Each module contributes the differences of its first complex against every other member, which
    spans the same subspace as all pairwise differences would.

    Metabolites already known to have ACR contribute their unit vector as an extra column (Remark S3-6
    of Langary et al., 2025). This is what lets the analysis make progress: once a metabolite is known
    to be robust, differences involving it become decidable, which can license merges that were not
    decidable before.

    Returned as a scipy matrix, integer for an integer Y and the nearest floats otherwise, for
    inspection and the benchmark scripts; free_directions works from the exact columns instead.
    """
    columns = _difference_columns(Y, modules, known_acr)
    integer = np.issubdtype(Y.dtype, np.integer)
    triples = [(i, k, int(v) if integer else float(v)) for k, col in enumerate(columns) for i, v in col.items()]
    return sparse.csc_matrix(([v for *_, v in triples], ([i for i, *_ in triples], [k for _, k, _ in triples])),
                             shape=(Y.shape[0], len(columns)), dtype=np.int64 if integer else np.float64)


def _difference_columns(Y, modules, known_acr):
    """The columns of YΔ as exact {metabolite: value} dicts.

    Exact because a float cannot hold them: biomass coefficients come with 15 significant digits, and
    a difference such as 9.86915907870449e-05 - 4.11265742135554e-05 needs more. The rounded value
    would enter the exact null space as if it were the true one.
    """
    Y = sparse.csc_matrix(Y)

    def column(c):
        return {int(i): Fraction(str(v)) for i, v in zip(Y.indices[Y.indptr[c]:Y.indptr[c + 1]],
                                                       Y.data[Y.indptr[c]:Y.indptr[c + 1]])}

    columns = []
    for module in modules:
        members = sorted(module)
        if len(members) < 2:
            continue
        reference = column(members[0])
        for other in members[1:]:
            difference = dict(reference)
            for i, v in column(other).items():
                difference[i] = difference.get(i, 0) - v
            columns.append({i: v for i, v in difference.items() if v != 0})
    return columns + [{met: Fraction(1)} for met in sorted(known_acr)]


def free_directions(Y, modules, known_acr=()):
    """L = ker(YΔᵀ) in row form, one hashable row per metabolite: the basis B of the iterated route.

    With no composition difference to span, nothing is fixed and L is the identity.
    """
    return _kernel_rows(_difference_columns(Y, modules, known_acr), Y.shape[0])


def _kernel_rows(columns, m):
    """Rows of an exact basis of ker(YΔᵀ) for the columns of YΔ; the identity when there are none."""
    if not columns:
        return [((i, Fraction(1)),) for i in range(m)]
    YDt = exact(((k, i, v) for k, col in enumerate(columns) for i, v in col.items()), (len(columns), m))
    return [tuple(sorted(r.items())) for r in null_space_rows(YDt)]


def _combined(Ycsc, rows, c):
    """Row c of Yᵀ B: the rows of B for complex c's metabolites, weighted by its composition."""
    acc = {}
    for met, coeff in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]],
                          Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
        for col, value in rows[int(met)]:
            acc[col] = acc.get(col, Fraction(0)) + Fraction(str(coeff)) * value
    return tuple(sorted((col, v) for col, v in acc.items() if v != 0))


def readout(Y, rows, skip=()):
    """Kinetic modules, ACR metabolites and ACRR groups from a basis B of the free directions.

    Args:
        Y: metabolites x complexes composition matrix.
        rows: B in row form, one hashable exact row per metabolite (sorted (column, value) pairs), as
            free_directions returns it for the iterated route and tangent.py for the tangent route.
        skip: complexes left out of the modules (the zero complex when uptakes are detached).

    Returns a :class:`Readout`: ``modules`` are the groups of complexes with equal rows of Yᵀ B, in
    order of their first complex, singletons included and skipped complexes left out; ``acr`` are
    the metabolites with a zero row of B; ``acrr`` are the groups of at least two non-ACR metabolites
    with equal rows.
    """
    Ycsc = sparse.csc_matrix(Y)
    # 1) Modules: complexes with equal rows of Yᵀ B, i.e. whose composition differences B annihilates
    groups = {}
    for c in range(Ycsc.shape[1]):
        if c not in skip:
            groups.setdefault(_combined(Ycsc, rows, c), set()).add(c)
    # 2) ACR and ACRR: zero and equal nonzero rows of B
    acr, acrr = _robust(rows)
    return Readout(modules=list(groups.values()), acr=acr, acrr=acrr)


def _robust(rows):
    """ACR metabolites (zero rows of B) and ACRR groups (equal nonzero rows of B).

    Pairs of ACR metabolites are left out of ACRR: their ratio is constant only because each
    concentration is, which is not an independent finding.
    """
    acr = [i for i, row in enumerate(rows) if not row]
    same = {}
    for i, row in enumerate(rows):
        if row:
            same.setdefault(row, []).append(i)
    return acr, [g for g in same.values() if len(g) > 1]


def coupling_signatures(Y, modules, known_acr=()):
    """Signature Lᵀ y_c of each module's representative, for the Proposition S3-4 merge.

    Two modules may be merged when the composition difference of their representatives lies in
    im(YΔ), which holds exactly when their signatures agree (Proposition S3-4 of Langary et al.,
    2025). kinetic.py reads the same keys for every complex through :func:`readout`; this
    per-representative form is kept for the benchmark scripts. None when YΔ has no columns.
    """
    columns = _difference_columns(Y, modules, known_acr)
    if not columns:
        return None
    rows, Ycsc = _kernel_rows(columns, Y.shape[0]), sparse.csc_matrix(Y)
    return [_combined(Ycsc, rows, min(module)) for module in modules]


def acr_acrr(Y, modules, metabolite_ids=None, known_acr=()):
    """ACR metabolites and ACRR metabolite pairs from the kinetic modules (the iterated route).

    Args:
        Y: metabolites x complexes composition matrix.
        modules: kinetic modules as iterables of complex indices.
        metabolite_ids: optional names; indices are returned when omitted.
        known_acr: metabolites already known to have ACR (Remark S3-6).

    This is the readout of :func:`readout` on the metabolites only, with ACRR given as pairs.
    """
    columns = _difference_columns(Y, modules, known_acr)
    names = metabolite_ids if metabolite_ids is not None else list(range(Y.shape[0]))
    if not columns:
        return AcrResult(acr=[], acrr=[], n_columns=0)
    acr, groups = _robust(_kernel_rows(columns, Y.shape[0]))
    acrr = [(names[a], names[b]) for g in groups for k, a in enumerate(g) for b in g[k + 1:]]
    return AcrResult(acr=[names[i] for i in acr], acrr=sorted(acrr), n_columns=len(columns))
