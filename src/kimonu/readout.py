"""The basis B and its readout: ACR, ACRR and kinetic modules as zero rows and equal rows (Fig. 1F).

Some combinations vᵀ ln x of the log concentrations take the same value at every positive steady state with the
same rate constants. They are called fixed combinations, and they form a linear subspace, the span of a matrix V.
A metabolite i has **absolute concentration robustness** (ACR) when ln x_i is fixed, a pair i, j has **ACRR** when
ln x_i - ln x_j is fixed, and two complexes c, d are **kinetically coupled** when ln ψ_c - ln ψ_d = (y_c - y_d)ᵀ ln x
is fixed (Eq. 2), with y_c the composition of complex c:

    metabolite i has ACR         <=>  e_i        in span(V)
    metabolites i, j have ACRR   <=>  e_i - e_j  in span(V)
    complexes c, d are coupled   <=>  y_c - y_d  in span(V)

Membership is decided through a basis B of the complementary directions, those in which ln x can differ between
two steady states, B = ker(Vᵀ), rather than by projection. v lies in span(V) exactly when Bᵀ v = 0, so

    ACR      <=>  row i of B is zero
    ACRR     <=>  rows i and j of B are equal
    coupled  <=>  rows c and d of Yᵀ B are equal

and :func:`readout` finds all three by grouping identical exact rows, as concordance groups the rows of M = A·K.
It costs one kernel computation instead of a span test per candidate, it finds every member of an ACRR group or a
kinetic module at once, and it is exact. A rank test in floating point needs a tolerance, and a tolerance is
precisely what decides borderline cases here.

The true B cannot be computed from the network directly. Two routes bound it with nullspaces (Eq. 4),

    ker((YΔ)ᵀ)  ⊇  B  ⊇  ker(W Yᵀ),

and both are read out here. The iterated route (:mod:`kimonu.iterated`) collects composition differences that it
has proven fixed as the columns of YΔ and takes ker((YΔ)ᵀ) (:func:`free_directions`). Every column of YΔ lies in
span(V), so ker((YΔ)ᵀ) contains B, and every zero row or pair of equal rows it shows is also one of B. The
tangent route (:mod:`kimonu.tangent`) takes ker(W Yᵀ), which is contained in B, so it shows every relation of B
and possibly false couplings. A larger nullspace has fewer zero rows and equal rows and reports fewer relations.

The reference-equivalent mode (:mod:`kimonu.reference`) reads ACR and ACRR by the pairwise Shinar–Feinberg rule
instead of from a basis. All three routes return the same :class:`RouteResult`.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np
from scipy import sparse

from kimonu.linalg import exact, null_space_rows

Readout = namedtuple("Readout", ["modules", "acr", "acrr_groups"])
Readout.__doc__ = """What :func:`readout` reads off a basis B, in indices.

modules      groups of complexes with equal rows of Yᵀ B, in the order of their first complex, singletons included
acr          metabolites with a zero row of B, in index order
acrr_groups  groups of two or more metabolites without ACR whose rows of B are equal, each in index order
"""

RouteResult = namedtuple("RouteResult", ["route", "modules", "module_of", "n_modules", "giant_module", "acr", "acrr",
                                         "B", "dimension", "flux", "rounds", "concordance"])
RouteResult.__doc__ = """The relations one route reports, in one layout for all three routes.

route         "iterated", "tangent" or "reference"
modules       kinetic modules as sets of complex indices. They partition the complexes, singletons included, and
              with detached uptakes the zero complex is a module of its own, listed last.
module_of     numpy array, the module of each complex, counted from 1
n_modules     number of modules with two or more complexes
giant_module  size of the largest module
acr           metabolites with ACR, as names when metabolite_ids is given and as indices otherwise
acrr          sorted pairs of metabolites with ACRR, as names or indices. Pairs of two ACR metabolites are left
              out, since their ratio is fixed only because each concentration is.
B             the basis the relations were read from, one row per metabolite, each row a sorted tuple of
              (basis vector, exact value). ker((YΔ)ᵀ) for the iterated route, the metabolite rows of ker(W Yᵀ)
              for the tangent route, None for the reference route.
dimension     number of basis vectors of that nullspace. With detached uptakes, ker(W Yᵀ) also has the source
              coordinates, so for the tangent route this can exceed the dimension of B in metabolite space.
              None for the reference route.
flux          the interior flux r* the route ran at (rate constants k = r* with the steady state x = 1, see
              :func:`kimonu.rate_constants`), None for the reference route
rounds        rounds of the iterated route, or passes of the overlap merge of the reference route; None for the
              tangent route, which needs no iteration
concordance   the reference route's :class:`kimonu.ConcordanceResult`, None for the other routes
"""


def composition_differences(Y, modules, known_acr=()):
    """YΔ, the composition differences within each module, as a scipy matrix for inspection.

    Δ (complexes x pairs) has one column e_c - e_d for each member c of a module other than its first member d, which
    spans the same subspace as all pairwise differences would. Metabolites already known to have ACR contribute
    their unit vector as an extra column, as in Langary et al. (2025). Once a metabolite is known to be fixed,
    differences involving it become decidable, which can license merges that were not decidable before.

    Integer for an integer Y and the nearest floats otherwise. :func:`free_directions` works from the exact columns.
    """
    columns = _difference_columns(Y, modules, known_acr)
    integer = np.issubdtype(Y.dtype, np.integer)
    triples = [(i, k, int(v) if integer else float(v)) for k, col in enumerate(columns) for i, v in col.items()]
    return sparse.csc_matrix(([v for *_, v in triples], ([i for i, *_ in triples], [k for _, k, _ in triples])),
                             shape=(Y.shape[0], len(columns)), dtype=np.int64 if integer else np.float64)


def _difference_columns(Y, modules, known_acr):
    """The columns of YΔ as exact {metabolite: value} dicts.

    Exact because a float cannot hold them. Biomass coefficients come with 15 significant digits, and a difference
    such as 9.86915907870449e-05 - 4.11265742135554e-05 needs more. The rounded value would enter the exact
    nullspace as if it were the true one.
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
    """ker((YΔ)ᵀ) in row form, one hashable row per metabolite, the basis B of the iterated route (Eq. 4, left).

    Each module contributes the differences of its first complex against every other member as columns of YΔ, and
    each metabolite in ``known_acr`` its unit vector (:func:`composition_differences`). v lies in the span of YΔ
    exactly when it is orthogonal to ker((YΔ)ᵀ), exactly, over the rationals. With no column to span, nothing is
    fixed and the basis is the identity.
    """
    columns, m = _difference_columns(Y, modules, known_acr), Y.shape[0]
    if not columns:
        return [((i, Fraction(1)),) for i in range(m)]
    YDt = exact(((k, i, v) for k, col in enumerate(columns) for i, v in col.items()), (len(columns), m))
    return [tuple(sorted(r.items())) for r in null_space_rows(YDt)]


def complex_directions(Y, B):
    """Yᵀ B in row form, one hashable row per complex: how ln ψ_c can differ between steady states.

    Row c is the rows of B for complex c's metabolites, weighted by its composition. Complexes with equal rows are
    kinetically coupled, and :func:`kimonu.kinetic_flux_coupling` lifts the same rows into flux space.
    """
    Ycsc = sparse.csc_matrix(Y)
    return [_combined(Ycsc, B, c) for c in range(Ycsc.shape[1])]


def _combined(Ycsc, rows, c):
    """Row c of Yᵀ B, the rows of B for complex c's metabolites weighted by its composition."""
    acc = {}
    for met, coeff in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]],
                          Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
        for col, value in rows[int(met)]:
            acc[col] = acc.get(col, Fraction(0)) + Fraction(str(coeff)) * value
    return tuple(sorted((col, v) for col, v in acc.items() if v != 0))


def readout(Y, B, skip=()):
    """Kinetic modules, ACR metabolites and ACRR groups from a basis B, by grouping identical rows.

    Args:
        Y: metabolites x complexes composition matrix.
        B: the basis in row form, one hashable exact row per metabolite (a sorted tuple of (basis vector, value)
            pairs), as :func:`free_directions` returns it for the iterated route and
            :func:`kimonu.tangent.tangent_space` for the tangent route.
        skip: complexes left out of the modules (the zero complex when uptakes are detached).

    Returns a :class:`Readout` in indices.
    """
    Ycsc = sparse.csc_matrix(Y)
    # 1) Modules: complexes with equal rows of Yᵀ B, whose composition differences B annihilates
    groups = {}
    for c in range(Ycsc.shape[1]):
        if c not in skip:
            groups.setdefault(_combined(Ycsc, B, c), set()).add(c)
    # 2) ACR and ACRR: zero rows and equal nonzero rows of B
    acr, acrr_groups = _robust(B)
    return Readout(modules=list(groups.values()), acr=acr, acrr_groups=acrr_groups)


def _robust(rows):
    """ACR metabolites (zero rows of B) and ACRR groups (equal nonzero rows of B).

    Pairs of ACR metabolites are left out of ACRR. Their ratio is constant only because each concentration is,
    which is not an independent finding.
    """
    acr = [i for i, row in enumerate(rows) if not row]
    same = {}
    for i, row in enumerate(rows):
        if row:
            same.setdefault(row, []).append(i)
    return acr, [g for g in same.values() if len(g) > 1]


def _pairs(groups):
    """All pairs (a, b), a before b, within each group, in the order of the groups."""
    return [(a, b) for g in groups for k, a in enumerate(g) for b in g[k + 1:]]


def _route_result(route, n_complexes, modules, acr, acrr, metabolite_ids=None, B=None, dimension=None, flux=None,
                  rounds=None, concordance=None):
    """The :class:`RouteResult` of a route, from its modules (indices), ACR (indices) and ACRR (index pairs)."""
    module_of = np.zeros(n_complexes, dtype=np.int64)
    for idx, group in enumerate(modules, start=1):
        for c in group:
            module_of[c] = idx
    name = (lambda i: i) if metabolite_ids is None else (lambda i: metabolite_ids[i])
    return RouteResult(
        route=route,
        modules=modules,
        module_of=module_of,
        n_modules=sum(1 for g in modules if len(g) > 1),
        giant_module=max((len(g) for g in modules), default=0),
        acr=[name(i) for i in sorted(acr)],
        acrr=sorted((name(a), name(b)) for a, b in acrr),
        B=B,
        dimension=dimension,
        flux=flux,
        rounds=rounds,
        concordance=concordance,
    )
