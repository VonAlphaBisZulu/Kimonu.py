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
has proven fixed as the columns of YΔ (:func:`composition_differences`) and takes ker((YΔ)ᵀ)
(:func:`free_directions`). Every column of YΔ lies in span(V), so ker((YΔ)ᵀ) contains B, and every zero row or pair
of equal rows it shows is also one of B. The tangent route (:mod:`kimonu.tangent`) takes ker(W Yᵀ), which is
contained in B, so it shows every relation of B and possibly false couplings. A larger nullspace has fewer zero
rows and equal rows and reports fewer relations.

The zero complex has no metabolites, so its row of Yᵀ B is zero, and the readout would couple it to every complex
made of ACR metabolites. With detached uptakes the readout leaves it out of the grouping, since no uptake draws
from it, and reports it as a module of its own, listed last. The modules then still partition the complexes.

The reference-equivalent mode (:mod:`kimonu.reference`) reads ACR and ACRR by the pairwise Shinar–Feinberg rule
instead of from a basis. All three routes return the same :class:`RouteResult` (:func:`route_result`).
"""
from collections import namedtuple

import numpy as np

from kimonu.io import zero_complexes
from kimonu.linalg import matrix_product, null_space_rows

Readout = namedtuple("Readout", ["modules", "acr", "acrr_groups"])
Readout.__doc__ = """What :func:`readout` reads off a basis B, in indices.

modules      kinetic modules, sets of complexes with equal rows of Yᵀ B, in the order of their first complex,
             singletons included. With detached uptakes the zero complex is a module of its own, listed last.
acr          metabolites with a zero row of B, in index order
acrr_groups  groups of two or more metabolites without ACR whose rows of B are equal, each in index order
"""

RouteResult = namedtuple("RouteResult", ["route", "modules", "module_of", "n_modules", "giant_module", "acr", "acrr",
                                         "B", "dimension", "flux", "rounds", "concordance"])
RouteResult.__doc__ = """The relations one route reports, in one layout for all three routes.

route         "iterated", "tangent" or "reference"
modules       kinetic modules as sets of complex indices. They partition the complexes, singletons included. With
              detached uptakes the zero complex is coupled to nothing and is a module of its own, listed last.
module_of     numpy array, the module of each complex, counted from 1
n_modules     number of modules with two or more complexes
giant_module  size of the largest module
acr           metabolites with ACR, named by crn.metabolite_ids, in index order
acrr          sorted pairs of metabolites with ACRR, named by crn.metabolite_ids. Pairs of two ACR metabolites are
              left out, since their ratio is fixed only because each concentration is.
B             the basis the relations were read from, in row form with one row per metabolite. ker((YΔ)ᵀ) for the
              iterated route, the metabolite rows of ker(W Yᵀ) for the tangent route, None for the reference route.
dimension     the dimension of the nullspace the relations were read from, its number of basis vectors.
              dim ker((YΔ)ᵀ) = m − rank(YΔ) for the iterated route, the directions in which ln x can still move by
              what the route has proven. dim ker(W Yᵀ) for the tangent route, the directions of the steady states
              through the reference to first order. With detached uptakes they include the source coordinates, so
              this can exceed the dimension of B in metabolite space. None for the reference route.
flux          the interior flux r* the route ran at (rate constants k = r* with the steady state x = 1, see
              :func:`kimonu.rate_constants`), None for the reference route
rounds        rounds of the iterated route, or passes of the overlap merge of the reference route; None for the
              tangent route, which needs no iteration
concordance   the reference route's :class:`kimonu.ConcordanceResult`, None for the other routes
"""


def composition_differences(crn, modules, known_acr=()):
    """(YΔ)ᵀ in row form, one row per column of YΔ: the composition differences within each module.

    Δ (complexes x pairs) has one column e_c - e_d for each member c of a module other than its first member d, which
    spans the same subspace as all pairwise differences would. Metabolites already known to have ACR contribute
    their unit vector as an extra column, as in Langary et al. (2025). Once a metabolite is known to be fixed,
    differences involving it become decidable, which can license merges that were not decidable before.

    The differences are exact, from the rows of crn.Yt. A float could not hold them, since biomass coefficients come
    with 15 significant digits and a difference such as 9.86915907870449e-05 - 4.11265742135554e-05 needs more.
    """
    rows = []
    for module in modules:
        if len(module) < 2:
            continue
        first, *others = sorted(module)
        reference = dict(crn.Yt[first])
        for other in others:
            difference = dict(reference)
            for i, v in crn.Yt[other]:
                difference[i] = difference.get(i, 0) - v
            rows.append(tuple(sorted((i, v) for i, v in difference.items() if v)))
    return rows + [((met, 1),) for met in sorted(known_acr)]


def free_directions(crn, modules, known_acr=()):
    """ker((YΔ)ᵀ) in row form, one row per metabolite, the basis B of the iterated route (Eq. 4, left).

    Each module contributes the differences of its first complex against every other member as columns of YΔ, and
    each metabolite in ``known_acr`` (indices) its unit vector (:func:`composition_differences`). v lies in the span
    of YΔ exactly when it is orthogonal to ker((YΔ)ᵀ), exactly, over the rationals. With no column to span, nothing
    is fixed and the basis is the identity.
    """
    return null_space_rows(composition_differences(crn, modules, known_acr), len(crn.metabolite_ids))


def complex_directions(crn, B):
    """Yᵀ B in row form, one row per complex: how ln ψ_c can differ between steady states.

    Row c is the rows of B for complex c's metabolites, weighted by its composition. Complexes with equal rows are
    kinetically coupled, and :func:`kimonu.kinetic_flux_coupling` lifts the same rows into flux space.
    """
    return matrix_product(crn.Yt, B)


def readout(crn, B, detach_uptakes=True):
    """Kinetic modules, ACR metabolites and ACRR groups from a basis B, by grouping identical rows.

    Args:
        crn: the network (:class:`kimonu.CRN`).
        B: the basis in row form, one row per metabolite, as :func:`free_directions` returns it for the iterated
            route and :func:`kimonu.tangent_space` for the tangent route.
        detach_uptakes: True when uptakes are detached, which leaves the zero complex out of the grouping and lists
            it last as a module of its own (see the module docstring).

    Returns a :class:`Readout` in indices.
    """
    skip = zero_complexes(crn) if detach_uptakes else set()
    # 1) Modules: complexes with equal rows of Yᵀ B, whose composition differences B annihilates
    groups = {}
    for c, row in enumerate(complex_directions(crn, B)):
        if c not in skip:
            groups.setdefault(row, set()).add(c)
    # 2) ACR: zero rows of B. ACRR: equal nonzero rows. Pairs of ACR metabolites are left out of ACRR, since their
    #    ratio is constant only because each concentration is, which is not an independent finding.
    same = {}
    for i, row in enumerate(B):
        same.setdefault(row, []).append(i)
    return Readout(modules=list(groups.values()) + [{c} for c in sorted(skip)], acr=same.pop((), []),
                   acrr_groups=[g for g in same.values() if len(g) > 1])


def route_result(route, crn, found, B=None, dimension=None, flux=None, rounds=None, concordance=None):
    """The :class:`RouteResult` of a route, from what it found, a :class:`Readout` in indices.

    Metabolites are named by crn.metabolite_ids, and every ACRR group gives all its pairs.
    """
    module_of = np.zeros(len(crn.Yt), dtype=np.int64)
    for idx, group in enumerate(found.modules, start=1):
        for c in group:
            module_of[c] = idx
    name = crn.metabolite_ids
    return RouteResult(
        route=route,
        modules=found.modules,
        module_of=module_of,
        n_modules=sum(1 for g in found.modules if len(g) > 1),
        giant_module=max((len(g) for g in found.modules), default=0),
        acr=[name[i] for i in sorted(found.acr)],
        acrr=sorted((name[a], name[b]) for g in found.acrr_groups for k, a in enumerate(g) for b in g[k + 1:]),
        B=B,
        dimension=dimension,
        flux=flux,
        rounds=rounds,
        concordance=concordance,
    )
