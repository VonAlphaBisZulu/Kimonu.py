"""Flux balance: balanced complexes and concordance, read off the nullspace K of N (Fig. 1C, E).

Every steady-state flux is a combination r = K·α of the columns of a basis K of ker(N). When rows i and j of K
are proportional, reactions i and j are fully coupled (Burgard et al., 2004). Concordance is flux coupling lifted
one matrix product up. With N = Y·A, the rows of M = A·K describe how each complex's net formation rate a = A·r moves
over the steady states. A complex is balanced when its row is zero, and two complexes are concordant when their rows
are proportional (Küken et al., 2022). Everything is exact linear algebra over the rationals, with K from
straindesign's sparse nullspace (the RREF basis, with 1 at each free reaction and 0 at the other free reactions).
The only LP is the one beforehand that finds blocked reactions (:func:`kimonu.prepare`).
"""
from collections import namedtuple

import straindesign as sd

from kimonu.io import build_crn, prepare
from kimonu.linalg import entries, exact, matrix_product

ConcordanceResult = namedtuple(
    "ConcordanceResult",
    ["n_complexes", "n_balanced", "n_modules", "n_concordant_pairs",
     "modules", "balanced", "complexes", "reaction_ids", "nullity", "K", "M"],
)
ConcordanceResult.__doc__ = """Balanced complexes and concordance modules, with the matrices they were read from.

modules             concordance modules, lists of mutually concordant complexes (proportional rows of M),
                    singletons included. Balanced complexes are in none of them.
balanced            complexes whose row of M is zero, so no steady state net produces or consumes them
n_modules           number of modules, singletons included
n_concordant_pairs  pairs of complexes within one module
complexes, reaction_ids  names of the rows of A and of the reactions, or indices
nullity             dim ker(N), the number of columns of K
K                   the basis of ker(N) in row form, one {basis vector: value} dict per reaction, the entries scaled
                    to integers per basis vector. Proportional rows are fully coupled reactions
                    (:func:`kimonu.linalg.proportional_classes`).
M                   A·K in row form, one {basis vector: value} dict per complex
"""


def concordance(model, preprocess=True, solver=None):
    """Compute the exact concordance partition of a cobra ``model``.

    A wrapper around :func:`concordance_from_matrices`. It runs :func:`kimonu.prepare` (split reversible reactions,
    remove blocked reactions, find pinned rates, build N = Y·A) and hands A, N and the pins on.

    Args:
        preprocess: if True (default), run :func:`kimonu.prepare` first. A model that has no blocked reactions and
            no pinned rates can be passed with ``preprocess=False``, which only builds N = Y·A.
        solver: solver for FVA, if the model needs it. None lets straindesign choose.

    Returns a :class:`ConcordanceResult`.
    """
    crn = prepare(model, solver) if preprocess else build_crn(model)
    return concordance_from_matrices(crn.A, crn.N, complexes=crn.complexes, reaction_ids=crn.reaction_ids,
                                     pins=crn.pins)


def _with_pins(N, reaction_ids, pins):
    """N with one row c_l·r_k − c_k·r_l = 0 for every pinned rate l after the first one, k.

    The partition is exact over the feasible steady states P. Their linear span is ker(N) unless some rate is
    constant across P, at zero for a blocked reaction, which :func:`kimonu.prepare` removes, or pinned at a nonzero
    value. Rates pinned at c_k and c_l keep that ratio along the whole span of P, which ker(N) does not show, so two
    independent pins couple reactions that ker(N) does not. One pin or none adds nothing.
    """
    if len(pins) < 2:
        return N
    triples, (m, q) = entries(N)
    col = {r: j for j, r in enumerate(reaction_ids)}
    (k, ck), *rest = pins.items()
    for i, (l, cl) in enumerate(rest):
        triples += [(m + i, col[k], cl), (m + i, col[l], -ck)]
    return exact(triples, (m + len(rest), q))


def concordance_from_matrices(A, N, complexes=None, reaction_ids=None, pins=None):
    """Concordance partition from the incidence matrix ``A`` and stoichiometry ``N`` directly.

    For networks that already carry a complex decomposition, such as the elementary-step models distributed with the
    reference implementation of Langary et al. (2025), there is nothing to enumerate, and rebuilding a cobra model
    only to take it apart again would risk changing the decomposition being analysed.

    Args:
        A: complexes x reactions incidence matrix.
        N: metabolites x reactions stoichiometric matrix with integer or rational entries. The kernel is computed
            from it exactly.
        complexes, reaction_ids: optional names, returned in the result.
        pins: optional {reaction: exact value} of pinned rates, keyed like ``reaction_ids`` (or by index), as
            :func:`kimonu.prepare` finds them. They enter N as rows.

    Returns a :class:`ConcordanceResult`.
    """
    if pins:
        N = _with_pins(N, reaction_ids if reaction_ids is not None else list(range(A.shape[1])), pins)
    # 1) K = ker(N): every steady-state flux vector is a combination K·α of its columns.
    K_matrix = sd.sparse_nullspace(N)
    K = [{} for _ in range(K_matrix.shape[0])]
    for i, j, v in entries(K_matrix)[0]:
        K[i][j] = v
    # 2) M = A·K: each column of M is a pattern of complex formation and consumption that some steady state
    #    produces, and row c follows complex c's net formation rate across all of them.
    M = matrix_product(A, K)
    # 3) Compare rows: complexes whose rows are proportional are always formed or consumed at the same relative rate.
    balanced = []
    groups = {}                          # canonical key -> [complex indices]
    for c, row in enumerate(M):
        if not row:
            # balanced: no steady state ever net produces or consumes the complex, so its row is zero
            balanced.append(c)
            continue
        # concordant: rows proportional, of any sign. Normalized so the first nonzero is 1, identical keys mean the
        # same 1-D span, and grouping equal keys is transitive by construction.
        nz = sorted(row.items())
        key = tuple((col, v / nz[0][1]) for col, v in nz)
        groups.setdefault(key, []).append(c)

    modules = list(groups.values())
    n_pairs = sum(len(g) * (len(g) - 1) // 2 for g in modules)

    return ConcordanceResult(
        n_complexes=A.shape[0],
        n_balanced=len(balanced),
        n_modules=len(modules),
        n_concordant_pairs=n_pairs,
        modules=modules,
        balanced=balanced,
        complexes=complexes if complexes is not None else list(range(A.shape[0])),
        reaction_ids=reaction_ids if reaction_ids is not None else list(range(A.shape[1])),
        nullity=K_matrix.shape[1],
        K=K,
        M=M,
    )
