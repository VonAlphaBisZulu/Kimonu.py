"""Flux balance: balanced complexes and concordance, read off the nullspace K of N (Fig. 1C, E).

Every steady-state flux is a combination r = K·α of the columns of a basis K of ker(N). When rows i and j of K
are proportional, row j = λ · row i, reactions i and j are fully coupled, r_j = λ r_i at every steady state (Burgard
et al., 2004, :func:`kimonu.mass_balance_flux_coupling`). Concordance is flux coupling lifted one matrix product
up. With N = Y·A, the rows of M = A·K describe how each complex's activity, its net formation rate a = A·r, moves
over the steady states. A complex is balanced when its row is zero, and two complexes are concordant when their rows
are proportional (Küken et al., 2022), so that a_c = λ a_d at every steady state. Everything is exact linear algebra
over the rationals, with K straindesign's RREF basis of ker(N) (:func:`kimonu.linalg.nullspace`), and both readouts
are :func:`kimonu.readout.row_readout`, of K for fully coupled reactions and of M for concordant complexes. The only
LP is the one beforehand that finds blocked reactions (:func:`kimonu.preprocess`).

Pinned rates. Concordance is exact over the feasible steady states P. Their linear span is ker(N) unless some rate is
constant across P, at zero for a blocked reaction, which :func:`kimonu.preprocess` removes, or pinned by a bound at
a nonzero value. Rates pinned at c_k and c_l keep that ratio along the whole span of P, which ker(N) does not show,
so two pins couple reactions that ker(N) does not. Each pinned rate l after the first one, k, therefore adds the
homogeneous row c_l·r_k − c_k·r_l = 0 to N. A single pin adds no row, since it fixes only the scale of a flux and
concordance compares directions.
"""
from collections import namedtuple

from kimonu.linalg import entries, exact, matrix_product, nullspace
from kimonu.network import as_crn
from kimonu.readout import row_readout

ConcordanceResult = namedtuple(
    "ConcordanceResult",
    ["n_complexes", "n_balanced", "n_modules", "n_concordant_pairs",
     "modules", "activity_ratio", "balanced", "complexes", "reaction_ids", "nullity", "K", "M"],
)
ConcordanceResult.__doc__ = """Balanced complexes and concordance modules, with the matrices they were read from.

modules             concordance modules, lists of mutually concordant complexes (proportional rows of M), in the order
                    of their first complex, singletons included. Balanced complexes are in none of them.
activity_ratio      {complex: λ} for every complex that is not balanced, with a_c = λ a_d at every steady state for
                    the first complex d of its module (row c of M = λ · row d), so λ = 1 for the first complex
balanced            complexes whose row of M is zero, so no steady state net produces or consumes them
n_modules           number of modules, singletons included
n_concordant_pairs  pairs of complexes within one module
complexes, reaction_ids  the names of the network (crn.complexes, crn.reaction_ids)
nullity             dim ker(N), the number of columns of K
K                   the RREF basis of ker(N) in row form, one row per reaction. Proportional rows are fully coupled
                    reactions (:func:`kimonu.mass_balance_flux_coupling`).
M                   A·K in row form, one row per complex
"""


def concordance(network, preprocess=True, solver=None):
    """The exact concordance partition of a network: balanced complexes and concordance modules.

    Args:
        network: a :class:`kimonu.CRN`, or a cobra model, which is preprocessed first (:func:`kimonu.as_crn`). Pass
            the CRN when calling several functions, so that preprocessing runs once.
        preprocess: for a cobra model, run :func:`kimonu.preprocess` (default), which splits reversible reactions,
            removes blocked reactions and finds pinned rates. A model that has no blocked reactions and no pinned
            rates can be passed with ``preprocess=False``, which only builds N = Y·A (:func:`kimonu.build_crn`).
        solver: solver for FVA, if the model needs it. None lets straindesign choose.

    The pins of the CRN enter N as rows that fix the ratio of every two pinned rates (see the module docstring).

    Returns a :class:`ConcordanceResult`.
    """
    crn = as_crn(network, preprocess, solver)
    # 1) K = ker(N): every steady-state flux vector is a combination K·α of its columns.
    K = nullspace(_with_pins(crn))
    # 2) M = A·K: each column of M is a pattern of complex formation and consumption that some steady state
    #    produces, and row c follows complex c's net formation rate across all of them.
    M = matrix_product(crn.A, K)
    # 3) Compare rows: a zero row is a balanced complex, and complexes whose rows are proportional, of any sign, are
    #    always formed or consumed at the same relative rate, the ratio of their rows.
    balanced, modules, activity_ratio = row_readout(M)
    return ConcordanceResult(
        n_complexes=len(M),
        n_balanced=len(balanced),
        n_modules=len(modules),
        n_concordant_pairs=sum(len(g) * (len(g) - 1) // 2 for g in modules),
        modules=modules,
        activity_ratio=activity_ratio,
        balanced=balanced,
        complexes=crn.complexes,
        reaction_ids=crn.reaction_ids,
        nullity=len({j for row in K for j, _ in row}),
        K=K,
        M=M,
    )


def _with_pins(crn):
    """N with the homogeneous row c_l·r_k − c_k·r_l = 0 for every pinned rate l after the first one, k.

    One pin or none adds no row and returns N itself (see the module docstring).
    """
    if len(crn.pins) < 2:
        return crn.N
    triples, (m, q) = entries(crn.N)
    col = {r: j for j, r in enumerate(crn.reaction_ids)}
    (k, ck), *rest = crn.pins.items()
    for i, (l, cl) in enumerate(rest):
        triples += [(m + i, col[k], cl), (m + i, col[l], -ck)]
    return exact(triples, (m + len(rest), q))
