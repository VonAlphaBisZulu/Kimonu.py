"""Balanced complexes and concordance, read off the null space of N.

Concordance is flux coupling lifted one matrix product up. With N = Y·A and K = ker(N), the rows of
M = A·K describe how each complex's activity moves over the steady states: a complex is balanced when
its row is zero, and two complexes are concordant when their rows are proportional. Everything is exact
linear algebra over the rationals, with K from straindesign's sparse null space; the only LP is the
FVA that removes blocked reactions beforehand.
"""
from collections import namedtuple

import straindesign as sd

from kimonu.io import build_crn
from kimonu.linalg import matrix_product

ConcordanceResult = namedtuple(
    "ConcordanceResult",
    ["n_complexes", "n_balanced", "n_modules", "n_concordant_pairs",
     "modules", "balanced", "complexes", "reaction_ids", "nullity"],
)


def concordance(model, preprocess=True, solver=None):
    """Compute the exact concordance partition of a cobra ``model``.

    A wrapper around :func:`concordance_from_matrices`: it removes blocked reactions (FVA), builds the
    complex decomposition N = Y·A from the model, and hands A and N on.

    Args:
        preprocess: if True (default), remove blocked reactions via FVA first (scenario i). Pass a
            model that is already FVA-unblocked with ``preprocess=False``.
        solver: solver for that FVA; None lets straindesign choose.

    Returns a :class:`ConcordanceResult`. ``modules`` are groups of mutually concordant complexes
    (proportional M-rows, singletons included); ``balanced`` are complexes whose activity is
    identically zero over the steady-state space.
    """
    if preprocess:
        # The partition is exact over the feasible steady-state set, whose linear span is ker(N) once
        # the reactions that cannot carry flux are gone. This FVA is the only LP in the analysis.
        model = model.copy()
        ranges = sd.fva(model, solver=solver) if solver else sd.fva(model)
        blocked = [r for r in ranges.index
                   if abs(ranges.at[r, "minimum"]) <= 1e-9 and abs(ranges.at[r, "maximum"]) <= 1e-9]
        if blocked:
            model.remove_reactions(blocked, remove_orphans=True)
    crn = build_crn(model)
    return concordance_from_matrices(crn.A, crn.N_exact, complexes=crn.complexes,
                                     reaction_ids=crn.reaction_ids)


def concordance_from_matrices(A, N, complexes=None, reaction_ids=None):
    """Concordance partition from the incidence matrix ``A`` and stoichiometry ``N`` directly.

    For networks that already carry a complex decomposition — the elementary-step-split models
    distributed with the reference implementation, for instance — there is nothing to enumerate,
    and rebuilding a cobra model only to take it apart again would risk changing the decomposition
    being analysed.

    ``N`` must have integer or rational entries; the kernel is computed from it exactly.
    """
    # 1) K = ker(N): every steady-state flux vector is a combination K·alpha of its columns.
    K = sd.sparse_nullspace(N)
    # 2) M = A·K: each column of M is a pattern of complex formation and consumption that some
    #    steady state produces; row c follows complex c's net formation rate across all of them.
    M = matrix_product(A, K)
    # 3) Compare rows: complexes whose rows are proportional are always formed or consumed at the
    #    same relative rate.

    balanced = []
    groups = {}                          # canonical key -> [complex indices]
    for c, row in enumerate(M):
        if not row:
            # balanced: no steady state ever net produces or consumes the complex, so its row is zero
            balanced.append(c)
            continue
        # concordant: rows proportional, of any sign. Normalize so the first nonzero == 1; identical
        # keys <=> same 1-D span <=> concordant, and grouping equal keys is transitive by construction
        nz = sorted(row.items())
        key = tuple((col, v / nz[0][1]) for col, v in nz)
        groups.setdefault(key, []).append(c)

    modules = list(groups.values())
    pairs = sum(len(g) * (len(g) - 1) // 2 for g in modules)

    return ConcordanceResult(
        n_complexes=A.shape[0],
        n_balanced=len(balanced),
        n_modules=len(modules),
        n_concordant_pairs=pairs,
        modules=modules,
        balanced=balanced,
        complexes=complexes if complexes is not None else list(range(A.shape[0])),
        reaction_ids=reaction_ids if reaction_ids is not None else list(range(A.shape[1])),
        nullity=K.shape[1],
    )
