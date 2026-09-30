"""Exact concordance partition via the sparse exact nullspace.

Concordance is flux coupling lifted one matrix product up: with N = Y·A and K = ker(N), two
complexes are concordant iff the corresponding rows of M = A·K are proportional, and balanced iff
the M-row is zero. Pure exact linear algebra over the rationals — no LP, no sampling. K comes from
straindesign. Ported from the validated Julia proof-of-concept (concordance_nullspace.jl).
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

    Args:
        preprocess: if True (default), remove blocked reactions via FVA first (scenario i). Pass a
            model that is already FVA-unblocked with ``preprocess=False``.
        solver: solver for that FVA; None lets straindesign choose.

    Returns a :class:`ConcordanceResult`. ``modules`` are groups of mutually concordant complexes
    (proportional M-rows, singletons included, matching the Julia PoC); ``balanced`` are complexes
    whose activity is identically zero over the steady-state space.
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
    # Every steady-state flux is K·alpha, so row c of M = A·K says how complex c's activity -- its net
    # formation rate -- changes along each steady-state direction.
    K = sd.sparse_nullspace(N)
    M = matrix_product(A, K)

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
