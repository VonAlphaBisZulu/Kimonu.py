"""Exact concordance partition via the sparse exact nullspace.

Concordance is flux coupling lifted one matrix product up: with N = Y·A and K = ker(N), two
complexes are concordant iff the corresponding rows of M = A·K are proportional, and balanced iff
the M-row is zero. Pure exact linear algebra over the rationals — no LP, no sampling. K comes from
SENUS. Ported from the validated Julia proof-of-concept (concordance_nullspace.jl).
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np
from scipy import sparse

import straindesign as sd

from kimonu.io import build_crn

ConcordanceResult = namedtuple(
    "ConcordanceResult",
    ["n_complexes", "n_balanced", "n_modules", "n_concordant_pairs",
     "modules", "balanced", "complexes", "reaction_ids"],
)


def _M_rows(A, K):
    """Yield (row_index, [(col, Fraction value), ...]) for the nonzero rows of M = A·K.

    A is integer (complexes x reactions). K is the exact nullspace (reactions x nu), either a scipy
    integer CSR (common word-sized case) or a sd.ExactCOO (>int64 coefficients). Either way M is
    exact; a common denominator only scales rows and never affects proportionality, so it is dropped.
    """
    if isinstance(K, sd.ExactCOO):
        # K[r, c] = data / denom; build M = A·K exactly with Python ints (denom cancels for ratios).
        Kcols = {}
        for r, c, v in zip(K.rows, K.cols, K.data):
            Kcols.setdefault(r, []).append((c, v))
        Acsr = A.tocsr()
        for i in range(Acsr.shape[0]):
            row = Acsr.getrow(i)
            acc = {}
            for rxn, a in zip(row.indices, row.data):
                for c, v in Kcols.get(int(rxn), ()):
                    acc[c] = acc.get(c, 0) + int(a) * int(v)
            nz = sorted((c, Fraction(val)) for c, val in acc.items() if val != 0)
            yield i, nz
    else:
        M = (A @ K).tocsr()   # complexes x nu, integer
        for i in range(M.shape[0]):
            row = M.getrow(i)
            if row.nnz == 0:
                yield i, []
                continue
            order = np.argsort(row.indices)
            nz = [(int(row.indices[k]), Fraction(int(row.data[k]))) for k in order]
            yield i, nz


def _fva_unblocked(model):
    """Return a copy of ``model`` with blocked reactions removed (FVA-based).

    The concordance partition is exact over the *feasible* steady-state set, whose linear span is
    ker(N) restricted to reactions that can carry flux. This preprocessing (an LP-based FVA, done
    once) is the only place an LP appears — the concordance detection itself is LP-free. Matches
    COCOA's scenario (i) (irreversibility only).
    """
    from cobra.flux_analysis import find_blocked_reactions
    m = model.copy()
    blocked = find_blocked_reactions(m)
    if blocked:
        m.remove_reactions(blocked, remove_orphans=True)
    return m


def concordance(model, preprocess=True):
    """Compute the exact concordance partition of a cobra ``model``.

    Args:
        preprocess: if True (default), remove blocked reactions via FVA first (scenario i). Pass a
            model that is already FVA-unblocked with ``preprocess=False``.

    Returns a :class:`ConcordanceResult`. ``modules`` are groups of mutually concordant complexes
    (proportional M-rows, singletons included, matching the Julia PoC); ``balanced`` are complexes
    whose activity is identically zero over the steady-state space.
    """
    if preprocess:
        model = _fva_unblocked(model)
    crn = build_crn(model)
    K = sd.sparse_nullspace(crn.N)   # reactions x nu, exact

    balanced = []
    groups = {}                          # canonical key -> [complex indices]
    for i, nz in _M_rows(crn.A, K):
        if not nz:
            balanced.append(i)
            continue
        c0, v0 = nz[0]
        # normalize so the first nonzero == 1; identical keys <=> same 1-D span <=> concordant
        key = tuple((c, val / v0) for c, val in nz)
        groups.setdefault(key, []).append(i)

    modules = list(groups.values())
    pairs = sum(len(g) * (len(g) - 1) // 2 for g in modules)

    return ConcordanceResult(
        n_complexes=len(crn.complexes),
        n_balanced=len(balanced),
        n_modules=len(modules),
        n_concordant_pairs=pairs,
        modules=modules,
        balanced=balanced,
        complexes=crn.complexes,
        reaction_ids=crn.reaction_ids,
    )
