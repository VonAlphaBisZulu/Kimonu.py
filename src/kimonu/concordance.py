"""Balanced complexes and concordance, read off the null space of N.

Concordance is flux coupling lifted one matrix product up. With N = Y·A and K = ker(N), the rows of
M = A·K describe how each complex's activity moves over the steady states: a complex is balanced when
its row is zero, and two complexes are concordant when their rows are proportional. Everything is exact
linear algebra over the rationals, with K from straindesign's sparse null space; the only LP is the
one beforehand that finds blocked reactions. FVA is needed only when lower bounds force a flux, since
only then can a bound pin a rate.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np
import straindesign as sd
from cobra.util.array import create_stoichiometric_matrix

from kimonu.io import build_crn, irreversible
from kimonu.linalg import blocked_reactions, entries, exact, matrix_product

ConcordanceResult = namedtuple(
    "ConcordanceResult",
    ["n_complexes", "n_balanced", "n_modules", "n_concordant_pairs",
     "modules", "balanced", "complexes", "reaction_ids", "nullity"],
)


def concordance(model, preprocess=True, solver=None):
    """Compute the exact concordance partition of a cobra ``model``.

    A wrapper around :func:`concordance_from_matrices`: it splits reversible reactions into two
    (:func:`kimonu.io.irreversible`), removes blocked reactions and finds pinned rates (:func:`_preprocess`),
    builds the complex decomposition N = Y·A from the model, and hands A and N on.

    Args:
        preprocess: if True (default), do that first. Pass a model that has no blocked reactions and
            no pinned rates with ``preprocess=False``.
        solver: solver for FVA, if the model needs it; None lets straindesign choose.

    Returns a :class:`ConcordanceResult`. ``modules`` are groups of mutually concordant complexes
    (proportional M-rows, singletons included); ``balanced`` are complexes whose activity is
    identically zero over the steady-state space.
    """
    pins = {}
    if preprocess:
        # The partition is exact over the feasible steady states P. Their linear span is ker(N) unless
        # some rate is constant across P: at zero, a blocked reaction, which is removed; or pinned at a
        # nonzero value. Two independent pins couple reactions that ker(N) does not, so the pins enter
        # the kernel as extra rows of N (_with_pins).
        model = irreversible(model).copy()
        blocked, pins = _preprocess(model, solver)
        if blocked:
            model.remove_reactions(blocked, remove_orphans=True)
    crn = build_crn(model)
    return concordance_from_matrices(crn.A, _with_pins(crn.N_exact, crn.reaction_ids, pins),
                                     complexes=crn.complexes, reaction_ids=crn.reaction_ids)


def _preprocess(model, solver=None):
    """(ids of the blocked reactions, pinned rates) of an irreversible model.

    With every lower bound at zero, r = 0 is a steady state, so a rate that is constant over the steady
    states is zero: nothing is pinned, and finite upper bounds only scale, so they cannot block a reaction.
    One LP then settles it. A positive lower bound forces flux and needs FVA.
    """
    if all(r.lower_bound <= 0 for r in model.reactions):
        S = create_stoichiometric_matrix(model, array_type="lil")
        closed = [j for j, r in enumerate(model.reactions) if r.upper_bound <= 0]
        return [model.reactions[j].id for j in np.flatnonzero(blocked_reactions(S.tocsr(), closed))], {}
    ranges = sd.fva(model, solver=solver) if solver else sd.fva(model)
    blocked = [r for r in ranges.index
               if abs(ranges.at[r, "minimum"]) <= 1e-9 and abs(ranges.at[r, "maximum"]) <= 1e-9]
    return blocked, _pins(model, ranges)


def _pins(model, ranges, tol=1e-9):
    """{reaction id: exact value} for the rates that a bound holds at a nonzero value in every steady state.

    FVA shows such a rate as minimum == maximum. The bound it sits on gives the value exactly; a rate
    that is constant only because others are pinned follows from those and needs no entry.
    """
    pins = {}
    for r in model.reactions:
        lo, hi = ranges.at[r.id, "minimum"], ranges.at[r.id, "maximum"]
        if abs(hi) <= tol or hi - lo > tol * max(1, abs(hi)):
            continue
        for b in (r.lower_bound, r.upper_bound):
            if abs(hi - b) <= tol * max(1, abs(b)):
                pins[r.id] = Fraction(str(b))
                break
    return pins


def _with_pins(N, reaction_ids, pins):
    """N with one row c_l·r_k − c_k·r_l = 0 for every pinned rate l after the first one, k.

    Rates pinned at c_k and c_l keep that ratio along the whole linear span of the feasible steady
    states, and that span, not ker(N), is what the kernel has to be for the readout to be complete.
    One pin or none adds nothing.
    """
    if len(pins) < 2:
        return N
    triples, (m, q) = entries(N)
    col = {r: j for j, r in enumerate(reaction_ids)}
    (k, ck), *rest = pins.items()
    for i, (l, cl) in enumerate(rest):
        triples += [(m + i, col[k], cl), (m + i, col[l], -ck)]
    return exact(triples, (m + len(rest), q))


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
