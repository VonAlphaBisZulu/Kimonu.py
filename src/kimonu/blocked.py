"""Blocked-reaction detection without linear programming.

A reaction is blocked when it carries zero flux in every steady state. Kimonu needs to know which
reactions those are before it can use the kernel test at all: "activity zero over ker(N)" agrees
with "activity zero over the feasible set" exactly when no reaction is blocked, so the FVA that
finds them is the one LP in an otherwise LP-free pipeline.

Two causes, and they are not equally hard.

**Stoichiometrically blocked.** No steady-state flux vector at all, of either sign, gives the
reaction nonzero flux: ``v_j = 0`` for every ``v`` in ker(N). This is exactly a zero row of a basis
K of ker(N) — the same primitive the concordance core already uses, and free if K is already at
hand.

**Irreversibility blocked.** The reaction could carry flux in ker(N), but not once the sign
constraints are imposed. Deciding this exactly is an implicit-equality question for the polyhedral
cone ``{v ∈ ker(N) : v_j ≥ 0, j irreversible}``, and that is what an LP answers. What is available
without one is propagation: a metabolite that can only be produced, and never consumed, forces its
producers to zero, which removes reactions, which can strand further metabolites. Iterating that to
a fixed point is exact for what it finds but is not guaranteed to find everything.

So the honest split is: the stoichiometric half is exactly and cheaply decidable, the irreversible
half is decidable only in the direction of "certainly blocked". :func:`blocked_reactions` returns
both, and says which reactions are certain and whether anything might remain.

Finite bounds do not enter into it. Scaling a feasible flux vector towards zero stays feasible, so
blockedness depends only on which bounds are zero and which are signed — not on their magnitude.
That is what makes the normalised form of the problem the right one to reason about.
"""
from collections import namedtuple

import numpy as np
from scipy import sparse

import straindesign as sd

BlockedResult = namedtuple(
    "BlockedResult",
    ["blocked", "stoichiometric", "irreversibility", "complete", "iterations"],
)


def stoichiometrically_blocked(N, K=None):
    """Reactions with zero flux throughout ker(N) — the zero rows of a kernel basis."""
    if K is None:
        K = sd.sparse_nullspace(N)
    if isinstance(K, sd.ExactCOO):
        touched = set(int(r) for r in K.rows)
        return {j for j in range(K.shape[0]) if j not in touched}
    K = K.tocsr()
    return {j for j in range(K.shape[0]) if K.indptr[j] == K.indptr[j + 1]}


def _propagate_dead_ends(N, reversible, orientation, alive):
    """Remove reactions forced to zero because a metabolite they touch cannot balance.

    With only irreversible reactions left touching a metabolite in one direction, steady state
    forces every one of them to zero: nothing can consume what they make. Removing them can strand
    the next metabolite, so this runs to a fixed point.
    """
    N = N.tocsr()
    iterations = 0
    while True:
        iterations += 1
        doomed = set()
        for i in range(N.shape[0]):
            start, end = N.indptr[i], N.indptr[i + 1]
            producers, consumers, reversible_here = [], [], False
            for j, coeff in zip(N.indices[start:end], N.data[start:end]):
                j = int(j)
                if j not in alive:
                    continue
                if reversible[j]:
                    reversible_here = True
                    break
                # A reaction confined to negative flux runs backwards, so the sign of its
                # stoichiometric coefficient is the opposite of what it does to the metabolite.
                effective = coeff * orientation[j]
                (producers if effective > 0 else consumers).append(j)
            if reversible_here:
                continue
            if producers and not consumers:
                doomed.update(producers)
            elif consumers and not producers:
                doomed.update(consumers)
        if not doomed:
            return alive, iterations
        alive -= doomed


def _propagate_signs(N, lb, ub, orientation):
    """Narrow each reaction's sign domain using the metabolite balances, to a fixed point.

    Every metabolite balance sums to zero. If, given what is currently known, every term in one
    balance is non-negative, then every term is zero — nothing can offset them. And if all terms
    but one are non-negative, that one is non-positive, which narrows its domain and may close a
    different balance on the next pass.

    That second step is what the dead-end rule cannot do. In an anaerobic model the oxygen exchange
    can only remove o2_e, so the reversible transporter must run backwards to supply it; the o2_c
    balance then requires it to run forwards; only zero satisfies both. Neither balance decides
    anything alone, and the transporter is reversible, so a rule that skips metabolites touched by
    reversible reactions never sees it.
    """
    N = N.tocsr()
    n = N.shape[1]
    # sign domain per reaction: whether negative and positive flux remain possible
    neg = [bool(lb[j] < 0) for j in range(n)]
    pos = [bool(ub[j] > 0) for j in range(n)]

    changed = True
    iterations = 0
    while changed:
        changed = False
        iterations += 1
        for i in range(N.shape[0]):
            start, end = N.indptr[i], N.indptr[i + 1]
            terms = [(int(j), float(c)) for j, c in zip(N.indices[start:end], N.data[start:end])
                     if neg[int(j)] or pos[int(j)]]
            if not terms:
                continue

            # sign range each term can contribute
            def term_can(j, coeff, want_positive):
                if coeff > 0:
                    return pos[j] if want_positive else neg[j]
                return neg[j] if want_positive else pos[j]

            all_non_negative = not any(term_can(j, c, False) for j, c in terms)
            all_non_positive = not any(term_can(j, c, True) for j, c in terms)
            if all_non_negative or all_non_positive:
                for j, _ in terms:
                    if neg[j] or pos[j]:
                        neg[j] = pos[j] = False
                        changed = True
                continue

            # all but one non-negative forces the remaining one non-positive, and vice versa
            for j, coeff in terms:
                others = [(k, c) for k, c in terms if k != j]
                if not any(term_can(k, c, False) for k, c in others):
                    if coeff > 0 and pos[j]:
                        pos[j] = False
                        changed = True
                    elif coeff < 0 and neg[j]:
                        neg[j] = False
                        changed = True
                if not any(term_can(k, c, True) for k, c in others):
                    if coeff > 0 and neg[j]:
                        neg[j] = False
                        changed = True
                    elif coeff < 0 and pos[j]:
                        pos[j] = False
                        changed = True

    return {j for j in range(n) if not neg[j] and not pos[j]}, iterations


def blocked_reactions(N, lb, ub, K=None):
    """Reactions that certainly carry zero flux in every steady state, found without an LP.

    Args:
        N: metabolites x reactions stoichiometry, integer or rational.
        lb, ub: flux bounds; only their signs are consulted.
        K: a basis of ker(N), if already computed.

    Returns a :class:`BlockedResult`. ``complete`` is True when the answer is provably the whole
    set — which is the case when every reaction is reversible, since then there are no sign
    constraints for the cone to hide an implicit equality behind.
    """
    n = N.shape[1]
    lb = np.asarray(lb, dtype=float).ravel()
    ub = np.asarray(ub, dtype=float).ravel()
    reversible = [bool(lb[j] < 0 < ub[j]) for j in range(n)]
    orientation = [(-1.0 if (ub[j] <= 0 and lb[j] < 0) else 1.0) for j in range(n)]

    fixed_zero = {j for j in range(n) if lb[j] == 0 and ub[j] == 0}
    structural = stoichiometrically_blocked(N, K) | fixed_zero

    alive = set(range(n)) - structural
    alive, iterations = _propagate_dead_ends(N, reversible, orientation, alive)
    by_signs, sign_iterations = _propagate_signs(N, lb, ub, orientation)
    iterations += sign_iterations
    irreversibility = (set(range(n)) - structural - alive) | (by_signs - structural)

    complete = all(reversible[j] for j in range(n) if j not in structural)
    return BlockedResult(
        blocked=sorted(structural | irreversibility),
        stoichiometric=sorted(structural),
        irreversibility=sorted(irreversibility),
        complete=complete,
        iterations=iterations,
    )
