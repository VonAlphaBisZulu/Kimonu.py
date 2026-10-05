"""Kinetic modules from a null space in complex-concentration space.

Under mass action, reaction j runs at r_j = k_j * psi_s(j), where psi_c is the monomial of complex c
and s(j) the substrate complex of j. Any positive steady state (r*, psi*) fixes rate constants
k_j = r*_j / psi*_s(j), and every steady state for those constants satisfies

    W u = 0,   W = N diag(r*) S^T,   u = psi / psi*,

with S the substrate part of the incidence matrix. Complexes whose rows in a basis of ker(W) are
proportional keep a fixed psi ratio at every positive steady state: they are kinetically coupled.
Every rate-constant set that admits a positive steady state arises this way, and r* is drawn at
random from the interior of the steady-state cone, so the coupling holds for generic rate constants.

The coupled sets are then closed in log space. Within a set, ln psi_c - ln psi_d =
(y_c - y_d)^T ln x is constant, so any composition difference in the span of those is constant as
well, and sets whose representatives differ by such a vector merge (Proposition S3-4 of Langary et
al., 2025). Metabolites found to have ACR on the way extend that span (Remark S3-6).

A merge found in log space is a new linear fact in psi space: u_c = u_d for the merged complexes,
since their ratio is fixed and equals one at u = 1. Added to W, it can make further rows of the
kernel proportional, because the proportional-row readout alone ignores linear relations among
three or more coupled sets. On EnvZ-OmpR this feedback is what couples the complexes that differ
by Yp. Both steps are iterated until the partition no longer changes.

Complexes that are no reaction's substrate have a zero column in W, so their u is free in ker(W) and
they reach the closure as singletons. Their psi is still a monomial of x, and the closure joins them to a
module exactly when their composition differs from a member's by a fixed log-combination; this is the
coupling of Langary et al. (2025, Supplement S.1), which is not restricted to substrates. Joining adds no
new vector to the span, so it changes no ACR or ACRR, and the u_c = u_d it feeds back only fixes u_c.

The zero complex has psi = 1. Attached (default), its uptake reactions run at constant rate and
complexes coupled to it share its module: a fixed medium. Detached, each uptake reaction is given
its own source, so only couplings that hold whatever the uptake rates remain.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np

from kimonu.acr import acr_acrr, coupling_signatures
from kimonu.io import substrate_complexes, zero_complexes
from kimonu.linalg import entries, exact, null_space_rows, positive_flux, proportional_classes

KineticResult = namedtuple(
    "KineticResult",
    ["modules", "module_of", "n_modules", "giant_module", "known_acr", "rounds"],
)


def _psi_coupled(A, Y, N, flux, detach_uptakes, equal=()):
    """Complexes whose psi ratios are fixed by the steady-state equations W u = 0, W = N diag(r*) S^T.

    ``equal`` lists pairs (c, d) already known to satisfy u_c = u_d. Returns the coupled sets and the
    substrate complexes (the only complexes that appear in a rate law).
    """
    n = A.shape[0]
    substrate, zero = substrate_complexes(A), zero_complexes(Y)
    by_reaction = {}
    triples, (m, _) = entries(N)
    for i, j, v in triples:
        by_reaction.setdefault(j, []).append((i, v))
    # W: column c collects, for every reaction consuming complex c, its stoichiometry times its rate at
    # r*. A detached uptake gets a column of its own, a free source rather than the constant psi_0 = 1.
    W, extra = {}, {}
    for j, c in substrate.items():
        col = extra.setdefault(j, n + len(extra)) if detach_uptakes and c in zero else c
        for i, v in by_reaction.get(j, ()):
            W[(i, col)] = W.get((i, col), Fraction(0)) + flux[j] * v
    for k, (c, d) in enumerate(equal):                    # u_c - u_d = 0
        W[(m + k, c)], W[(m + k, d)] = Fraction(1), Fraction(-1)
    W = exact(((i, j, v) for (i, j), v in W.items()), (m + len(equal), n + len(extra)))
    # proportional rows of ker(W): these complexes keep fixed psi ratios at every steady state
    _, classes = proportional_classes(null_space_rows(W))
    drop = zero if detach_uptakes else set()
    return [{c for c in g if c < n and c not in drop} for g in classes], set(substrate.values()) - drop


def _closure(Y, sets):
    """Merge sets with equal signatures (Proposition S3-4), extended by ACR metabolites (Remark S3-6),
    until neither the sets nor the ACR metabolites change."""
    known, rounds = set(), 0
    while True:
        rounds += 1
        signatures = coupling_signatures(Y, sets, known)
        if signatures is None:
            return sets, known, rounds
        merged = {}
        for s, key in zip(sets, signatures):
            merged.setdefault(key, set()).update(s)
        merged = list(merged.values())
        found = known | set(acr_acrr(Y, merged, known_acr=known).acr)
        if len(merged) == len(sets) and found == known:
            return merged, known, rounds
        sets, known = merged, found


def kinetic_modules(A, Y, N, detach_uptakes=False, seed=0):
    """Kinetic modules of a network with no blocked reactions.

    Args:
        A: complexes x reactions incidence matrix, one -1 and one +1 per column.
        Y: metabolites x complexes composition matrix.
        N: metabolites x reactions stoichiometric matrix, exact (straindesign RationalMatrix) or
            integer; N = Y A.
        detach_uptakes: give every reaction out of the zero complex its own source, so couplings
            through the constant zero-complex monomial are not used.
        seed: seed of the random interior flux. Results hold for generic rate constants, so the
            seed changes them only with probability zero.

    Returns a :class:`KineticResult`; every complex belongs to exactly one module, singletons
    included. Complexes that are no reaction's substrate join a module when the log closure couples
    them to it.
    """
    n = A.shape[0]
    # 1) Generic rate constants: an interior steady-state flux r* fixes k_j = r*_j at psi* = 1.
    flux = positive_flux(N, seed)
    equal, rounds = [], 0
    while True:
        # 2) Linear part: every steady state for these k solves W u = 0; proportional rows of ker(W)
        #    are complexes with fixed psi ratios. Complexes on their own, substrates or not, still take
        #    part in the closure, which can attach them to a set.
        coupled, _ = _psi_coupled(A, Y, N, flux, detach_uptakes, equal)
        coupled = [s for s in coupled if s]
        # 3) Log part: within a set, differences of log compositions are fixed; merge sets whose
        #    representatives differ by such a vector, and fix metabolites that become ACR on the way.
        sets, known, r = _closure(Y, coupled)
        rounds += r
        # 4) Feedback: a merge found in log space is the linear fact u_c = u_d; add it to W and repeat
        #    until nothing merges.
        if len(sets) == len(coupled):
            break
        equal = [(min(s), c) for s in sets for c in s if c != min(s)]

    result = [s for s in sets if s]
    assigned = set().union(*result) if result else set()
    result += [{i} for i in range(n) if i not in assigned]

    module_of = np.zeros(n, dtype=np.int64)
    for idx, group in enumerate(result, start=1):
        for i in group:
            module_of[i] = idx

    return KineticResult(
        modules=result,
        module_of=module_of,
        n_modules=sum(1 for m in result if len(m) > 1),
        giant_module=max((len(m) for m in result), default=0),
        known_acr=sorted(known),
        rounds=rounds,
    )
