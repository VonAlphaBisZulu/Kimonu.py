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

The zero complex has psi = 1. Attached (default), its uptake reactions run at constant rate and
complexes coupled to it share its module: a fixed medium. Detached, each uptake reaction is given
its own source, so only couplings that hold whatever the uptake rates remain.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np
import straindesign as sd
from scipy import sparse
from scipy.optimize import linprog

from kimonu.acr import acr_acrr, coupling_signatures
from kimonu.concordance import concordance_from_matrices
from kimonu.linalg import entries

KineticResult = namedtuple(
    "KineticResult",
    ["modules", "module_of", "n_modules", "giant_module", "known_acr", "rounds"],
)


def positive_flux(N, seed=0):
    """An exact steady-state flux with every entry strictly positive.

    An LP point with random lower bounds r >= 1 + U(0, 1) is snapped onto the exact kernel through
    its free coordinates, the rows where the basis has a single entry. The steady-state set is a
    cone, so no upper bounds are needed. Raises ValueError when no such flux exists, which means
    the network still has blocked reactions.
    """
    triples, (m, q) = entries(N)
    Nf = sparse.csr_matrix(([float(v) for _, _, v in triples],
                            ([i for i, _, _ in triples], [j for _, j, _ in triples])), shape=(m, q))
    rng = np.random.default_rng(seed)
    lower = 1 + rng.random(q)
    lp = linprog(rng.random(q), A_eq=Nf, b_eq=np.zeros(m),
                 bounds=[(lo, None) for lo in lower], method="highs")
    if lp.status != 0:
        raise ValueError("no strictly positive steady-state flux; remove blocked reactions first")

    rows = {}
    for i, j, v in entries(sd.sparse_nullspace(N))[0]:
        rows.setdefault(i, []).append((j, v))
    free = {}
    for i, row in rows.items():
        if len(row) == 1:
            j, v = row[0]
            free.setdefault(j, Fraction(float(lp.x[i])) / v)
    flux = [sum((v * free.get(j, 0) for j, v in rows.get(i, ())), Fraction(0)) for i in range(q)]
    if not all(x > 0 for x in flux):
        raise ValueError("snapped flux is not strictly positive")
    return flux


def _coupled_sets(A, Y, N, flux, detach_uptakes, equal=()):
    """Groups of complexes with proportional rows in ker(N diag(r*) S^T), with u_c = u_d imposed for
    every pair (c, d) in ``equal``."""
    n = A.shape[0]
    Acsc = A.tocsc()
    substrate = {}
    for j in range(Acsc.shape[1]):
        start, end = Acsc.indptr[j], Acsc.indptr[j + 1]
        for i, a in zip(Acsc.indices[start:end], Acsc.data[start:end]):
            if a < 0:
                substrate[j] = int(i)
    zero = {c for c in range(n) if Y.tocsc()[:, c].nnz == 0}

    by_reaction = {}
    triples, (m, _) = entries(N)
    for i, j, v in triples:
        by_reaction.setdefault(j, []).append((i, v))
    W, extra = {}, {}
    for j, s in substrate.items():
        col = extra.setdefault(j, n + len(extra)) if detach_uptakes and s in zero else s
        for i, v in by_reaction.get(j, ()):
            W[(i, col)] = W.get((i, col), Fraction(0)) + flux[j] * v
    ncol = n + len(extra)
    for k, (c, d) in enumerate(equal):
        W[(m + k, c)], W[(m + k, d)] = Fraction(1), Fraction(-1)
    Wx = sd.RationalMatrix.from_fractions(((i, j, x) for (i, j), x in W.items() if x != 0),
                                          (m + len(equal), ncol))
    groups = concordance_from_matrices(sparse.identity(ncol, format="csr", dtype=np.int64), Wx).modules
    drop = zero if detach_uptakes else set()
    return [set(c for c in g if c < n and c not in drop) for g in groups], set(substrate.values()) - drop


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
    included. Complexes that are no reaction's substrate appear in no rate law and stay singletons.
    """
    n = A.shape[0]
    flux = positive_flux(N, seed)
    equal, rounds = [], 0
    while True:
        coupled, substrates = _coupled_sets(A, Y, N, flux, detach_uptakes, equal)
        coupled = [s for s in coupled if len(s) > 1 or (len(s) == 1 and next(iter(s)) in substrates)]
        sets, known, r = _closure(Y, coupled)
        rounds += r
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
