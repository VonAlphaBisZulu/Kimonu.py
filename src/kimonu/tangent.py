"""Kinetic modules and local ACR by the tangent route: one null space at one generic steady state.

The iterated route (kinetic.py) builds the fixed log-combinations up from relations it can prove at
every steady state. The tangent route reads them all off the linearization at a single steady state.
With rate constants fixed by an interior flux r* at x* = 1, as in kinetic.py, the steady-state equations
in log coordinates z = ln x are

    F(z) = N diag(r*) Sᵀ exp(Yᵀ z) = W u(z) = 0,

and since u = exp(Yᵀ z) has derivative Yᵀ at u = 1, the chain rule gives their Jacobian at z* = 0 as

    J = N diag(r*) Sᵀ Yᵀ = W Yᵀ.

T = ker(J) is the tangent space of the steady states at x*: the directions in which ln x can move
without leaving them. A log-combination vᵀ ln x is fixed exactly when v is orthogonal to T, so the
readout of acr.py with B = T gives the kinetic modules (equal rows of Yᵀ T), the metabolites with ACR
(zero rows of T) and the ACRR groups (equal nonzero rows). One null space, no iteration.

Why nothing is missed. A steady state z near z* is itself the point z* = 0 for the rate constants
fixed by its own flux r* ∘ exp(Sᵀ Yᵀ z), so the tangent spaces at nearby steady states are the T of
nearby interior fluxes. For generic r* the rank of J does not change between nearby interior fluxes, and
the constant-rank theorem then makes the steady states near x* a manifold with tangent space T: a
combination that is constant on them is orthogonal to T, so no locally fixed relation is missed.
Conversely, each candidate v (e_i, e_i - e_j or y_c - y_d: finitely many, the same for every r*) is
orthogonal to T exactly when a polynomial condition on r* holds; holding at a generic r*, it holds at
every interior flux nearby, so v is orthogonal to the tangent space at every steady state near x*, and
vᵀ ln x is constant there. Neither direction needs the steady state to be nondegenerate.

What is found is local ACR in the sense of Pascual-Escudero & Feliu (2022): constant on the steady
states near one of them. When the same rate constants admit several separate steady states, a
relation can hold near each of them with a different value. The cooperative kinase 0 <-> X, 2X -> C,
C -> 0, C + 2P -> 3P, P -> 0 has the two steady states X = 10, (C, P) = (1, 4) and (4, 1) for one choice
of rate constants; T = {0} at both, so this route calls X, C and P fixed, while only X is. The
iterated route claims neither C nor P, and misses X: for generic rate constants,
iterated ⊆ true ⊆ tangent.

The zero complex has ψ = 1. Detached (default), each uptake reaction gets a source coordinate of its own,
a free variable in T; only the metabolite rows of T are read, and the zero complex is a module of its own.
Attached, its uptake reactions run at constant rate and add nothing to J.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np
from scipy import sparse

from kimonu.acr import readout
from kimonu.io import zero_complexes
from kimonu.kinetic import steady_state_map
from kimonu.linalg import exact, null_space_rows, positive_flux

TangentResult = namedtuple(
    "TangentResult",
    ["modules", "module_of", "n_modules", "giant_module", "acr", "acrr", "dimension", "flux"],
)


def tangent_route(A, Y, N, detach_uptakes=True, seed=0, metabolite_ids=None):
    """Kinetic modules, local ACR and local ACRR from the tangent space of the steady states.

    Args:
        A: complexes x reactions incidence matrix, one -1 and one +1 per column.
        Y: metabolites x complexes composition matrix.
        N: metabolites x reactions stoichiometric matrix, exact (straindesign RationalMatrix) or
            integer; N = Y A. The network must have no blocked reactions.
        detach_uptakes: give every reaction out of the zero complex its own source, as in
            :func:`kimonu.kinetic_modules`.
        seed: seed of the random interior flux; results hold for generic rate constants.
        metabolite_ids: optional names for ``acr`` and ``acrr``; indices are returned when omitted.

    Returns a :class:`TangentResult`. ``modules`` partition the complexes like those of
    :func:`kimonu.kinetic_modules`; ``acr`` are the metabolites with local ACR and ``acrr`` the
    sorted pairs with local ACRR, both in the sense of Pascual-Escudero & Feliu (2022), so they can
    claim relations that differ between separate steady states (see the module docstring);
    ``dimension`` is dim T. ``flux`` is the interior flux r* the route ran at, so the relations are local
    around the steady state x = 1 for the rate constants k = r* (see :func:`kimonu.rate_constants`).
    """
    n, m = A.shape[0], Y.shape[0]
    skip = zero_complexes(Y) if detach_uptakes else set()
    # 1) Generic rate constants from an interior flux r*, made integral at scale 10^6: the exact null
    #    space of J is several times faster on integers, and at scale 100 the rounding was coarse
    #    enough to make r* non-generic (a spurious ACRR pair on iHN637); from 10^4 on the results agree.
    flux = positive_flux(N, seed, scale=10**6)
    # 2) J = W Yᵀ: the column of complex c spreads over c's metabolites; a detached source keeps its own
    #    column, after the m metabolites.
    W, (_, n_columns) = steady_state_map(A, Y, N, flux, detach_uptakes)
    Ycsc = sparse.csc_matrix(Y)
    J = {}
    for (i, c), w in W.items():
        if c >= n:
            J[(i, m + c - n)] = w
            continue
        for k, y in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
            J[(i, int(k))] = J.get((i, int(k)), Fraction(0)) + w * Fraction(str(y))
    # 3) T = ker(J), one row per coordinate; the metabolite rows are B
    rows = null_space_rows(exact(((i, k, v) for (i, k), v in J.items()), (m, m + n_columns - n)))
    found = readout(Y, [tuple(sorted(r.items())) for r in rows[:m]], skip)
    # 4) A detached zero complex is coupled to nothing and forms a module of its own.
    modules = found.modules + [{c} for c in sorted(skip)]
    module_of = np.zeros(n, dtype=np.int64)
    for idx, group in enumerate(modules, start=1):
        for i in group:
            module_of[i] = idx

    names = metabolite_ids if metabolite_ids is not None else list(range(m))
    acrr = [(names[a], names[b]) for g in found.acrr for k, a in enumerate(g) for b in g[k + 1:]]
    return TangentResult(
        modules=modules,
        module_of=module_of,
        n_modules=sum(1 for g in modules if len(g) > 1),
        giant_module=max((len(g) for g in modules), default=0),
        acr=[names[i] for i in found.acr],
        acrr=sorted(acrr),
        dimension=len({j for r in rows for j in r}),
        flux=flux,
    )
