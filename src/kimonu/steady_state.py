"""The reference steady state: an exact interior flux r* and the mass balance W in its units (Eq. 3).

Under mass action, reaction j runs at r_j = k_j ψ_σ(j), with ψ_c the concentration of complex c (the product of
its metabolite concentrations, each raised to its coefficient, so ln ψ = Yᵀ ln x) and σ(j) the substrate complex
of j. The rate constants k are fixed but unknown. A reference steady state takes their place. Let r* be any
strictly positive flux in ker(N). The rate constants k = r* have the steady state x* = 1 with exactly these fluxes
and ψ* = 1, and every complex concentration is measured in units of this reference, u_c = ψ_c / ψ*_c. Every rate
then changes by the factor of its substrate complex, r_j = r*_j u_σ(j), and the mass balance becomes

    dx/dt = N diag(k) Sᵀ ψ = N diag(r*) Sᵀ u = W u = 0,                                         (Eq. 3)

with S (complexes x reactions) the substrate part of the incidence matrix, S_cj = 1 when c = σ(j). Column c of W is
the sum of the columns of N for the reactions that consume complex c, each weighted by its reference flux. W u = 0
is therefore the balance N r = 0, except that reactions that consume the same complex share one substrate
concentration and scale by the same factor u_c.

Both routes start from W (:func:`steady_state_map`). The iterated route reads coupled complexes from ker(W), and
the tangent route linearizes W u around u = 1. Uptakes are either attached, so they consume the zero complex with
ψ_∅ = 1 and run at the constant rate k_j, or detached (the default), so each draws from a source of its own whose
concentration is free and adds a column of its own to W (:func:`rate_columns`).

:func:`interior_flux` draws r* at random, which gives the result for generic rate constants, that is, for all rate
constants that admit a positive steady state except finely tuned ones. With the default uniform draw a
Schwartz–Zippel bound limits the probability that r* is finely tuned (Schwartz, 1980; Zippel, 1979), and the
keyword ``bits`` sets that bound.
"""
from fractions import Fraction
from math import ceil, lcm

import numpy as np
from scipy.optimize import linprog

from kimonu.io import substrate_complexes, zero_complexes
from kimonu.linalg import entries, null_space_rows, to_float

_LP_SCALE = 10**6           # the LP point's free coordinates are rounded to multiples of 1/_LP_SCALE


def interior_flux(N, seed=0, draw="uniform", bits=64):
    """An exact, strictly positive steady-state flux r* in ker(N), as a list of integers (Fractions), drawn at random.

    Both draws start from an LP point. HiGHS (bundled with scipy) solves N r = 0 with random lower bounds
    r >= 1 + U(0, 1) and a random objective. The steady states form a cone, so no upper bounds are needed. Forcing
    every rate to at least 1 spans rates over about eight orders of magnitude on genome-scale models, and HiGHS
    solved it on every model and seed tried, where SCIP, GLPK and CPLEX through straindesign often did not. The
    basis K of ker(N) (the RREF basis of :func:`kimonu.linalg.null_space_rows`) has 1 at each free coordinate and 0
    at the other free coordinates, so the LP point's coefficients in K are its entries at the free coordinates.
    They are rounded to multiples of 10^-6 (the scale is raised tenfold while a rate rounds to zero or below), the
    flux is formed exactly as K times these coefficients, and its common denominator is cleared. N r* = 0 therefore
    holds by construction, and r* > 0 is checked in exact arithmetic.

    Args:
        N: metabolites x reactions stoichiometric matrix, exact (straindesign RationalMatrix) or integer.
        seed: seed of the LP's random bounds and objective and of the uniform draw.
        draw: "uniform" (default) or "lp".
            "uniform". r* = K (T c + s), with c the integer free coordinates of the LP point, s drawn uniformly
            from the integers [0, 2^bits) in each free coordinate, and
            T = ceil(2 · 2^bits · max_i Σ_j |K_ij| / min r_c), with r_c = K c the LP point, large enough that r*
            stays strictly positive whatever s is. A shift by the fixed T c keeps the degree of a polynomial, so a
            nonzero polynomial of degree D in the free coordinates vanishes at r* with probability at most
            D / 2^bits (Schwartz–Zippel), and the bound applies to r* itself.
            "lp". r* = K c, the LP point. Its integers are smaller, so the tangent route runs several times faster,
            but an LP optimum is not a random draw, and no probability bound applies to it.
            On all 33 genome-scale networks both draws gave the same kinetic modules, ACR and ACRR for both routes.
        bits: size of the uniform draw, s in [0, 2^bits), used only for draw="uniform" (default 64). The bound
            D / 2^bits halves with every extra bit. Fewer bits make the tangent route faster, because the entries
            of r* are smaller.

    Raises ValueError when N has no strictly positive steady state, which means the network still has blocked
    reactions.
    """
    if draw not in ("uniform", "lp"):
        raise ValueError(f"draw must be 'uniform' or 'lp', not {draw!r}")
    if not (isinstance(bits, int) and bits >= 1):
        raise ValueError(f"bits must be a positive integer, not {bits!r}")
    Nf = to_float(N); m, q = Nf.shape
    rng = np.random.default_rng(seed)
    lower = 1 + rng.random(q)
    lp = linprog(rng.random(q), A_eq=Nf, b_eq=np.zeros(m), bounds=[(lo, None) for lo in lower], method="highs")
    if lp.status != 0:
        raise ValueError("no strictly positive steady-state flux; remove blocked reactions first")
    rows = null_space_rows(N)                    # K, one {basis vector: value} dict per reaction

    def free_coordinates(values):
        """{basis vector j: b_j} read off a flux in ker(N) at the rows of K with a single entry"""
        free = {}
        for i, row in enumerate(rows):
            if len(row) == 1:
                (j, v), = row.items()
                free.setdefault(j, values[i] / v)
        return free

    def combine(b):
        """K b, exactly"""
        return [sum((v * b.get(j, 0) for j, v in row.items()), Fraction(0)) for row in rows]

    # 1) the LP point, its free coordinates rounded at scale 10^6, as an integer flux
    scale = _LP_SCALE
    while True:
        flux = combine(free_coordinates([Fraction(round(x * scale)) for x in lp.x]))
        if all(x > 0 for x in flux):
            break
        scale *= 10
        if scale >= 10**12:
            raise ValueError("rounding did not give a strictly positive flux")
    denominator = lcm(*(x.denominator for x in flux))
    flux = [x * denominator for x in flux]
    if draw == "lp":
        return flux

    # 2) the uniform draw around it, b = T c + s with s in [0, 2^bits)
    c = free_coordinates(flux)
    draws = np.random.default_rng(_LP_SCALE + seed)
    s = {j: Fraction(int.from_bytes(draws.bytes((bits + 7) // 8), "big") % 2**bits) for j in c}
    reach = max(sum(abs(v) for v in row.values()) for row in rows)      # |K s|_i <= 2^bits · reach
    T = ceil(2 * 2**bits * reach / min(flux))
    flux = combine({j: T * c[j] + s[j] for j in c})
    if not all(x > 0 for x in flux):
        raise ValueError("the uniform draw did not give a strictly positive flux")
    denominator = lcm(*(x.denominator for x in flux))
    return [x * denominator for x in flux]


def rate_columns(A, Y, detach_uptakes=True):
    """The column of u that carries each reaction's rate, {reaction: column}, and the number of columns.

    Reaction j runs at r*_j u_c for its column c. That is its substrate complex σ(j), or, for a detached uptake, a
    source of its own after the n complexes, numbered in the order of the reactions.
    """
    n = A.shape[0]
    zero = zero_complexes(Y) if detach_uptakes else set()
    column, n_columns = {}, n
    for j, c in substrate_complexes(A).items():
        if c in zero:
            column[j], n_columns = n_columns, n_columns + 1
        else:
            column[j] = c
    return column, n_columns


def steady_state_map(A, Y, N, flux, detach_uptakes=True):
    """W = N diag(r*) Sᵀ of Eq. 3 as {(metabolite, column): value}, and its shape, so that dx/dt = W u.

    Column c collects, for every reaction that consumes complex c, its column of N times its rate at r*. A detached
    uptake gets a column of its own after the n complexes (:func:`rate_columns`), a free source rather than the
    constant ψ_∅ = 1. Both routes, :func:`kimonu.relation_values` (with k in place of r*) and
    :func:`kimonu.kinetic_flux_coupling` build on this map.
    """
    column, n_columns = rate_columns(A, Y, detach_uptakes)
    by_reaction = {}
    triples, (m, _) = entries(N)
    for i, j, v in triples:
        by_reaction.setdefault(j, []).append((i, v))
    W = {}
    for j, col in column.items():
        for i, v in by_reaction.get(j, ()):
            W[(i, col)] = W.get((i, col), Fraction(0)) + flux[j] * v
    return W, (m, n_columns)
