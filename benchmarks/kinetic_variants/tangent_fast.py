"""Tangent-space readout with an integer interior flux.

Same readout as tangent_one_shot.one_shot_exact, but r* is integral: the LP point is scaled, its free
coordinates (the rows where the kernel basis has a single entry) are rounded to integers, their common
denominator against the basis is cleared, and r* = K b is then an integer vector in ker N. It is checked
to be strictly positive and N r* = 0 exactly. The previous snapping kept the LP's binary fractions,
whose entries are already small (146 bits on iNJ661) but not integral; a first attempt at small decimal
denominators made them larger (1542 bits), since the common denominator of many small denominators
grows. Integers (64 bits on iNJ661) cut the tangent computation there from 99 s to 13 s.

Building the tangent space from the sparse block system [W 0; I -Y^T] instead of J = W Y^T was also
tried and was 1.3-5x slower (iIS312_Trypomastigote, iNJ661), so J is formed as before.

    one_shot_fast(A, Y, N, detach, seed=0, want_groups=False)
"""
from fractions import Fraction
from math import lcm
import numpy as np, scipy.sparse as sp, straindesign as sd
from scipy.optimize import linprog
from kimonu.kinetic import _exact_entries, _kernel_entries
import kinetic_variants.tangent_one_shot as _t


def integer_positive_flux(N, seed=0, scale=100):
    entries, (m, q) = _exact_entries(N)
    Nf = sp.csr_matrix(([float(v) for *_, v in entries], ([i for i, *_ in entries], [j for _, j, _ in entries])), shape=(m, q))
    rng = np.random.default_rng(seed)
    lp = linprog(rng.random(q), A_eq=Nf, b_eq=np.zeros(m), bounds=[(lo, None) for lo in 1 + rng.random(q)], method="highs")
    if lp.status != 0:
        raise ValueError("no strictly positive steady-state flux; remove blocked reactions first")
    rows = {}
    for i, j, v in _kernel_entries(sd.sparse_nullspace(N)):
        rows.setdefault(i, []).append((j, v))
    while scale < 1e12:
        free = {}
        for i, row in rows.items():
            if len(row) == 1:
                j, v = row[0]; free.setdefault(j, Fraction(round(lp.x[i] * scale)) / v)
        L = 1
        for b in free.values(): L = lcm(L, b.denominator)
        flux = [sum((v * free.get(j, 0) * L for j, v in rows.get(i, ())), Fraction(0)) for i in range(q)]
        if all(x > 0 for x in flux):
            res = {}
            for i, j, v in entries: res[i] = res.get(i, Fraction(0)) + v * flux[j]
            assert not any(res.values()), "integer flux is not an exact steady state"
            return flux
        scale *= 10
    raise ValueError("rounding did not give a strictly positive flux")


def one_shot_fast(A, Y, N, detach, seed=0, want_groups=False):
    saved = _t.positive_flux
    _t.positive_flux = integer_positive_flux
    try:
        return _t.one_shot_exact(A, Y, N, detach, seed=seed, want_groups=want_groups)
    finally:
        _t.positive_flux = saved
