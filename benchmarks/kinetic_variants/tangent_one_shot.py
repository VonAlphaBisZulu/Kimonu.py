"""Exact one-shot: T = ker(J), J = N diag(r*) [G | E], G rows = y_substrate(j)^T, E = detached sources.
Coupled complexes: equal rows of Y^T T_xi; ACR: zero rows of T_xi; ACRR: equal rows of T_xi."""
import numpy as np, scipy.sparse as sp, straindesign as sd
from fractions import Fraction
from kimonu import interior_flux
from kimonu.linalg import entries


def positive_flux(N, seed=0):
    """r* of the one-shot readout. tangent_fast.py swaps in its own integer flux."""
    return interior_flux(N, seed)

def one_shot_exact(A, Y, N, detach, seed=0, want_groups=False):
    n = A.shape[0]; Ac = sp.csc_matrix(A); Yc = sp.csc_matrix(Y)
    flux = positive_flux(N, seed)
    sub = {}
    for j in range(Ac.shape[1]):
        for i, a in zip(Ac.indices[Ac.indptr[j]:Ac.indptr[j+1]], Ac.data[Ac.indptr[j]:Ac.indptr[j+1]]):
            if a < 0: sub[j] = int(i)
    zero = {c for c in range(n) if Yc[:, c].nnz == 0}
    ycol = {c: [(int(i), Fraction(str(v))) for i, v in zip(Yc.indices[Yc.indptr[c]:Yc.indptr[c+1]], Yc.data[Yc.indptr[c]:Yc.indptr[c+1]])] for c in range(n)}
    triples, (m, q) = entries(N)
    byr = {}
    for i, j, v in triples: byr.setdefault(j, []).append((i, v))
    J = {}; extra = {}
    for j, s in sub.items():
        cols = [(k, w) for k, w in ycol[s]]
        if detach and s in zero: cols = [(m + extra.setdefault(j, len(extra)), Fraction(1))]
        for i, v in byr.get(j, ()):
            for k, w in cols:
                J[(i, k)] = J.get((i, k), Fraction(0)) + flux[j] * v * w
    ncol = m + len(extra)
    Jx = sd.RationalMatrix.from_fractions(((i, k, x) for (i, k), x in J.items() if x != 0), (m, ncol))
    K = sd.sparse_nullspace(Jx)
    rows = {}
    for i, j, v in entries(K)[0]:
        if i < m: rows.setdefault(i, {})[j] = v
    key = lambda d: tuple(sorted((c, v) for c, v in d.items() if v != 0))
    acr = [i for i in range(m) if not rows.get(i)]
    groups = {}
    for c in sorted(set(sub.values())):
        if detach and c in zero: continue
        acc = {}
        for i, w in ycol[c]:
            for j, v in rows.get(i, {}).items(): acc[j] = acc.get(j, 0) + w * v
        groups.setdefault(key(acc), set()).add(c)
    mets = {}
    for i in range(m):
        if rows.get(i): mets.setdefault(key(rows[i]), []).append(i)
    acrr_groups = [g for g in mets.values() if len(g) > 1]
    acrr = sum(len(g) * (len(g) - 1) // 2 for g in acrr_groups)
    if want_groups:
        return list(groups.values()), acr, acrr_groups, K.shape[1]
    return list(groups.values()), acr, acrr, K.shape[1]
