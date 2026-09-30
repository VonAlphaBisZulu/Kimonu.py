"""Iterated psi route with a Descartes step (benchmark variant of kimonu.kinetic_modules).

Besides the psi-kernel readout and the log closure, every metabolite balance sum_j N_ij k_j psi_s(j) = 0 is
examined against the directions F already known to be fixed (composition differences within modules,
plus vectors proven fixed earlier). In u = psi / psi* (reference steady state from the interior flux r*),
if every substrate complex of the balance varies along one common direction t modulo F, then
u_c = exp(lambda_c t) exactly, and the balance reads sum_l a_l exp(l t) = 0 with a_l the sum of
N_ij r*_j over reactions whose substrate has lambda_c = l. If the nonzero a_l change sign exactly once
along increasing l, the sum has exactly one real root (Descartes' rule of signs, valid for real
exponents), which is t = 0: t -- and y_c for any term with lambda_c != 0 -- is fixed at every positive
steady state. The vector is added to F and the closure repeats; merges feed back into
the psi kernel as in the default route. Every step is a theorem, so the route stays globally sound.

    kinetic_modules_descartes(A, Y, N, detach_uptakes=False, seed=0) -> (modules, acr, acrr_groups, firings)
"""
from fractions import Fraction
import numpy as np
from scipy import sparse
from kimonu.kinetic import positive_flux, _psi_coupled
from kimonu.linalg import entries
from kimonu.acr import composition_differences
from kimonu.linalg import left_null_rows


def _rows(Y, sets, extra):
    YD = composition_differences(Y, sets)
    cols = [YD] + [sparse.csc_matrix(v.reshape(-1, 1)) for v in extra]
    YD = sparse.hstack(cols).tocsc() if extra else YD
    return left_null_rows(YD) if YD.shape[1] else None


def _signature(Ycsc, rows, c):
    acc = {}
    for met, coeff in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
        for col, v in (rows[int(met)] if rows is not None else ((int(met), Fraction(1)),)):
            acc[col] = acc.get(col, Fraction(0)) + Fraction(str(coeff)) * v
    return {k: v for k, v in acc.items() if v != 0}


def _closure(Y, sets, extra):
    """merge sets with equal signatures modulo F = span(YΔ, extra) until stable"""
    Ycsc = Y.tocsc()
    while True:
        rows = _rows(Y, sets, extra)
        groups = {}
        for s in sets:
            key = tuple(sorted(_signature(Ycsc, rows, min(s)).items()))
            groups.setdefault(key, set()).update(s)
        merged = list(groups.values())
        if len(merged) == len(sets):
            return merged, rows
        sets = merged


def _descartes(A, Y, N, rows, extra, flux, detach_uptakes):
    """one pass over the metabolite balances; returns newly fixed vectors"""
    Ycsc, Acsc = Y.tocsc(), A.tocsc()
    substrate = {}
    for j in range(Acsc.shape[1]):
        for i, a in zip(Acsc.indices[Acsc.indptr[j]:Acsc.indptr[j + 1]], Acsc.data[Acsc.indptr[j]:Acsc.indptr[j + 1]]):
            if a < 0: substrate[j] = int(i)
    triples, (m, _) = entries(N); by_row = {}
    for i, j, v in triples: by_row.setdefault(i, []).append((j, v))
    zero = {c for c in range(Ycsc.shape[1]) if Ycsc.indptr[c] == Ycsc.indptr[c + 1]}
    found = []
    for i, terms in by_row.items():
        # detached: every uptake has its own free source, so an uptake term is not a constant
        if detach_uptakes and any(substrate[j] in zero for j, _ in terms): continue
        coef = {}                                      # substrate complex -> sum of N_ij r*_j
        for j, v in terms:
            c = substrate[j]; coef[c] = coef.get(c, Fraction(0)) + v * flux[j]
        sig = {c: _signature(Ycsc, rows, c) for c in coef}
        ref = next((s for s in sig.values() if s), None)
        if ref is None: continue
        k0, v0 = next(iter(ref.items()))
        lam = {}
        for c, s in sig.items():
            if not s: lam[c] = Fraction(0); continue
            if set(s) != set(ref) or any(s[k] * v0 != ref[k] * s[k0] for k in s): lam = None; break
            lam[c] = s[k0] / v0
        if lam is None or len(set(lam.values())) < 2: continue
        by_l = {}
        for c, l in lam.items(): by_l[l] = by_l.get(l, Fraction(0)) + coef[c]
        signs = [1 if by_l[l] > 0 else -1 for l in sorted(by_l) if by_l[l] != 0]
        if len({l for l in by_l if by_l[l] != 0 and l != 0}) == 0: continue
        if sum(1 for a, b in zip(signs, signs[1:]) if a != b) == 1:
            c = next(c for c, l in lam.items() if l != 0)
            v = np.zeros(Y.shape[0])
            for met, coeff in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
                v[int(met)] = coeff
            found.append((i, v))
    return found


def kinetic_modules_descartes(A, Y, N, detach_uptakes=False, seed=0):
    n = A.shape[0]
    flux = positive_flux(N, seed)
    equal, extra, firings = [], [], []
    while True:
        coupled, subs = _psi_coupled(A, Y, N, flux, detach_uptakes, equal)
        coupled = [s for s in coupled if len(s) > 1 or (len(s) == 1 and next(iter(s)) in subs)]
        sets, rows = _closure(Y, coupled, extra)
        new = _descartes(A, Y, N, rows, extra, flux, detach_uptakes)
        if new:
            i, v = new[0]; extra.append(v); firings.append(i)
            continue
        if len(sets) == len(coupled):
            break
        equal = [(min(s), c) for s in sets for c in s if c != min(s)]
    rows = _rows(Y, sets, extra)
    m = Y.shape[0]
    acr = [i for i in range(m) if rows is not None and not rows[i]]
    groups = {}
    for i in range(m):
        if rows is not None and rows[i]: groups.setdefault(rows[i], []).append(i)
    assigned = set().union(*sets) if sets else set()
    modules = [s for s in sets if s] + [{i} for i in range(n) if i not in assigned]
    return modules, acr, [g for g in groups.values() if len(g) > 1], firings
