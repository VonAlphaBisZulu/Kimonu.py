"""Toric (binomial) check on the psi kernel, before and after the iterated route's feedback.

The steady-state equations are W u = 0 with W = N diag(r*) S^T (plus the rows u_c = u_d the iterated
route adds). If ker W has a basis with pairwise disjoint supports -- equivalently, its proportional-row
classes are linearly independent, so dim ker W equals the number of classes -- then W's row space is
spanned by binomials, the positive steady states are log-linear (toric, Perez Millan et al. 2012) and
connected, and every relation the tangent space finds is global.

    python kinetic_variants/toric_check.py <model|envz|langary|cubic> [...]
"""
import json, os, sys
from fractions import Fraction
import numpy as np, scipy.sparse as sp, straindesign as sd
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src")); sys.path.insert(0, os.path.join(os.path.dirname(B), "tests"))
import kimonu
from kimonu.kinetic import positive_flux, _exact_entries, _closure
from kimonu.concordance import concordance_from_matrices

def psi_kernel(A, Y, N, flux, equal=()):
    n = A.shape[0]; Acsc = A.tocsc(); substrate = {}
    for j in range(Acsc.shape[1]):
        for i, a in zip(Acsc.indices[Acsc.indptr[j]:Acsc.indptr[j + 1]], Acsc.data[Acsc.indptr[j]:Acsc.indptr[j + 1]]):
            if a < 0: substrate[j] = int(i)
    entries, (m, _) = _exact_entries(N); by = {}
    for i, j, v in entries: by.setdefault(j, []).append((i, v))
    W = {}
    for j, s in substrate.items():
        for i, v in by.get(j, ()): W[(i, s)] = W.get((i, s), Fraction(0)) + flux[j] * v
    for k, (c, d) in enumerate(equal): W[(m + k, c)], W[(m + k, d)] = Fraction(1), Fraction(-1)
    Wx = sd.RationalMatrix.from_fractions(((i, j, x) for (i, j), x in W.items() if x != 0), (m + len(equal), n))
    res = concordance_from_matrices(sp.identity(n, format="csr", dtype=np.int64), Wx)
    return res, set(substrate.values())

def check(A, Y, N):
    flux = positive_flux(N)
    out = []
    equal = []
    while True:
        res, subs = psi_kernel(A, Y, N, flux, equal)
        out.append((res.nullity, len(res.modules)))
        coupled = [set(g) for g in res.modules]
        coupled = [s for s in coupled if len(s) > 1 or next(iter(s)) in subs]
        sets, known, _ = _closure(Y, coupled)
        if len(sets) == len(coupled): break
        equal = [(min(s), c) for s in sets for c in s if c != min(s)]
    return out

if __name__ == "__main__":
    from run_kimonu_crn import triplets_to_sparse, to_exact
    from test_kinetic import _crn, ENVZ, ENVZ_SPECIES
    TOY = {"cubic": ([("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")], ["A"]), "envz": (ENVZ, ENVZ_SPECIES),
           "langary": ([('A','B'),('B','A'),('B','C'),('C','B'),('C','D'),('D+E','F'),('F','D+E'),('F','B+G'),('C+G','H'),('H','C+G'),('H','C+E'),('A+G','I'),('I','A+G'),('I','A+E')], list("ABCDEFGHI"))}
    for name in sys.argv[1:]:
        if name in TOY:
            A, Y, N, _ = _crn(*TOY[name])
        else:
            crn = json.load(open(f"{B}/unblocked/{name}_fixed.crn.json")); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
            A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc); N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
        rounds = check(A, Y, N)
        first, last = rounds[0], rounds[-1]
        tor = lambda r: "toric" if r[0] == r[1] else f"not toric ({r[1] - r[0]} relations among classes)"
        print(f"{name:24s} initial W: dim {first[0]:5d}, classes {first[1]:5d} -> {tor(first):40s} | after feedback: dim {last[0]:5d}, classes {last[1]:5d} -> {tor(last)}", flush=True)
