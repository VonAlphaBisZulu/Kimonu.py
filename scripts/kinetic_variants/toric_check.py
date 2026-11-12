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
import records                                         # paths of the records (KIMONU_RECORDS, KIMONU_OUT)
import kimonu
from kimonu import free_directions, interior_flux, row_readout
from kimonu.linalg import entries, matrix_product, nullspace

def psi_kernel(A, Y, N, flux, equal=()):
    n = A.shape[0]; Acsc = A.tocsc(); substrate = {}
    for j in range(Acsc.shape[1]):
        for i, a in zip(Acsc.indices[Acsc.indptr[j]:Acsc.indptr[j + 1]], Acsc.data[Acsc.indptr[j]:Acsc.indptr[j + 1]]):
            if a < 0: substrate[j] = int(i)
    triples, (m, _) = entries(N); by = {}
    for i, j, v in triples: by.setdefault(j, []).append((i, v))
    W = {}
    for j, s in substrate.items():
        for i, v in by.get(j, ()): W[(i, s)] = W.get((i, s), Fraction(0)) + flux[j] * v
    for k, (c, d) in enumerate(equal): W[(m + k, c)], W[(m + k, d)] = Fraction(1), Fraction(-1)
    Wx = sd.RationalMatrix.from_fractions(((i, j, x) for (i, j), x in W.items() if x != 0), (m + len(equal), n))
    K = nullspace(Wx)                             # groups of proportional nonzero rows of ker(W)
    return len({j for row in K for j, _ in row}), row_readout(K).groups, set(substrate.values())

def check(A, Y, N):
    net = kimonu.crn_from_matrices(A, Y)
    flux = interior_flux(net)
    out = []
    equal = []
    while True:
        nullity, classes, subs = psi_kernel(A, Y, N, flux, equal)
        out.append((nullity, len(classes)))
        coupled = [set(g) for g in classes]
        coupled = [s for s in coupled if len(s) > 1 or next(iter(s)) in subs]
        YtB = matrix_product(net.Yt, free_directions(net, coupled))      # Yᵀ B, read out only for the sets
        groups = {}
        for c in sorted(set().union(*coupled)):
            groups.setdefault(YtB[c], set()).add(c)
        sets = list(groups.values())
        if len(sets) == len(coupled): break
        equal = [(min(s), c) for s in sets for c in s if c != min(s)]
    return out

if __name__ == "__main__":
    from run_kimonu_crn import triplets_to_sparse, to_exact
    from networks import matrices, ENVZ, ENVZ_SPECIES
    TOY = {"cubic": ([("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")], ["A"]), "envz": (ENVZ, ENVZ_SPECIES),
           "langary": ([('A','B'),('B','A'),('B','C'),('C','B'),('C','D'),('D+E','F'),('F','D+E'),('F','B+G'),('C+G','H'),('H','C+G'),('H','C+E'),('A+G','I'),('I','A+G'),('I','A+E')], list("ABCDEFGHI"))}
    for name in sys.argv[1:]:
        if name in TOY:
            A, Y, N, _ = matrices(*TOY[name])
        else:
            crn = json.load(open(records.network(name))); nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
            A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc); N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
        rounds = check(A, Y, N)
        first, last = rounds[0], rounds[-1]
        tor = lambda r: "toric" if r[0] == r[1] else f"not toric ({r[1] - r[0]} relations among classes)"
        print(f"{name:24s} initial W: dim {first[0]:5d}, classes {first[1]:5d} -> {tor(first):40s} | after feedback: dim {last[0]:5d}, classes {last[1]:5d} -> {tor(last)}", flush=True)
