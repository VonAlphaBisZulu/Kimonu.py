"""Check: if the mass-action deficiency dim(ker Y ∩ im L) is one, ker(W) alone couples all nonterminal
complexes (Langary et al. Theorem S4-6 / Proposition S4-8 would add nothing).
Random small networks with a positive steady-state flux; exact arithmetic. See proof_upstream_in_psi.md.

    python deficiency_one_check.py
"""
import random
import numpy as np, sympy as sp
from scipy.sparse.csgraph import connected_components
from scipy import sparse

def positive_flux(N, rng):
    K = sp.Matrix.hstack(*N.nullspace()) if N.nullspace() else None
    if K is None: return None
    # random positive combination of extreme-ish vectors: sample until all positive
    for _ in range(200):
        w = sp.Matrix([sp.Rational(rng.randint(-5, 9), rng.randint(1, 5)) for _ in range(K.cols)])
        r = K * w
        if all(x > 0 for x in r): return r
    return None

rng = random.Random(1)
hits = ok = 0
big = []
tried = 0
while hits < 300 and tried < 200000:
    tried += 1
    ns, nc, nr = rng.randint(2, 4), rng.randint(3, 7), rng.randint(3, 8)
    comps = set()
    while len(comps) < nc:
        comps.add(tuple(rng.randint(0, 2) for _ in range(ns)))
    comps = list(comps)
    reacs = set()
    for _ in range(nr):
        a, b = rng.sample(range(nc), 2); reacs.add((a, b))
    reacs = list(reacs)
    used = sorted({c for r in reacs for c in r})
    if len(used) < nc: continue
    A = sp.zeros(nc, len(reacs))
    for j, (a, b) in enumerate(reacs): A[a, j] = -1; A[b, j] = 1
    Y = sp.Matrix([[comps[c][i] for c in range(nc)] for i in range(ns)])
    N = Y * A
    r = positive_flux(N, rng)
    if r is None: continue
    S = sp.zeros(nc, len(reacs))
    for j, (a, b) in enumerate(reacs): S[a, j] = 1
    L = A * sp.diag(*r) * S.T
    W = Y * L
    imL = L.columnspace()
    # dim(ker Y ∩ im L) = dim im L - dim Y(im L)
    ML = sp.Matrix.hstack(*imL)
    dk = ML.rank() - (Y * ML).rank()
    if dk != 1: continue
    hits += 1
    # terminal strong classes
    adj = sparse.csr_matrix(np.array([[1 if (a, b) in reacs else 0 for b in range(nc)] for a in range(nc)]))
    k, lab = connected_components(adj, directed=True, connection="strong")
    term = set()
    for s in range(k):
        mem = [c for c in range(nc) if lab[c] == s]
        if not any((a, b) in reacs and lab[b] != s for a in mem for b in range(nc)):
            term |= set(mem)
    nonterm = [c for c in range(nc) if c not in term]
    Kw = sp.Matrix.hstack(*W.nullspace())
    rows = [Kw.row(c) for c in nonterm]
    big.append(len(nonterm))
    if all(sp.Matrix.vstack(rows[0], x).rank() == 1 for x in rows[1:]):
        ok += 1
    else:
        print("COUNTEREXAMPLE", comps, reacs)
print(f"tried {tried}, deficiency-one networks {hits}, nonterminal complexes all coupled by ker(W): {ok}")
import collections
print("nonterminal complexes per network:", sorted(collections.Counter(big).items()))
