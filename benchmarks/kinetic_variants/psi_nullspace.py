"""Kinetic modules from a nullspace in complex-concentration space, against the upstream algorithm.

Mass action: v_r = k_r ψ_s(r), with ψ_c the monomial of complex c. Fix any positive steady state
(v*, ψ*) and write ψ = ψ*·u. The rate constants k_r = v*_r / ψ*_s(r) then admit ψ* as a steady state,
and every steady state for those k satisfies Y·A·diag(v*)·Sᵀ·u = N·diag(v*)·Sᵀ·u = 0. So u lies in
ker(W), W = N·diag(v*)·Sᵀ (metabolites x complexes), and u = 1 is in it. Two complexes whose rows of a
basis of ker(W) are proportional have a fixed ψ ratio at every positive steady state of those k: a
sound coupling test with no graph traversal. v* is drawn at random from the interior of the flux
cone (every reaction positive), exact over the rationals; two independent draws are compared.

    python kinetic_variants/psi_nullspace.py <model> [...]   -> kinetic_variants/psi_<model>.json
"""
import json, sys, os, time, random
from fractions import Fraction
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B); sys.path.insert(0, os.path.join(os.path.dirname(B), "src"))
import kimonu, straindesign as sd
from compat import upstream
from run_kimonu_crn import triplets_to_sparse, to_exact

def blocks(groups): return {frozenset(g) for g in groups if len(g) > 1}
def relation(P, Q):
    if P == Q: return "identical"
    ref = lambda X, Z: all(any(b <= c for c in Z) for b in X)
    return "finer" if ref(P, Q) else "coarser" if ref(Q, P) else "neither"

def exact_positive_flux(crn, Sx, Kx, seed):
    """Exact v in ker(N), every entry > 0: an LP point with v >= 1, snapped onto the exact kernel
    through its free coordinates (the basis carries an identity block on them)."""
    nm, nr = crn["n_metabolites"], crn["n_reactions"]
    rng = np.random.default_rng(seed)
    Nf = triplets_to_sparse(crn["S"], nm, nr)
    lo = 1 + rng.random(nr)          # the steady-state set is a cone, so upper bounds only get in the way
    r = linprog(rng.random(nr), A_eq=Nf, b_eq=np.zeros(nm), bounds=[(l, None) for l in lo], method="highs")
    assert r.status == 0, r.message
    cols = {}
    for i, j, v in zip(Kx.rows, Kx.cols, Kx.data): cols.setdefault(int(i), []).append((int(j), int(v)))
    denom = int(getattr(Kx, "denom", 1) or 1)
    b = {}
    for i, entries in cols.items():
        if len(entries) == 1:
            j, c = entries[0]
            b.setdefault(j, Fraction(float(r.x[i])) * denom / c)
    v = [sum((Fraction(c, denom) * b[j] for j, c in cols.get(i, [])), Fraction(0)) for i in range(nr)]
    assert all(x > 0 for x in v), "snapped flux not strictly positive"
    return v

DETACH = '--detach' in sys.argv
for M in [a for a in sys.argv[1:] if not a.startswith('--')]:
  try:
      crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
      nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
      A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
      Sx = to_exact(triplets_to_sparse(crn["S"], nm, nr))
      t0 = time.time()
      res = kimonu.concordance_from_matrices(A, Sx)
      k_over = upstream.kinetic_modules(A, res.modules, res.balanced, Y=None)
      k_def = upstream.kinetic_modules(A, res.modules, res.balanced, Y=Y, nullity=res.nullity)
      t_kim = time.time() - t0
      Kx = sd.sparse_nullspace(Sx)
      if not isinstance(Kx, sd.ExactCOO):
          Kc = Kx.tocoo(); Kx = type("K", (), dict(rows=Kc.row, cols=Kc.col, data=Kc.data, denom=1))()
      sub = {}
      for i, j, v in zip(*sp.find(A)):
          if v < 0: sub[int(j)] = int(i)
      Nent = {}
      for i, j, v in crn["S"]: Nent.setdefault(int(j) - 1, []).append((int(i) - 1, Fraction(str(v))))   # 1-based triplets
      out = dict(model=M, n_complexes=nc); P = {}
      for draw in (1, 2):
          v = exact_positive_flux(crn, Sx, Kx, draw)
          # detached: every reaction out of the zero complex gets its own source column, so no
          # complex can be coupled through the constant monomial of the zero complex
          zero = [c for c in range(nc) if Y[:, c].nnz == 0]
          extra = {}
          W = {}
          for r_, s in sub.items():
              col = s
              if DETACH and s in zero:
                  col = extra.setdefault(r_, nc + len(extra))
              for i, a in Nent.get(r_, []):
                  W[(i, col)] = W.get((i, col), Fraction(0)) + v[r_] * a
          ncol = nc + len(extra)
          Wx = sd.RationalMatrix.from_fractions(((i, j, x) for (i, j), x in W.items() if x != 0), (nm, ncol))
          t0 = time.time()
          psi = kimonu.concordance_from_matrices(sp.identity(ncol, format="csr", dtype=np.int64), Wx)
          psi = psi._replace(modules=[[c for c in m if c < nc and not (DETACH and c in zero)] for m in psi.modules])
          out[f"seconds_psi_{draw}"] = time.time() - t0
          P[draw] = blocks(psi.modules)
          out[f"draw{draw}"] = dict(n_modules=len(P[draw]), giant=max((len(x) for x in P[draw]), default=0), n_zero=len(psi.balanced),
                                    vs_overlap=relation(P[draw], blocks(k_over.modules)), vs_default=relation(P[draw], blocks(k_def.modules)))
      out["draws_agree"] = P[1] == P[2]
      O = blocks(k_over.modules)
      out.update(overlap=dict(n_modules=len(O), giant=k_over.giant_module), default=dict(n_modules=len(blocks(k_def.modules)), giant=k_def.giant_module),
                 seconds_kimonu=t_kim, only_psi=sorted((len(x) for x in P[1] - O), reverse=True)[:15], only_overlap=sorted((len(x) for x in O - P[1]), reverse=True)[:15])
      json.dump(out, open(f"{B}/kinetic_variants/psi_{M}{'_detached' if DETACH else ''}.json", "w"), indent=1)
      d = out["draw1"]
      print(f"{M:26s}{' det' if DETACH else '    '} n={nc:5d} psi mods {d['n_modules']:5d} giant {d['giant']:4d} zero {d['n_zero']:3d} | overlap mods {len(O):5d} giant {k_over.giant_module:4d} -> {d['vs_overlap']:9s} | default giant {k_def.giant_module:4d} -> {d['vs_default']:9s} | draws agree {out['draws_agree']} | psi {out['seconds_psi_1']:.1f}s kimonu {t_kim:.1f}s", flush=True)
  except Exception as e:
    print(f'{M:26s} FAILED {type(e).__name__}: {e}', flush=True)
