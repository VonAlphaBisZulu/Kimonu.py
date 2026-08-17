"""Build the chemical-reaction-network matrices from a cobra model.

A metabolic network factorizes as N = Y·A:
  N : metabolites x reactions   stoichiometric matrix
  Y : metabolites x complexes   composition (which species, and how many, form each complex)
  A : complexes  x reactions    incidence (-1 on a reaction's substrate complex, +1 on its product)

A *complex* is one side of a reaction — a multiset of species (e.g. A+B, 2B, C). Boundary/exchange
reactions have an empty complex on one side, which contributes no incidence entry (the "zero complex",
as in COCOA's `:empty`). Complexes are enumerated in first-seen order over the reactions, matching
COCOA's deterministic ordering.
"""
from fractions import Fraction
from math import gcd

import numpy as np
from scipy import sparse


def _exact(coeff):
    """Exact rational value of a stoichiometric coefficient.

    ``Fraction(str(x))`` reproduces the decimal the model declares (1.496 -> 1496/1000), whereas
    ``Fraction(x)`` would recover the binary approximation of it. Coefficients are frequently
    fractional — biomass reactions and half-oxygen terms such as 0.5 O2 in CYTBD — so this is the
    difference between analysing the model and analysing a rounded version of it.
    """
    return Fraction(str(coeff))


def build_crn(model):
    """Return a namespace with A, Y, N and the id lists.

    Returns an object with attributes:
      A            scipy.csr int  (complexes x reactions), incidence
      Y            scipy.csr float (metabolites x complexes), composition
      N            scipy.csr float (metabolites x reactions) = Y @ A
      N_integer    scipy.csr int  (metabolites x reactions), N with every row scaled to integers
      complexes    list of canonical complexes (tuple of (met_id, Fraction coeff)), in A's row order
      reaction_ids, metabolite_ids

    ``N_integer`` is what the exact nullspace is computed from. Scaling a row of N by a positive
    constant leaves ker(N) unchanged, so clearing denominators row by row is exact and keeps the
    kernel computation in integer arithmetic.
    """
    met_ids = [m.id for m in model.metabolites]
    met_idx = {m: i for i, m in enumerate(met_ids)}
    rxn_ids = [r.id for r in model.reactions]

    complex_index = {}   # canonical key -> complex row index
    complex_comp = []    # complex row index -> {met_id: coeff (>0)}

    def _complex(comp):
        """Return the complex index for a multiset {met_id: coeff>0}, or None if empty."""
        if not comp:
            return None
        key = tuple(sorted(comp.items()))
        idx = complex_index.get(key)
        if idx is None:
            idx = len(complex_comp)
            complex_index[key] = idx
            complex_comp.append(comp)
        return idx

    Ai, Aj, Av = [], [], []
    for j, r in enumerate(model.reactions):
        substrate = {m.id: -_exact(c) for m, c in r.metabolites.items() if c < 0}   # |coeff| on the substrate side
        product = {m.id: _exact(c) for m, c in r.metabolites.items() if c > 0}      # coeff on the product side
        si = _complex(substrate)
        pi = _complex(product)
        if si is not None:
            Ai.append(si); Aj.append(j); Av.append(-1)
        if pi is not None:
            Ai.append(pi); Aj.append(j); Av.append(+1)

    ncplx = len(complex_comp)
    A = sparse.csr_matrix((Av, (Ai, Aj)), shape=(ncplx, len(rxn_ids)), dtype=np.int64)

    Yi, Yj, Yv = [], [], []
    for c, comp in enumerate(complex_comp):
        for m, coeff in comp.items():
            Yi.append(met_idx[m]); Yj.append(c); Yv.append(coeff)
    Y = sparse.csr_matrix(([float(v) for v in Yv], (Yi, Yj)),
                          shape=(len(met_ids), ncplx), dtype=np.float64)

    N = (Y @ A).tocsr()  # == stoichiometric matrix, by construction

    # Exact integer N for the kernel: accumulate the rational entries, then clear denominators one
    # metabolite row at a time.
    exact = {}
    for c, comp in enumerate(complex_comp):
        for m, coeff in comp.items():
            i = met_idx[m]
            for k in range(A.indptr[c], A.indptr[c + 1]):
                j = A.indices[k]
                exact[(i, j)] = exact.get((i, j), Fraction(0)) + Fraction(int(A.data[k])) * coeff

    Ni, Nj, Nv = [], [], []
    row_terms = {}
    for (i, j), v in exact.items():
        if v:
            row_terms.setdefault(i, []).append((j, v))
    for i, terms in row_terms.items():
        scale = 1
        for _, v in terms:
            scale = scale * v.denominator // gcd(scale, v.denominator)
        for j, v in terms:
            Ni.append(i); Nj.append(j); Nv.append(int(v * scale))
    N_integer = sparse.csr_matrix((Nv, (Ni, Nj)), shape=(len(met_ids), len(rxn_ids)), dtype=np.int64)

    complexes = [tuple(sorted(comp.items())) for comp in complex_comp]

    class CRN:
        pass
    out = CRN()
    out.A, out.Y, out.N, out.N_integer = A, Y, N, N_integer
    out.complexes = complexes
    out.reaction_ids = rxn_ids
    out.metabolite_ids = met_ids
    return out
