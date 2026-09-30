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

import numpy as np
from scipy import sparse

import straindesign as sd


def _exact(coeff):
    """Exact rational value of a stoichiometric coefficient.

    ``Fraction(str(x))`` reproduces the decimal the model declares (1.496 -> 1496/1000), whereas
    ``Fraction(x)`` would recover the binary approximation of it. Coefficients are frequently
    fractional — biomass reactions and half-oxygen terms such as 0.5 O2 in CYTBD — so this is the
    difference between analysing the model and analysing a rounded version of it.
    """
    return Fraction(str(coeff))


def exact_matrix(entries, shape):
    """A straindesign RationalMatrix from ``(row, col, value)`` triples, holding them exactly.

    Clearing denominators to int64 instead overflows on models whose coefficients are written as
    15-significant-digit decimals: iYS854's -7.73333333333333e-08 has a denominator of 10**22, and
    Chlorella_variabilis_iAJ526 mixes seven such terms in one row.
    """
    return sd.RationalMatrix.from_fractions(
        ((int(i), int(j), v if isinstance(v, Fraction) else Fraction(str(v))) for i, j, v in entries),
        shape)


def build_crn(model):
    """Return a namespace with A, Y, N and the id lists.

    Returns an object with attributes:
      A            scipy.csr int  (complexes x reactions), incidence
      Y            scipy.csr float (metabolites x complexes), composition
      N            scipy.csr float (metabolites x reactions) = Y @ A
      N_exact      straindesign.RationalMatrix (metabolites x reactions), N with exact entries
      complexes    list of canonical complexes (tuple of (met_id, Fraction coeff)), in A's row order
      reaction_ids, metabolite_ids

    ``N_exact`` is what the exact nullspace is computed from.
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

    # Exact N for the kernel, accumulated from the rational compositions.
    exact = {}
    for c, comp in enumerate(complex_comp):
        for m, coeff in comp.items():
            i = met_idx[m]
            for k in range(A.indptr[c], A.indptr[c + 1]):
                j = A.indices[k]
                exact[(i, j)] = exact.get((i, j), Fraction(0)) + Fraction(int(A.data[k])) * coeff

    N_exact = exact_matrix(((i, j, v) for (i, j), v in exact.items() if v),
                           (len(met_ids), len(rxn_ids)))

    complexes = [tuple(sorted(comp.items())) for comp in complex_comp]

    class CRN:
        pass
    out = CRN()
    out.A, out.Y, out.N, out.N_exact = A, Y, N, N_exact
    out.complexes = complexes
    out.reaction_ids = rxn_ids
    out.metabolite_ids = met_ids
    return out


def substrate_complexes(A):
    """The substrate complex of every reaction, {reaction: complex}."""
    A = sparse.csc_matrix(A)
    return {j: int(A.indices[k]) for j in range(A.shape[1])
            for k in range(A.indptr[j], A.indptr[j + 1]) if A.data[k] < 0}


def zero_complexes(Y):
    """Complexes with no species: the zero complex of uptake and export reactions."""
    Y = sparse.csc_matrix(Y)
    return {c for c in range(Y.shape[1]) if Y.indptr[c] == Y.indptr[c + 1]}

