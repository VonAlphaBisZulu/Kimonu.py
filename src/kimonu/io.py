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
import numpy as np
from scipy import sparse


def build_crn(model):
    """Return a namespace with A, Y, N and the id lists.

    Returns an object with attributes:
      A            scipy.csr int  (complexes x reactions), incidence
      Y            scipy.csr int  (metabolites x complexes), composition
      N            scipy.csr      (metabolites x reactions) = Y @ A
      complexes    list of canonical complexes (tuple of (met_id, coeff)), in A's row order
      reaction_ids, metabolite_ids
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
        substrate = {m.id: -c for m, c in r.metabolites.items() if c < 0}   # |coeff| on the substrate side
        product = {m.id: c for m, c in r.metabolites.items() if c > 0}      # coeff on the product side
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
    Y = sparse.csr_matrix((Yv, (Yi, Yj)), shape=(len(met_ids), ncplx), dtype=np.int64)

    N = (Y @ A).tocsr()  # == stoichiometric matrix, by construction

    complexes = [tuple(sorted(comp.items())) for comp in complex_comp]

    class CRN:
        pass
    out = CRN()
    out.A, out.Y, out.N = A, Y, N
    out.complexes = complexes
    out.reaction_ids = rxn_ids
    out.metabolite_ids = met_ids
    return out
