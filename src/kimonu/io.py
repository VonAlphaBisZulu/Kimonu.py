"""Networks and preprocessing: from a cobra model to the matrices N = Y·A.

A metabolic network factorizes as N = Y·A:
  N : metabolites x reactions   stoichiometric matrix
  Y : metabolites x complexes   composition (which metabolites, and how many, form each complex)
  A : complexes  x reactions    incidence (-1 on a reaction's substrate complex, +1 on its product)

A *complex* is one side of a reaction, a multiset of metabolites such as D + E, 2B or C. Uptake and export
reactions have the zero complex on one side, a complex with no metabolites (an all-zero column of Y, as in COCOA's
`:empty`). It gets its own row of A, so every reaction has a substrate and a product complex. Mass action needs
that, because an uptake from the zero complex runs at rate k_j ψ_∅ with ψ_∅ = 1. Complexes are numbered in the
order in which the reactions first use them, as in COCOA.

Mass action gives every reaction one substrate complex and one direction, so a reaction that can run both ways
becomes two (:func:`irreversible`). Each reaction is otherwise taken as one mass-action step. Kimonu does not split
enzyme-catalyzed reactions into elementary steps, and networks that need them have to carry them already.

:func:`prepare` is the preprocessing every analysis of a cobra model shares. It splits reversible reactions,
removes blocked reactions, finds pinned rates and builds the factorization (:func:`build_crn`). Blocked reactions
follow from one LP (:func:`kimonu.linalg.blocked_reactions`) when every lower bound is zero, and from FVA only when
a positive lower bound forces flux, since only then can a bound pin a rate.
"""
from collections import namedtuple
from fractions import Fraction

import cobra

import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from scipy import sparse

import straindesign as sd

from kimonu.linalg import blocked_reactions

CRN = namedtuple("CRN", ["A", "Y", "N", "N_float", "complexes", "reaction_ids", "metabolite_ids", "pins"])
CRN.__doc__ = """The matrices of N = Y·A for one network, with the names of their rows and columns.

A               scipy.sparse csr, int64 (complexes x reactions), the incidence matrix
Y               scipy.sparse csr, float64 (metabolites x complexes), the composition. scipy.sparse cannot hold
                Fractions, so exact readers convert each entry with Fraction(str(v)), which recovers the decimal
                the model declares. ``complexes`` holds the same compositions as Fractions.
N               straindesign RationalMatrix (metabolites x reactions), the stoichiometric matrix with exact
                entries. Every nullspace is computed from it.
N_float         scipy.sparse csr, float64, N = Y @ A in floating point, for inspection
complexes       one canonical complex per row of A, a tuple of (metabolite id, Fraction coefficient)
reaction_ids    reaction ids in the column order of A and N
metabolite_ids  metabolite ids in the row order of Y and N
pins            {reaction id: exact value} of the rates that a bound holds at a nonzero value in every steady
                state. Empty from :func:`build_crn`, found by :func:`prepare`.
"""


def _exact(coeff):
    """Exact rational value of a stoichiometric coefficient.

    ``Fraction(str(x))`` reproduces the decimal the model declares (1.496 -> 1496/1000), whereas ``Fraction(x)``
    would recover the binary approximation of it. Coefficients are frequently fractional, as in biomass reactions
    and half-oxygen terms such as 0.5 O2 in CYTBD, so this is the difference between analysing the model and
    analysing a rounded version of it.
    """
    return Fraction(str(coeff))


def exact_matrix(entries, shape):
    """A straindesign RationalMatrix from ``(row, col, value)`` triples, holding them exactly.

    Clearing denominators to int64 instead overflows on models whose coefficients are written as
    15-significant-digit decimals. iYS854's -7.73333333333333e-08 has a denominator of 10**22, and
    Chlorella_variabilis_iAJ526 mixes seven such terms in one row.
    """
    return sd.RationalMatrix.from_fractions(
        ((int(i), int(j), v if isinstance(v, Fraction) else Fraction(str(v))) for i, j, v in entries),
        shape)


def irreversible(model):
    """``model`` with every reaction running in its written direction only.

    A reaction with a negative lower bound gets a reverse reaction ``<id>_rev`` with the opposite stoichiometry and
    mirrored bounds. It keeps its own forward direction if its upper bound is positive and is replaced by the
    reverse reaction otherwise. The input is not changed. A model that needs no split is returned as it is, any
    other as a changed copy.
    """
    if all(r.lower_bound >= 0 for r in model.reactions):
        return model
    out = model.copy()
    for r in list(out.reactions):
        if r.lower_bound >= 0:
            continue
        reverse = cobra.Reaction(r.id + "_rev", lower_bound=max(0, -r.upper_bound), upper_bound=-r.lower_bound)
        out.add_reactions([reverse])
        reverse.add_metabolites({m: -c for m, c in r.metabolites.items()})
        if r.upper_bound > 0:
            r.lower_bound = 0
        else:
            out.remove_reactions([r])
    return out


def build_crn(model):
    """The :class:`CRN` of ``model`` after :func:`irreversible`, with no reaction removed and no pins.

    Use :func:`prepare` for a model that may have blocked reactions or pinned rates.
    """
    model = irreversible(model)
    met_ids = [m.id for m in model.metabolites]
    met_idx = {m: i for i, m in enumerate(met_ids)}
    rxn_ids = [r.id for r in model.reactions]

    complex_index = {}   # canonical key -> complex row index
    complex_comp = []    # complex row index -> {met_id: coeff (>0)}

    def _complex(comp):
        """Return the complex index for a multiset {met_id: coeff>0}; {} is the zero complex."""
        key = tuple(sorted(comp.items()))
        idx = complex_index.get(key)
        if idx is None:
            idx = len(complex_comp)
            complex_index[key] = idx
            complex_comp.append(comp)
        return idx

    # 1) A: one substrate and one product complex per reaction
    Ai, Aj, Av = [], [], []
    for j, r in enumerate(model.reactions):
        substrate = {m.id: -_exact(c) for m, c in r.metabolites.items() if c < 0}   # |coeff| on the substrate side
        product = {m.id: _exact(c) for m, c in r.metabolites.items() if c > 0}      # coeff on the product side
        si = _complex(substrate)
        pi = _complex(product)
        Ai.append(si); Aj.append(j); Av.append(-1)
        Ai.append(pi); Aj.append(j); Av.append(+1)

    ncplx = len(complex_comp)
    A = sparse.csr_matrix((Av, (Ai, Aj)), shape=(ncplx, len(rxn_ids)), dtype=np.int64)

    # 2) Y: the composition of every complex
    Yi, Yj, Yv = [], [], []
    for c, comp in enumerate(complex_comp):
        for m, coeff in comp.items():
            Yi.append(met_idx[m]); Yj.append(c); Yv.append(coeff)
    Y = sparse.csr_matrix(([float(v) for v in Yv], (Yi, Yj)),
                          shape=(len(met_ids), ncplx), dtype=np.float64)

    N_float = (Y @ A).tocsr()  # == stoichiometric matrix, by construction

    # 3) N = Y·A exactly, accumulated from the rational compositions
    exact = {}
    for c, comp in enumerate(complex_comp):
        for m, coeff in comp.items():
            i = met_idx[m]
            for k in range(A.indptr[c], A.indptr[c + 1]):
                j = A.indices[k]
                exact[(i, j)] = exact.get((i, j), Fraction(0)) + Fraction(int(A.data[k])) * coeff

    N = exact_matrix(((i, j, v) for (i, j), v in exact.items() if v), (len(met_ids), len(rxn_ids)))

    return CRN(A=A, Y=Y, N=N, N_float=N_float, complexes=[tuple(sorted(comp.items())) for comp in complex_comp],
               reaction_ids=rxn_ids, metabolite_ids=met_ids, pins={})


def prepare(model, solver=None):
    """The :class:`CRN` of ``model`` ready for analysis: split, without blocked reactions, with its pinned rates.

    1. Reversible reactions are split into two (:func:`irreversible`).
    2. Blocked reactions are found and removed, together with metabolites that no reaction uses any more.
    3. Rates that a bound holds at a nonzero value in every steady state (pins) are found.
    4. The factorization N = Y·A is built (:func:`build_crn`), with the pins attached.

    Concordance is exact over the feasible steady states only when no reaction is blocked and the pins enter N as
    rows (:func:`kimonu.concordance` adds them). The mass-action routes need a network without blocked reactions,
    since they start from a strictly positive steady-state flux, and do not use the pins, since under mass action
    the rate constants, not the bounds, set the rates.

    Args:
        model: a cobra model. It is not changed.
        solver: solver for FVA, which is needed only when a positive lower bound forces flux; None lets
            straindesign choose.
    """
    model = irreversible(model).copy()
    blocked, pins = _blocked_and_pinned(model, solver)
    if blocked:
        model.remove_reactions(blocked, remove_orphans=True)
    return build_crn(model)._replace(pins=pins)


def _blocked_and_pinned(model, solver=None):
    """(ids of the blocked reactions, pinned rates) of an irreversible model.

    With every lower bound at zero, r = 0 is a steady state, so a rate that is constant over the steady states is
    zero. Nothing is pinned then, and finite upper bounds only scale, so they cannot block a reaction. One LP then
    settles it. A positive lower bound forces flux and needs FVA.
    """
    if all(r.lower_bound <= 0 for r in model.reactions):
        S = create_stoichiometric_matrix(model, array_type="lil")
        closed = [j for j, r in enumerate(model.reactions) if r.upper_bound <= 0]
        return [model.reactions[j].id for j in np.flatnonzero(blocked_reactions(S.tocsr(), closed))], {}
    ranges = sd.fva(model, solver=solver) if solver else sd.fva(model)
    blocked = [r for r in ranges.index
               if abs(ranges.at[r, "minimum"]) <= 1e-9 and abs(ranges.at[r, "maximum"]) <= 1e-9]
    return blocked, _pins(model, ranges)


def _pins(model, ranges, tol=1e-9):
    """{reaction id: exact value} for the rates that a bound holds at a nonzero value in every steady state.

    FVA shows such a rate as minimum == maximum. The bound it sits on gives the value exactly. A rate that is
    constant only because others are pinned follows from those and needs no entry.
    """
    pins = {}
    for r in model.reactions:
        lo, hi = ranges.at[r.id, "minimum"], ranges.at[r.id, "maximum"]
        if abs(hi) <= tol or hi - lo > tol * max(1, abs(hi)):
            continue
        for b in (r.lower_bound, r.upper_bound):
            if abs(hi - b) <= tol * max(1, abs(b)):
                pins[r.id] = Fraction(str(b))
                break
    return pins


def substrate_complexes(A):
    """The substrate complex σ(j) of every reaction j, as {reaction: complex}. It defines S, S_cj = 1 when c = σ(j)."""
    A = sparse.csc_matrix(A)
    return {j: int(A.indices[k]) for j in range(A.shape[1])
            for k in range(A.indptr[j], A.indptr[j + 1]) if A.data[k] < 0}


def zero_complexes(Y):
    """Complexes with no metabolites, the zero complex of uptake and export reactions."""
    Y = sparse.csc_matrix(Y)
    return {c for c in range(Y.shape[1]) if Y.indptr[c] == Y.indptr[c + 1]}
