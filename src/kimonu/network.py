"""Networks and preprocessing: from a cobra model, or from given matrices, to the network N = Y·A.

A metabolic network factorizes as N = Y·A:
  N : metabolites x reactions   stoichiometric matrix
  Y : metabolites x complexes   composition (which metabolites, and how many, form each complex)
  A : complexes  x reactions    incidence (-1 on a reaction's substrate complex, +1 on its product)

A *complex* is one side of a reaction, a multiset of metabolites such as D + E, 2B or C. Column c of Y is its
composition y_c. Uptake and export reactions have the zero complex on one side, a complex with no metabolites (an
all-zero column of Y, as in COCOA's `:empty`). It gets its own row of A, so every reaction has a substrate and a
product complex. Mass action needs that, because an uptake from the zero complex runs at rate k_j ψ_∅ with ψ_∅ = 1.
Complexes are numbered in the order in which the reactions first use them, as in COCOA.

Mass action gives every reaction one substrate complex and one direction, so a reaction that can run both ways
becomes two (:func:`make_irreversible`). Each reaction is otherwise taken as one mass-action step. Kimonu does not
split enzyme-catalyzed reactions into elementary steps, and networks that need them have to carry them already.

The network is a :class:`CRN`, and every analysis takes it as its first argument. :func:`preprocess` builds it from a
cobra model with the preprocessing every analysis shares. It splits reversible reactions, removes blocked
reactions, finds pinned rates and builds the factorization (:func:`build_crn`). Blocked reactions follow from one
LP (:func:`_blocked_reactions`) when every lower bound is zero, and from FVA only when a positive lower bound forces
flux, since only then can a bound pin a rate. :func:`crn_from_matrices` builds the CRN of a network that already
comes as A and Y, such as the elementary-step models of Langary et al. (2025). Every analysis also takes a cobra
model and preprocesses it first (:func:`as_crn`). When several analyses run on one model, preprocess it once and
pass the CRN to each.

Detached uptakes. With detached uptakes (the default of every analysis), each reaction out of the zero complex draws
from a source of its own, whose concentration is free. :func:`uptake_columns` gives the source of each such
reaction, and :func:`with_sources` extends the network by one coordinate per source, for the computations that need
the sources as coordinates.
"""
from collections import namedtuple

import cobra

import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from scipy import sparse
from scipy.optimize import linprog

import straindesign as sd

from kimonu.linalg import exact, fraction, matrix_product, row_form, to_float

CRN = namedtuple("CRN", ["A", "Yt", "N", "complexes", "reaction_ids", "metabolite_ids", "pins"])
CRN.__doc__ = """One network, N = Y·A, with the names of the rows and columns. Every analysis takes it first.

A               scipy.sparse csr, int64 (complexes x reactions), the incidence matrix
Yt              Yᵀ (complexes x metabolites) in the row form of :mod:`kimonu.linalg`, exact. Row c is the composition
                y_c of complex c, a sorted tuple of (metabolite, coefficient) pairs, empty for the zero complex. Every
                computation reads Y one complex at a time, Yᵀ B as matrix_product(Yt, B) and y_c - y_d from two rows.
N               straindesign RationalMatrix (metabolites x reactions), the stoichiometric matrix N = Y·A, exact.
                Every nullspace of N is computed from it.
complexes       one complex per row of A, as a sorted tuple of (metabolite id, coefficient) pairs, for reading
reaction_ids    reaction ids in the column order of A and N
metabolite_ids  metabolite ids in the row order of N. Results name metabolites by them.
pins            {reaction id: exact value} of the rates that a bound holds at a nonzero value in every steady
                state. Empty from :func:`build_crn` and :func:`crn_from_matrices`, found by :func:`preprocess`.
"""


def _crn(A, Yt, metabolite_ids, reaction_ids):
    """The CRN of the incidence matrix A and the compositions Yt, with N = Y·A computed exactly.

    Column j of N is y_product(j) - y_substrate(j), the rows of Yt weighted by column j of A, so the rows of Nᵀ are
    Aᵀ Yᵀ (:func:`kimonu.linalg.matrix_product`).
    """
    Nt = matrix_product(A.T, Yt)
    N = exact(((i, j, v) for j, row in enumerate(Nt) for i, v in row), (len(metabolite_ids), len(reaction_ids)))
    complexes = [tuple(sorted((metabolite_ids[i], v) for i, v in y)) for y in Yt]
    return CRN(A=A, Yt=Yt, N=N, complexes=complexes, reaction_ids=list(reaction_ids),
               metabolite_ids=list(metabolite_ids), pins={})


def crn_from_matrices(A, Y, metabolite_ids=None, reaction_ids=None):
    """The :class:`CRN` of a network that comes as its incidence matrix ``A`` and composition ``Y``.

    For networks that already carry a complex decomposition, such as the elementary-step models distributed with the
    reference implementation of Langary et al. (2025), there is nothing to enumerate, and rebuilding a cobra model
    only to take it apart again would risk changing the decomposition being analysed.

    Args:
        A: complexes x reactions incidence matrix, one -1 (the substrate complex) and one +1 (the product complex)
            per column, as any scipy or numpy matrix.
        Y: metabolites x complexes composition matrix, as any scipy or numpy matrix. Every entry is read exactly, a
            float as the decimal it prints as (:func:`kimonu.linalg.fraction`).
        metabolite_ids, reaction_ids: names of the rows of Y and of the columns of A. Indices when omitted, so
            results then name metabolites by their index.

    The network is taken as it is. It must have no blocked reactions for the mass-action routes and for exact
    concordance, as :func:`preprocess` ensures for a cobra model.
    """
    A = sparse.csr_matrix(A).astype(np.int64)
    Yt = row_form(sparse.csr_matrix(Y).T)
    m = sparse.csr_matrix(Y).shape[0]
    return _crn(A, Yt, metabolite_ids if metabolite_ids is not None else list(range(m)),
                reaction_ids if reaction_ids is not None else list(range(A.shape[1])))


def make_irreversible(model):
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
    """The :class:`CRN` of ``model`` after :func:`make_irreversible`, with no reaction removed and no pins.

    Use :func:`preprocess` for a model that may have blocked reactions or pinned rates.
    """
    model = make_irreversible(model)
    met_ids = [m.id for m in model.metabolites]
    met_idx = {m: i for i, m in enumerate(met_ids)}

    complex_index = {}   # composition (sorted (metabolite index, coefficient) pairs) -> complex row index

    def _complex(side):
        """The row index of the complex with composition {metabolite: coefficient > 0}; {} is the zero complex."""
        y = tuple(sorted((met_idx[m.id], fraction(abs(c))) for m, c in side.items()))
        return complex_index.setdefault(y, len(complex_index))

    # A: one substrate and one product complex per reaction, numbered in order of first use
    Ai, Aj, Av = [], [], []
    for j, r in enumerate(model.reactions):
        Ai += [_complex({m: c for m, c in r.metabolites.items() if c < 0}),
               _complex({m: c for m, c in r.metabolites.items() if c > 0})]
        Aj += [j, j]
        Av += [-1, 1]
    A = sparse.csr_matrix((Av, (Ai, Aj)), shape=(len(complex_index), len(model.reactions)), dtype=np.int64)
    # Yᵀ: the composition of every complex, in the order of the rows of A
    return _crn(A, list(complex_index), met_ids, [r.id for r in model.reactions])


def preprocess(model, solver=None):
    """The :class:`CRN` of ``model`` ready for analysis: split, without blocked reactions, with its pinned rates.

    1. Reversible reactions are split into two (:func:`make_irreversible`).
    2. Blocked reactions are found and removed, together with metabolites that no reaction uses any more.
    3. Rates that a bound holds at a nonzero value in every steady state (pins) are found.
    4. The factorization N = Y·A is built (:func:`build_crn`), with the pins attached.

    Concordance is exact over the feasible steady states only when no reaction is blocked and the pins enter N.
    :func:`kimonu.concordance` adds them as rows that fix the ratio of every two pinned rates. The mass-action routes
    need a network without blocked reactions, since they start from a strictly positive steady-state flux, and do
    not use the pins, since under mass action the rate constants, not the bounds, set the rates.

    Args:
        model: a cobra model. It is not changed.
        solver: solver for FVA, which is needed only when a positive lower bound forces flux; None lets
            straindesign choose.
    """
    model = make_irreversible(model).copy()
    blocked, pins = _find_blocked_and_pinned(model, solver)
    if blocked:
        model.remove_reactions(blocked, remove_orphans=True)
    return build_crn(model)._replace(pins=pins)


_preprocess = preprocess          # as_crn has a keyword of the same name


def as_crn(network, preprocess=True, solver=None):
    """The :class:`CRN` of ``network``, which every analysis calls on its first argument.

    Args:
        network: a :class:`CRN`, which is returned as it is, or a cobra model.
        preprocess: for a cobra model, run :func:`preprocess` (default), which splits reversible reactions, removes
            blocked reactions and finds pinned rates. A model that has no blocked reactions and no pinned rates can
            be passed with ``preprocess=False``, which only builds N = Y·A (:func:`build_crn`).
        solver: solver for FVA, if the model needs it (:func:`preprocess`). None lets straindesign choose.

    A cobra model is preprocessed anew in every call. When several analyses run on one model, call
    :func:`preprocess` once and pass the CRN to each.
    """
    if isinstance(network, CRN):
        return network
    return _preprocess(network, solver) if preprocess else build_crn(network)


def _find_blocked_and_pinned(model, solver=None):
    """(ids of the blocked reactions, pinned rates) of an irreversible model.

    With every lower bound at zero, r = 0 is a steady state, so a rate that is constant over the steady states is
    zero. Nothing is pinned then, and finite upper bounds only scale, so they cannot block a reaction. One LP then
    settles it (:func:`_blocked_reactions`). A positive lower bound forces flux and needs FVA.
    """
    if all(r.lower_bound <= 0 for r in model.reactions):
        S = create_stoichiometric_matrix(model, array_type="lil")
        closed = [j for j, r in enumerate(model.reactions) if r.upper_bound <= 0]
        return [model.reactions[j].id for j in np.flatnonzero(_blocked_reactions(S.tocsr(), closed))], {}
    ranges = sd.fva(model, solver=solver) if solver else sd.fva(model)
    blocked = [r for r in ranges.index
               if abs(ranges.at[r, "minimum"]) <= 1e-9 and abs(ranges.at[r, "maximum"]) <= 1e-9]
    return blocked, _find_pins(model, ranges)


def _blocked_reactions(N, closed=()):
    """Boolean array, True for the reactions of N that no steady state with r >= 0 can run.

    One LP: max sum(z) subject to N r = 0, r >= 0, 0 <= z <= 1 and z <= r. The steady states form a cone,
    and two steady states add up to one that runs the reactions of both, so a single r runs every reaction
    that can run at all. Scaling it brings each of them to at least 1, which makes z_j = 1 for those and
    z_j = 0 for the blocked ones. Maximizing sum(r) under upper bounds instead stops at a vertex and leaves
    reactions that can run at zero. ``closed`` lists the indices of reactions with an upper bound of zero.
    """
    Nf = to_float(N); m, q = Nf.shape
    closed = set(closed)
    bounds = [(0, 0 if j in closed else None) for j in range(q)] + [(0, 1)] * q
    A_eq = sparse.hstack([Nf, sparse.csr_matrix((m, q))], format="csr")
    A_ub = sparse.hstack([-sparse.identity(q), sparse.identity(q)], format="csr")
    lp = linprog(np.r_[np.zeros(q), -np.ones(q)], A_ub=A_ub, b_ub=np.zeros(q), A_eq=A_eq, b_eq=np.zeros(m),
                 bounds=bounds, method="highs")
    if lp.status != 0:
        raise ValueError(f"blocked-reaction LP failed: {lp.message}")
    return lp.x[q:] <= 0.5


def _find_pins(model, ranges, tol=1e-9):
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
                pins[r.id] = fraction(b)
                break
    return pins


def substrate_complexes(crn):
    """The substrate complex σ(j) of every reaction j, as {reaction: complex}. It defines S, S_cj = 1 when c = σ(j)."""
    A = sparse.csc_matrix(crn.A)
    return {j: int(A.indices[k]) for j in range(A.shape[1])
            for k in range(A.indptr[j], A.indptr[j + 1]) if A.data[k] < 0}


def zero_complexes(crn):
    """Complexes with no metabolites, the zero complex of uptake and export reactions."""
    return {c for c, y in enumerate(crn.Yt) if not y}


def uptake_columns(crn, detach_uptakes=True):
    """The column of u that carries each reaction's rate, {reaction: column}, and the number of columns.

    Reaction j runs at r*_j u_c for its column c. That is its substrate complex σ(j), or, for a detached uptake, a
    source of its own after the n complexes, numbered in the order of the reactions. These are the columns of
    W (:func:`kimonu.steady_state_map`).
    """
    n = len(crn.Yt)
    zero = zero_complexes(crn) if detach_uptakes else set()
    column, n_columns = {}, n
    for j, c in substrate_complexes(crn).items():
        if c in zero:
            column[j], n_columns = n_columns, n_columns + 1
        else:
            column[j] = c
    return column, n_columns


def with_sources(crn, detach_uptakes=True):
    """The network with every detached source as a species of its own, for the routes that give sources coordinates.

    A detached source s (column n + s of W, :func:`uptake_columns`) has a free concentration. As a species of its
    own, with composition e_(m+s), it is one more row of Yᵀ after the n complexes, the Ŷᵀ that the tangent route
    multiplies W by (W Ŷᵀ) and flux coupling reads (Ŷᵀ B), and one more coordinate after the m metabolites, named
    "source of <reaction id>". Only Yt and metabolite_ids are extended, since a source is not balanced and adds no
    row to N. With attached uptakes, or without uptakes, the network is returned as it is.
    """
    column, n_columns = uptake_columns(crn, detach_uptakes)
    m, n = len(crn.metabolite_ids), len(crn.Yt)
    if n_columns == n:
        return crn
    uptakes = sorted((c, j) for j, c in column.items() if c >= n)
    return crn._replace(Yt=crn.Yt + [((m + s, 1),) for s in range(n_columns - n)],
                        metabolite_ids=list(crn.metabolite_ids) + [f"source of {crn.reaction_ids[j]}"
                                                                   for _, j in uptakes])
