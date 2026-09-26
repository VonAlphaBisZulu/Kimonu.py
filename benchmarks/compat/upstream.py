"""Kinetic modules from the concordance partition, as the published tools compute them.

Kept for comparability with COCOA.jl and the reference implementation; Kimonu itself uses the
psi null space (src/kimonu/kinetic.py). Selected in run_kimonu_crn.py with --route upstream.

A kinetic module is a set of complexes whose reactions are kinetically coupled. Following Langary
et al. (2025), each concordance module is taken together with the balanced complexes and reduced to
its *upstream set* on the complex graph — one vertex per complex, one edge per reaction, from its
substrate complex to its product complex. The upstream algorithm (Remark S2-1) has two phases:

1. remove entry complexes — complexes produced by a reaction whose substrate lies outside the set —
   repeatedly, which leaves a set nothing outside can feed;
2. remove terminal strongly connected components, which leaves a set that feeds itself.

The upstream sets are then merged, and the merging is where most of the work is. Three rules apply,
and all three are needed: sets that share a complex merge; sets whose representatives are *coupled*
merge (Proposition S3-4); and if the mass-action deficiency is one, every unbalanced concordance
module merges (Lemma S4-5). Merges feed back into the concordance modules, so the whole thing is
iterated to a fixed point, with metabolites found to have ACR augmenting the coupling matrix as they
are discovered (Remark S3-6).

The feedback matters. On EnvZ-OmpR the coupling rule alone cannot merge the set containing XD with
the set containing XD + Yp: doing so requires e_Yp to lie in the span of the coupling matrix, which
is the very question the merge would settle. The deficiency rule breaks that circle.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np

KineticResult = namedtuple(
    "KineticResult",
    ["modules", "module_of", "n_modules", "giant_module", "known_acr", "iterations", "deficiency"],
)


def reaction_complexes(A):
    """Per reaction, its substrate and product complexes."""
    A = A.tocsc()
    out = []
    for j in range(A.shape[1]):
        start, end = A.indptr[j], A.indptr[j + 1]
        rows, vals = A.indices[start:end], A.data[start:end]
        out.append(([int(r) for r, v in zip(rows, vals) if v < 0],
                    [int(r) for r, v in zip(rows, vals) if v > 0]))
    return out


class _Topology:
    """Reaction adjacency in the form the upstream algorithm needs, built once."""

    __slots__ = ("succ", "pred", "producers", "substrates", "products", "n")

    def __init__(self, A):
        rxns = reaction_complexes(A)
        n = A.shape[0]
        self.n = n
        self.succ = [set() for _ in range(n)]
        self.pred = [set() for _ in range(n)]
        self.producers = [[] for _ in range(n)]      # complex -> reactions producing it
        self.substrates = [s for s, _ in rxns]       # reaction -> substrate complexes
        self.products = [p for _, p in rxns]         # reaction -> product complexes
        for j, (subs, prods) in enumerate(rxns):
            for p in prods:
                self.producers[p].append(j)
            for sc in subs:
                for p in prods:
                    self.succ[sc].add(p)
                    self.pred[p].add(sc)


def complex_graph(A):
    """Successors and predecessors of each complex in the reaction graph."""
    n = A.shape[0]
    succ = [set() for _ in range(n)]
    pred = [set() for _ in range(n)]
    for substrates, products in reaction_complexes(A):
        for s in substrates:
            for p in products:
                succ[s].add(p)
                pred[p].add(s)
    return succ, pred


def _sccs(succ, scope):
    """Tarjan's strongly connected components restricted to ``scope``."""
    index, low, on_stack, stack, result = {}, {}, set(), [], []
    counter = 0
    for root in sorted(scope):
        if root in index:
            continue
        work = [(root, iter(sorted(succ[root] & scope)))]
        index[root] = low[root] = counter
        counter += 1
        stack.append(root)
        on_stack.add(root)
        while work:
            v, it = work[-1]
            descended = False
            for w in it:
                if w not in index:
                    index[w] = low[w] = counter
                    counter += 1
                    stack.append(w)
                    on_stack.add(w)
                    work.append((w, iter(sorted(succ[w] & scope))))
                    descended = True
                    break
                if w in on_stack:
                    low[v] = min(low[v], index[w])
            if descended:
                continue
            work.pop()
            if work:
                low[work[-1][0]] = min(low[work[-1][0]], low[v])
            if low[v] == index[v]:
                component = []
                while True:
                    w = stack.pop()
                    on_stack.discard(w)
                    component.append(w)
                    if w == v:
                        break
                result.append(set(component))
    return result


def upstream_set(candidates, topology, succ=None):
    """The upstream set of a candidate complex set (Remark S2-1)."""
    current = set(candidates)
    succ = topology.succ if succ is None else succ

    # phase I: drop entry complexes, to a fixed point
    #
    # A complex is an entry complex when some reaction producing it draws on a substrate outside
    # the set.
    substrates = topology.substrates
    producers = topology.producers
    while True:
        entry = {c for c in current
                 if any(any(sub not in current for sub in substrates[j]) for j in producers[c])}
        if not entry:
            break
        current -= entry

    # phase II: drop terminal strongly connected components
    #
    # Terminality is decided against the set as it stands after phase I, and every terminal
    # component is removed together. Removing them one at a time and re-testing against the
    # shrinking set cascades: each removal can strand the component that fed it, and the whole set
    # unravels.
    if current:
        components = _sccs(succ, current)
        # Terminality is judged on the full network, not on the candidate set: a component that
        # feeds something outside the candidate is still feeding, so it is not terminal. Testing
        # only against `current` calls such a component terminal and deletes it, which strands its
        # complexes in singleton modules and loses every ACR/ACRR relation among them.
        terminal = [c for c in components if not any(succ[v] - c for v in c)]
        for component in terminal:
            current -= component
    return current


def _weak_components(succ, pred, scope):
    seen, out = set(), []
    for start in sorted(scope):
        if start in seen:
            continue
        stack, comp = [start], []
        seen.add(start)
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in (succ[v] | pred[v]) & scope:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        out.append(set(comp))
    return out


def is_weakly_reversible(succ, pred, n):
    """True when every linkage class is a single strongly connected component."""
    scope = set(range(n))
    for component in _weak_components(succ, pred, scope):
        if len(_sccs(succ, component)) > 1:
            return False
    return True


def classical_deficiency(Y, A, succ, pred, nullity=None):
    """delta = n - l - rank(N): complexes, minus linkage classes, minus the rank of N = Y*A.

    The rank comes from the exact nullspace, not from an SVD. A float rank test would put a
    numerical tolerance in the middle of a discrete rule — the deficiency decides whether every
    unbalanced module merges, and the difference between delta = 1 and delta = 2 changes the answer.
    It is also the difference between an O(n^3) dense factorisation and a sparse exact one.

    ``nullity`` may be supplied when the kernel has already been computed for the concordance
    partition, in which case no factorisation happens here at all.
    """
    n = A.shape[0]
    linkage = len(_weak_components(succ, pred, set(range(n))))
    if nullity is None:
        import straindesign as sd
        nullity = sd.sparse_nullspace((Y @ A).tocsr()).shape[1]
    rank = A.shape[1] - nullity
    return n - linkage - rank


def kinetic_modules(A, modules, balanced, Y=None, max_iterations=20, nullity=None):
    """Kinetic modules from the concordance partition.

    Args:
        A: complexes x reactions incidence matrix.
        modules: concordance modules as iterables of complex indices (singletons included).
        balanced: indices of balanced complexes.
        Y: metabolites x complexes composition matrix. Required for the coupling and deficiency
            merges; without it only overlapping upstream sets are merged, which reproduces the
            reference implementation's behaviour and can miss ACR.
        max_iterations: cap on the outer loop; convergence is normally reached in a few passes.
        nullity: dimension of ker(N), when already known from the concordance step. Supplying it
            skips the only factorisation this function would otherwise perform.

    Returns a :class:`KineticResult`.
    """
    from kimonu.acr import acr_acrr, coupling_signatures

    n = A.shape[0]
    topology = _Topology(A)
    succ, pred = topology.succ, topology.pred
    balanced = set(balanced)

    concordance = [set(m) - balanced for m in modules]
    concordance = [m for m in concordance if m]
    initial_delta = (classical_deficiency(Y, A, succ, pred, nullity=nullity)
                     if Y is not None else -1)
    weakly_reversible = is_weakly_reversible(succ, pred, n)

    known_acr = set()
    n_merges = 0
    previous = None
    iterations = 0
    current_modules = []

    while iterations < max_iterations:
        iterations += 1
        upstream_sets = []
        origin = []
        for idx, module in enumerate(concordance):
            s = upstream_set(module | balanced, topology)
            if s:
                upstream_sets.append(s)
                origin.append(idx)

        if not upstream_sets:
            break

        sets = _group_sets(upstream_sets, origin, Y, known_acr)

        # a merge of upstream sets is a merge of the concordance modules behind them
        new_concordance = []
        for _, sources in sets:
            union = set()
            for i in sources:
                union |= concordance[i]
            new_concordance.append(union)
        for idx, module in enumerate(concordance):
            if not any(idx in sources for _, sources in sets):
                new_concordance.append(module)

        n_merges += len(concordance) - len(new_concordance)
        current_modules = [s for s, _ in sets]

        if Y is not None:
            found = acr_acrr(Y, current_modules, known_acr=known_acr)
            known_acr |= set(found.acr)

        state = (frozenset(frozenset(s) for s in current_modules),
                 frozenset(frozenset(m) for m in new_concordance))
        concordance = new_concordance
        if state == previous:
            break
        previous = state

    # Lemma S4-5: with mass-action deficiency one, every unbalanced module belongs together
    deficiency_merged = False
    if Y is not None and initial_delta >= 0 and len(concordance) > 1:
        lower = 0 if weakly_reversible else 1
        upper = initial_delta - n_merges
        if lower == 1 and upper == 1:
            concordance = [set().union(*concordance)]
            deficiency_merged = True
            upstream_sets = []
            for module in concordance:
                s = upstream_set(module | balanced, topology)
                if s:
                    upstream_sets.append(s)
            sets = _group_sets(upstream_sets, list(range(len(upstream_sets))), Y, known_acr)
            current_modules = [s for s, _ in sets]

    result = list(current_modules)

    # linkage classes made up entirely of balanced complexes stand on their own
    assigned = set().union(*result) if result else set()
    for comp in _weak_components(succ, pred, set(range(n))):
        if comp <= balanced and not (comp & assigned):
            result.append(comp)
            assigned |= comp

    for i in range(n):
        if i not in assigned:
            result.append({i})

    module_of = np.zeros(n, dtype=np.int64)
    for idx, group in enumerate(result, start=1):
        for i in group:
            module_of[i] = idx

    return KineticResult(
        modules=result,
        module_of=module_of,
        n_modules=sum(1 for m in result if len(m) > 1),
        giant_module=max((len(m) for m in result), default=0),
        known_acr=sorted(known_acr),
        iterations=iterations,
        deficiency=initial_delta,
    )


def _group_sets(upstream_sets, origin, Y, known_acr):
    """Merge upstream sets by overlap and, when Y is given, by coupling (Proposition S3-4).

    Returns a list of (merged set, set of originating concordance module indices).
    """
    from kimonu.acr import coupling_signatures

    items = [(set(s), {o}) for s, o in zip(upstream_sets, origin)]

    # overlap
    merged = []
    for group, sources in items:
        touching = [x for x in merged if x[0] & group]
        for x in touching:
            merged.remove(x)
            group = group | x[0]
            sources = sources | x[1]
        merged.append((group, sources))

    if Y is None or len(merged) < 2:
        return merged

    signatures = coupling_signatures(Y, [g for g, _ in merged], known_acr)
    if signatures is None:
        return merged

    by_signature = {}
    for (group, sources), sig in zip(merged, signatures):
        key = by_signature.setdefault(sig, [set(), set()])
        key[0] |= group
        key[1] |= sources
    return [(g, s) for g, s in by_signature.values()]
