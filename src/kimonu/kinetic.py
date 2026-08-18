"""Kinetic modules from the concordance partition.

A kinetic module is a set of complexes whose reactions are kinetically coupled. Following the
upstream algorithm of Langary et al. (2025), each concordance module is tested together with the
balanced complexes for *autonomy* on the directed complex graph — one vertex per complex, one edge
per reaction, from its substrate complex to its product complex — and the surviving complexes are
what the four phases below retain.

Ported from the reference implementation's `code_kineticModule_analysis.R`. Two things there are
worth recording because they are not obvious from the paper:

* The reference deletes edges whose `Elementary_Kinetic_Param` is zero. That field is absent from
  the committed result files, so in R the expression selects nothing and the step is inert. This
  implementation keeps every edge, which reproduces the published behaviour.
* Phase IV computes an induced subgraph and then takes strongly connected components of the *full*
  graph rather than of that subgraph. Reproducing the published numbers requires the full graph, so
  that is what ``phase_iv_on_full_graph`` does by default; pass ``False`` for the reading the
  induced subgraph suggests.
"""
from collections import namedtuple
from fractions import Fraction

import numpy as np

KineticResult = namedtuple("KineticResult", ["modules", "module_of", "n_modules", "giant_module"])


def complex_graph(A):
    """Directed graph on complexes: an edge per reaction, substrate complex -> product complex.

    Returns (successors, predecessors) as lists of sets, and the edge list.
    """
    A = A.tocsc()
    n = A.shape[0]
    succ = [set() for _ in range(n)]
    pred = [set() for _ in range(n)]
    edges = []
    for j in range(A.shape[1]):
        start, end = A.indptr[j], A.indptr[j + 1]
        rows = A.indices[start:end]
        vals = A.data[start:end]
        substrates = [int(r) for r, v in zip(rows, vals) if v < 0]
        products = [int(r) for r, v in zip(rows, vals) if v > 0]
        for s in substrates:
            for p in products:
                succ[s].add(p)
                pred[p].add(s)
                edges.append((s, p))
    return succ, pred, edges


def _strongly_connected(succ, vertices=None):
    """Tarjan's strongly connected components, iterative."""
    n = len(succ)
    scope = set(range(n)) if vertices is None else set(vertices)
    index = {}
    low = {}
    on_stack = set()
    stack = []
    result = []
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
            advanced = False
            for w in it:
                if w not in index:
                    index[w] = low[w] = counter
                    counter += 1
                    stack.append(w)
                    on_stack.add(w)
                    work.append((w, iter(sorted(succ[w] & scope))))
                    advanced = True
                    break
                if w in on_stack:
                    low[v] = min(low[v], index[w])
            if advanced:
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
                result.append(component)
    return result


def _non_terminal(succ, scope, vertices):
    """Complexes in a strongly connected component that has an edge leaving it."""
    out = set()
    for component in _strongly_connected(succ, scope):
        members = set(component)
        leaves = any((succ[v] - members) for v in members)
        if leaves:
            out |= members & set(vertices)
    return out


def autonomous(candidates, succ, pred):
    """The four phases of the upstream algorithm on a candidate complex set.

    Phase I removes complexes fed from outside the set, iterated to a fixed point. Phase II removes
    complexes with no outgoing reaction. Phase III removes, again to a fixed point, complexes that
    feed outside the set; the ones removed here are kept. Phase IV keeps the non-terminal survivors.
    """
    x = set(candidates)

    # phase I: no producer outside the set
    while True:
        drop = {v for v in x if pred[v] - x}
        if not drop:
            break
        x -= drop
        if not x:
            return set()

    # phase II: must consume
    x = {v for v in x if succ[v]}
    if not x:
        return set()

    # phase III: no consumer outside the set
    x1 = set(x)
    while True:
        drop = {v for v in x if succ[v] - x}
        if not drop:
            break
        x -= drop
        if not x:
            break
    phase_iii = x1 - x

    # phase IV: non-terminal survivors
    phase_iv = _non_terminal(succ, set(range(len(succ))), x) if x else set()
    return phase_iii | phase_iv


def _weak_components(succ, pred, n):
    seen = [False] * n
    out = []
    for start in range(n):
        if seen[start]:
            continue
        stack, comp = [start], []
        seen[start] = True
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in succ[v] | pred[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
        out.append(comp)
    return out


def _coupling_merge(groups, Y):
    """Merge upstream sets whose representatives are coupled (Proposition S3-4).

    Two sets merge when the composition difference of their representatives lies in the span of the
    coupling matrix YD built from the sets themselves. Sharing a complex is not required, and this
    is the step that distinguishes the exhaustive path from the fast one: without it, sets that are
    coupled but disjoint stay apart, and a metabolite whose ACR is witnessed by one complex in each
    of them is missed.

    Span membership is decided through the left nullspace, as in :mod:`kimonu.acr`: with
    L = ker(YD'), a difference lies in the span exactly when the two representatives have the same
    signature L'Y[:, c]. So the merge is again a grouping of equal exact vectors, with no tolerance.
    """
    from kimonu.acr import composition_differences, _rows_of_left_nullspace

    sized = [g for g in groups if len(g) >= 2]
    if len(groups) < 2 or not sized:
        return groups

    YD = composition_differences(Y, sized)
    rows = _rows_of_left_nullspace(YD)
    if rows is None:
        return groups

    Ycsc = Y.tocsc()

    def signature(complex_index):
        start, end = Ycsc.indptr[complex_index], Ycsc.indptr[complex_index + 1]
        entries = list(zip(Ycsc.indices[start:end], Ycsc.data[start:end]))
        acc = {}
        for met, coeff in entries:
            for col, value in rows[int(met)]:
                acc[col] = acc.get(col, Fraction(0)) + Fraction(str(coeff)) * value
        return tuple(sorted((c, v) for c, v in acc.items() if v != 0))

    # a deterministic representative, rather than an arbitrary set element
    by_signature = {}
    for group in groups:
        by_signature.setdefault(signature(min(group)), []).append(group)

    merged = []
    for members in by_signature.values():
        union = set()
        for g in members:
            union |= g
        merged.append(union)
    return merged


def kinetic_modules(A, modules, balanced, Y=None):
    """Kinetic modules from the concordance partition.

    Args:
        A: complexes x reactions incidence matrix.
        modules: concordance modules as lists of complex indices (singletons included).
        balanced: indices of balanced complexes.
        Y: metabolites x complexes composition matrix. When given, coupled upstream sets are
            merged as well as overlapping ones (Proposition S3-4). Omitting it reproduces the
            reference implementation's overlap-only merging, which is coarser-grained nowhere and
            can therefore miss ACR.

    Returns a :class:`KineticResult`.
    """
    n = A.shape[0]
    succ, pred, _ = complex_graph(A)
    balanced = set(balanced)

    found = []
    for module in modules:
        module = set(module) - balanced
        if not module:
            continue
        survivors = autonomous(module | balanced, succ, pred)
        if survivors:
            found.append(survivors)

    # linkage classes made up entirely of balanced complexes are kinetic modules in their own right
    for comp in _weak_components(succ, pred, n):
        if comp and set(comp) <= balanced:
            found.append(set(comp))

    # merge modules that share a complex
    merged = []
    for group in found:
        overlapping = [m for m in merged if m & group]
        for m in overlapping:
            merged.remove(m)
        for m in overlapping:
            group = group | m
        merged.append(group)

    if Y is not None:
        merged = _coupling_merge(merged, Y)

    assigned = set().union(*merged) if merged else set()
    for i in range(n):
        if i not in assigned:
            merged.append({i})

    module_of = np.zeros(n, dtype=np.int64)
    for idx, group in enumerate(merged, start=1):
        for i in group:
            module_of[i] = idx

    return KineticResult(
        modules=merged,
        module_of=module_of,
        n_modules=sum(1 for m in merged if len(m) > 1),
        giant_module=max((len(m) for m in merged), default=0),
    )
