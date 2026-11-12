"""The reference-equivalent mode: the definitions of the reference, computed from nullspaces (Fig. 4E).

"The reference" is the implementation of Langary et al. (2025), in MATLAB and R, which the COCOA.jl authors call
the reference. It tests complex pairs with linear programs. This mode keeps its definitions and replaces the
pairwise tests by one nullspace, so that Kimonu can be compared with it like for like. It reproduces the results of
the exhaustive reference on all 11 networks on which that finished.

1. Concordance from M = A·K (:func:`kimonu.concordance`), without pinned rates, which mass action does not see.
2. Seeds. Each concordance module with two or more complexes, together with the balanced complexes. Like the
   reference, the mode does not seed singletons.
3. Upstream sets (:func:`upstream_set`). Within a seed, an upstream set is a set of complexes that no reaction
   enters from the rest of the network and that each of its complexes can leave, and its complexes are coupled
   (Langary et al., 2025). It is found on the complex graph (:func:`complex_graph`) in two phases. Entry complexes,
   produced by a reaction whose substrate lies outside the set, are removed until none is left. Then the terminal
   strongly connected components, which no reaction leaves, are removed.
4. Overlap merge (:func:`upstream_modules`). Upstream sets that share a complex merge, and so do the concordance
   modules behind them. The upstream sets of the merged modules are taken again, until nothing changes. The
   reference applies no other merge.
5. ACR and ACRR by the pairwise Shinar–Feinberg rule within each kinetic module (:func:`shinar_feinberg`), which
   derives them from coupled complexes that differ in one or two metabolites (Shinar & Feinberg, 2010).

The mode has no reference flux and no random point, and no uptake setting. The zero complex is a vertex of the
complex graph like any other. It takes about as long as the iterated route.
"""
from collections import namedtuple

from scipy import sparse

from kimonu.concordance import concordance
from kimonu.network import as_crn
from kimonu.readout import Readout, route_result

ComplexGraph = namedtuple("ComplexGraph", ["succ", "pred", "producers", "substrates"])
ComplexGraph.__doc__ = """The complex graph, one vertex per complex and one edge per reaction from substrate to product.

succ, pred   per complex, the set of complexes one reaction leads to, and comes from
producers    per complex, the reactions that produce it
substrates   per reaction, its substrate complexes (one for a mass-action network)
"""


def complex_graph(crn):
    """The :class:`ComplexGraph` of the network, read off its incidence matrix A."""
    A = sparse.csc_matrix(crn.A)
    n = A.shape[0]
    succ, pred = [set() for _ in range(n)], [set() for _ in range(n)]
    producers, substrates = [[] for _ in range(n)], []
    for j in range(A.shape[1]):
        rows, values = A.indices[A.indptr[j]:A.indptr[j + 1]], A.data[A.indptr[j]:A.indptr[j + 1]]
        subs = [int(r) for r, v in zip(rows, values) if v < 0]
        prods = [int(r) for r, v in zip(rows, values) if v > 0]
        substrates.append(subs)
        for p in prods:
            producers[p].append(j)
        for s in subs:
            for p in prods:
                succ[s].add(p)
                pred[p].add(s)
    return ComplexGraph(succ=succ, pred=pred, producers=producers, substrates=substrates)


def strong_components(succ, scope):
    """The strongly connected components of the complex graph restricted to ``scope`` (Tarjan's algorithm)."""
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


def linkage_classes(succ, pred, scope):
    """The linkage classes (weakly connected components) of the complex graph restricted to ``scope``."""
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


def upstream_set(candidates, graph):
    """The upstream set of a set of candidate complexes, a seed of balanced complexes and one concordance module.

    Phase 1 removes entry complexes, complexes that some reaction produces from a substrate outside the set, until
    none is left. A single pass is not enough, since each removal can make the next complex an entry complex.
    Phase 2 removes the terminal strongly connected components of what is left, all at once. Terminality is judged
    on the full complex graph. A component that feeds a complex outside the candidates still feeds, so it stays.
    Judging it within the candidates instead would delete it and strand its complexes as singletons.
    """
    current = set(candidates)
    # phase 1: drop entry complexes, to a fixed point
    while True:
        entry = {c for c in current
                 if any(any(sub not in current for sub in graph.substrates[j]) for j in graph.producers[c])}
        if not entry:
            break
        current -= entry
    # phase 2: drop the terminal strongly connected components, judged against the set after phase 1 and on the
    # full graph, all together. Removing them one at a time and testing again against the shrinking set cascades,
    # since each removal can strand the component that fed it.
    if current:
        for component in strong_components(graph.succ, current):
            if not any(graph.succ[v] - component for v in component):
                current -= component
    return current


def upstream_modules(crn, modules, balanced, max_rounds=20):
    """Kinetic modules from upstream sets merged by overlap, as the reference computes them, and the rounds taken.

    Args:
        crn: the network (:class:`kimonu.CRN`).
        modules: the concordance modules to seed with, as iterables of complex indices. :func:`reference_route`
            passes those with two or more complexes, as the reference does.
        balanced: the balanced complexes, added to every seed.
        max_rounds: cap on the rounds of upstream sets and merging, which normally end after a few.

    Each round takes the upstream set of every concordance module together with the balanced complexes, merges the
    upstream sets that share a complex, and merges the concordance modules behind them. The rounds end when the
    upstream sets and the concordance modules no longer change. Linkage classes made up entirely of balanced
    complexes that no upstream set holds form modules of their own, and every complex left over is a singleton.

    Returns (modules, rounds). The modules partition the complexes, singletons included.
    """
    n = crn.A.shape[0]
    graph = complex_graph(crn)
    balanced = set(balanced)
    concordance = [set(m) - balanced for m in modules]
    concordance = [m for m in concordance if m]
    previous, rounds, current = None, 0, []
    while rounds < max_rounds:
        rounds += 1
        # 1) the upstream set of every concordance module with the balanced complexes
        upstream_sets, origin = [], []
        for idx, module in enumerate(concordance):
            s = upstream_set(module | balanced, graph)
            if s:
                upstream_sets.append(s)
                origin.append(idx)
        if not upstream_sets:
            break
        # 2) upstream sets that share a complex merge, with the concordance modules behind them
        merged = []
        for group, sources in ((set(s), {o}) for s, o in zip(upstream_sets, origin)):
            for x in [x for x in merged if x[0] & group]:
                merged.remove(x)
                group, sources = group | x[0], sources | x[1]
            merged.append((group, sources))
        new_concordance = [set().union(*(concordance[i] for i in sources)) for _, sources in merged]
        new_concordance += [module for idx, module in enumerate(concordance)
                            if not any(idx in sources for _, sources in merged)]
        current = [s for s, _ in merged]
        # 3) stop when neither the upstream sets nor the concordance modules change
        state = (frozenset(frozenset(s) for s in current), frozenset(frozenset(m) for m in new_concordance))
        concordance = new_concordance
        if state == previous:
            break
        previous = state

    result = list(current)
    # linkage classes made up entirely of balanced complexes stand on their own
    assigned = set().union(*result) if result else set()
    for comp in linkage_classes(graph.succ, graph.pred, set(range(n))):
        if comp <= balanced and not (comp & assigned):
            result.append(comp)
            assigned |= comp
    result += [{c} for c in range(n) if c not in assigned]
    return result, rounds


def shinar_feinberg(crn, modules):
    """ACR metabolites and ACRR pairs by the pairwise Shinar–Feinberg rule within each kinetic module, as indices.

    Two complexes of one module whose compositions differ in exactly one metabolite give that metabolite ACR
    (Shinar & Feinberg, 2010). Two complexes that differ in exactly two metabolites and have exactly one further
    metabolite in common, the A + X against A + Y pattern, give the pair ACRR. As in the reference, coefficients and
    signs are not compared.

    Returns (acr, acrr), a set of metabolites and a set of sorted metabolite pairs.
    """
    acr, acrr = set(), set()
    for module in modules:
        members = sorted(module)
        if len(members) < 2:
            continue
        cols = {c: dict(crn.Yt[c]) for c in members}
        for k, c in enumerate(members):
            ca = cols[c]
            for d in members[k + 1:]:
                cb = cols[d]
                diff = [k for k in set(ca) | set(cb) if ca.get(k, 0) != cb.get(k, 0)]
                if len(diff) == 1:
                    acr.add(diff[0])
                elif len(diff) == 2 and len(set(ca) - set(diff)) == 1 and len(set(cb) - set(diff)) == 1:
                    acrr.add(tuple(sorted(diff)))
    return acr, acrr


def reference_route(network):
    """Kinetic modules, ACR and ACRR in the reference-equivalent mode, the definitions of Langary et al. (2025).

    Args:
        network: a :class:`kimonu.CRN` without blocked reactions, or a cobra model, which is preprocessed first
            (:func:`kimonu.as_crn`). Pass the CRN when calling several functions, so that preprocessing runs once.

    Returns a :class:`kimonu.RouteResult` with ``rounds`` the rounds of :func:`upstream_modules` and
    ``concordance`` the :class:`kimonu.ConcordanceResult` it started from. ``B`` and ``flux`` are None, since the
    mode reads no basis and needs no reference flux.
    """
    crn = as_crn(network)
    # 1) concordance from M = A·K, without the pins, since under mass action the bounds do not set the rates
    conc = concordance(crn._replace(pins={}))
    # 2) to 4) seeds, upstream sets and the overlap merge
    modules, rounds = upstream_modules(crn, [g for g in conc.modules if len(g) > 1], conc.balanced)
    # 5) ACR and ACRR by the Shinar–Feinberg rule. Each ACRR pair enters as a group of two, and route_result joins
    #    pairs that share a metabolite into the groups of acrr_groups.
    acr, acrr = shinar_feinberg(crn, modules)
    return route_result("reference", crn, Readout(modules, sorted(acr), [list(p) for p in sorted(acrr)]),
                        rounds=rounds, concordance=conc)
