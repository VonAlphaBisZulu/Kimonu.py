"""Kinetic modules from upstream sets with the extra merges of Kimonu's September route (superseded).

The reference-equivalent mode is kimonu.reference_route (src/kimonu/reference.py), the upstream sets of Langary et
al. (2025) merged by overlap only. This module adds the merges that Kimonu's route of September applied on top of
them when Y is given, and is kept for two one-off checks with no number in the paper,
matlab_reference/compare_kinetic_partitions.py and validate_acr_published.py. With Y=None, kinetic_modules is
kimonu.reference.upstream_modules.

With Y, three rules merge upstream sets. Sets that share a complex merge. Sets whose representatives are coupled
merge, which is decided by their coupling signatures, the rows of Yᵀ ker((YΔ)ᵀ) (coupling_signatures). If the
mass-action deficiency is one, every unbalanced concordance module merges. Merges feed back into the concordance
modules, so the whole thing is iterated to a fixed point, with the metabolites found to have ACR added to YΔ as they
are discovered.

On EnvZ-OmpR the coupling rule alone cannot merge the set containing XD with the set containing XD + Yp. Doing so
requires e_Yp to lie in the span of YΔ, which is the very question the merge would settle. The deficiency rule
breaks that circle. The iterated route of kimonu.iterated_route needs neither rule, since its feedback rows couple
these complexes in ker(W).
"""
from collections import namedtuple

import numpy as np

import kimonu
from kimonu.readout import complex_directions
from kimonu.reference import complex_graph, linkage_classes, strong_components, upstream_modules, upstream_set

KineticResult = namedtuple(
    "KineticResult",
    ["modules", "module_of", "n_modules", "giant_module", "known_acr", "iterations", "deficiency"],
)


def coupling_signatures(crn, modules, known_acr=()):
    """The row of Yᵀ ker((YΔ)ᵀ) of each module's first complex. Modules with equal signatures are coupled.

    YΔ holds the composition differences within the modules and the unit vectors of ``known_acr``. None when it has
    no column.
    """
    if not known_acr and not any(len(m) > 1 for m in modules):
        return None
    rows = complex_directions(crn, kimonu.free_directions(crn, modules, known_acr))
    return [rows[min(module)] for module in modules]


def is_weakly_reversible(succ, pred, n):
    """True when every linkage class is a single strongly connected component."""
    return all(len(strong_components(succ, component)) == 1
               for component in linkage_classes(succ, pred, set(range(n))))


def classical_deficiency(crn, succ, pred, nullity=None):
    """delta = n - l - rank(N): complexes, minus linkage classes, minus the rank of N = Y*A.

    The rank comes from the exact nullspace, not from an SVD. A float rank test would put a numerical tolerance in
    the middle of a discrete rule. The deficiency decides whether every unbalanced module merges, and the difference
    between delta = 1 and delta = 2 changes the answer. ``nullity`` may be supplied when the kernel is known from the
    concordance partition.
    """
    n, q = crn.A.shape
    linkage = len(linkage_classes(succ, pred, set(range(n))))
    if nullity is None:
        nullity = kimonu.concordance(crn).nullity
    rank = q - nullity
    return n - linkage - rank


def _result(result, known_acr, iterations, deficiency):
    n = sum(len(m) for m in result)
    module_of = np.zeros(n, dtype=np.int64)
    for idx, group in enumerate(result, start=1):
        for i in group:
            module_of[i] = idx
    return KineticResult(modules=result, module_of=module_of, n_modules=sum(1 for m in result if len(m) > 1),
                         giant_module=max((len(m) for m in result), default=0), known_acr=sorted(known_acr),
                         iterations=iterations, deficiency=deficiency)


def kinetic_modules(crn, modules, balanced, merges=True, max_iterations=20, nullity=None):
    """Kinetic modules from the concordance partition, with the September merges.

    Args:
        crn: the network (kimonu.CRN).
        modules: concordance modules as iterables of complex indices (singletons included).
        balanced: indices of balanced complexes.
        merges: the coupling and deficiency merges of September (default). Without them only overlapping upstream
            sets merge, as in kimonu.reference.upstream_modules.
        max_iterations: cap on the outer loop.
        nullity: dimension of ker(N), when known from the concordance step.

    Returns a :class:`KineticResult`.
    """
    if not merges:
        result, rounds = upstream_modules(crn, modules, balanced, max_iterations)
        return _result(result, (), rounds, -1)

    n = crn.A.shape[0]
    graph = complex_graph(crn)
    succ, pred = graph.succ, graph.pred
    balanced = set(balanced)
    concordance = [set(m) - balanced for m in modules]
    concordance = [m for m in concordance if m]
    initial_delta = classical_deficiency(crn, succ, pred, nullity=nullity)
    weakly_reversible = is_weakly_reversible(succ, pred, n)

    known_acr, n_merges, previous, iterations, current_modules = set(), 0, None, 0, []
    while iterations < max_iterations:
        iterations += 1
        upstream_sets, origin = [], []
        for idx, module in enumerate(concordance):
            s = upstream_set(module | balanced, graph)
            if s:
                upstream_sets.append(s)
                origin.append(idx)
        if not upstream_sets:
            break
        sets = _group_sets(upstream_sets, origin, crn, known_acr)
        # a merge of upstream sets is a merge of the concordance modules behind them
        new_concordance = [set().union(*(concordance[i] for i in sources)) for _, sources in sets]
        new_concordance += [module for idx, module in enumerate(concordance)
                            if not any(idx in sources for _, sources in sets)]
        n_merges += len(concordance) - len(new_concordance)
        current_modules = [s for s, _ in sets]
        known_acr |= set(kimonu.readout(crn, kimonu.free_directions(crn, current_modules, known_acr), False).acr)
        state = (frozenset(frozenset(s) for s in current_modules),
                 frozenset(frozenset(m) for m in new_concordance))
        concordance = new_concordance
        if state == previous:
            break
        previous = state

    # deficiency one: every unbalanced module belongs together (Langary et al., 2025)
    if initial_delta >= 0 and len(concordance) > 1:
        lower = 0 if weakly_reversible else 1
        upper = initial_delta - n_merges
        if lower == 1 and upper == 1:
            concordance = [set().union(*concordance)]
            upstream_sets = [s for s in (upstream_set(m | balanced, graph) for m in concordance) if s]
            sets = _group_sets(upstream_sets, list(range(len(upstream_sets))), crn, known_acr)
            current_modules = [s for s, _ in sets]

    result = list(current_modules)
    # linkage classes made up entirely of balanced complexes stand on their own
    assigned = set().union(*result) if result else set()
    for comp in linkage_classes(succ, pred, set(range(n))):
        if comp <= balanced and not (comp & assigned):
            result.append(comp)
            assigned |= comp
    result += [{i} for i in range(n) if i not in assigned]
    return _result(result, known_acr, iterations, initial_delta)


def _group_sets(upstream_sets, origin, crn, known_acr):
    """Merge upstream sets by overlap and by coupling signature. Returns (merged set, origin indices) pairs."""
    merged = []
    for group, sources in ((set(s), {o}) for s, o in zip(upstream_sets, origin)):
        for x in [x for x in merged if x[0] & group]:
            merged.remove(x)
            group, sources = group | x[0], sources | x[1]
        merged.append((group, sources))
    if len(merged) < 2:
        return merged
    signatures = coupling_signatures(crn, [g for g, _ in merged], known_acr)
    if signatures is None:
        return merged
    by_signature = {}
    for (group, sources), sig in zip(merged, signatures):
        key = by_signature.setdefault(sig, [set(), set()])
        key[0] |= group
        key[1] |= sources
    return [(g, s) for g, s in by_signature.values()]
