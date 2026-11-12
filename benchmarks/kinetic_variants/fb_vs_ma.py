"""Flux balance against mass action: what the rate law adds to flux coupling and to concordance.

Flux balance looks at all r >= 0 with N r = 0. Under mass action with fixed rate constants only some of
these flux vectors are steady states, so every relation that holds over the flux cone still holds, and more
can hold. This script counts both sides on one network, for reactions and for complexes.

Reactions. Under flux balance two reactions are fully coupled when their rows of K = ker(N) are proportional.
Under mass action r_j = r*_j u_s(j), with u the monomial of the substrate complex s(j) in units of the
reference steady state, so two rates keep a fixed ratio exactly when u_s(i) / u_s(j) is fixed: the two
reactions leave the same complex, or two complexes of one kinetic module. A detached uptake leaves a source
of its own, which is no complex and belongs to no module, so the classes are read over all columns of W:

    iterated route   columns with proportional rows in ker(W), W carrying the feedback rows u_c = u_d of
                     the final modules. Every coupling holds at every positive steady state (lower bound).
    tangent route    columns with equal rows of [Y^T T ; T_sources], T = ker(W Y^T) with one extra
                     coordinate per source. No coupling is missed (upper bound).

A reaction has a fixed rate when the monomial of its substrate complex is fixed (zero row of Y^T B).
Blocked, directional and partial coupling have no counterpart, since every rate is positive at a positive
steady state.

Complexes. Under flux balance two complexes are concordant when their rows of M = A K are proportional, and
balanced when the row is zero. Under mass action they are coupled when they share a kinetic module. The first
compares activities, the second complex concentrations, so neither contains the other, and the script counts the pairs in
each and in both.

Checks of this run: flux coupling refines the iterated classes, which refine the tangent classes; the
classes of ker(W) restricted to the complexes are the kinetic modules; and the modules of both routes
reproduce the record of uptake_modes.py; and kimonu.kinetic_flux_coupling, the same lift in the package,
reproduces the classes of both routes. Its fixed rates also count detached uptakes whose source is coupled
to a fixed complex (field fixed_rate_src). The tangent route runs on every model with detached uptakes (the
sources need it) and, with attached uptakes, only where that record reports a gap.

    python kinetic_variants/fb_vs_ma.py <model> [attached]   -> kinetic_variants/fb_vs_ma/<model>.json
"""
import json, os, sys, time
from collections import Counter
from fractions import Fraction
import numpy as np
from scipy import sparse
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.readout import complex_directions, free_directions, readout
from kimonu.io import substrate_complexes, zero_complexes
from kimonu.linalg import exact, null_space_rows, proportional_classes
from kimonu.steady_state import interior_flux, steady_state_map
from run_kimonu_crn import triplets_to_sparse, to_exact


def n_pairs(sizes):
    return sum(s * (s - 1) // 2 for s in sizes)


def labels(groups, n):
    """One label per item from a list of groups; items in no group get a label of their own."""
    out = [-1 - i for i in range(n)]
    for k, g in enumerate(groups):
        for i in g:
            out[i] = k
    return out


def join(a, b):
    """The finest partition that both labelings refine (union-find over their classes)."""
    parent = list(range(len(a)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for lab in (a, b):
        first = {}
        for i, l in enumerate(lab):
            parent[find(i)] = find(first.setdefault(l, i))
    return [find(i) for i in range(len(a))]


def refines(a, b):
    """True when every class of labeling a lies inside one class of labeling b."""
    return len(set(zip(a, b))) == len(set(a))


def summary(lab):
    sizes = sorted((s for s in Counter(lab).values() if s > 1), reverse=True)
    return dict(pairs=n_pairs(sizes), groups=len(sizes), members=sum(sizes), largest=max(sizes, default=1))


def both(a, b):
    """Pairs that share a class in both labelings."""
    return n_pairs(Counter(zip(a, b)).values())


def module_record(modules):
    sizes = [len(g) for g in modules if len(g) > 1]
    return dict(modules=len(sizes), giant=max(sizes, default=1), pairs=n_pairs(sizes))


def iterated(A, Y, N, detach, column):
    """Kinetic modules by the iterated route, the class of every column of W, and the fixed complexes."""
    n = A.shape[0]
    it = kimonu.iterated_route(A, Y, N, detach_uptakes=detach)
    W, (m, n_columns) = steady_state_map(A, Y, N, it.flux, detach)
    equal = [(min(g), c) for g in it.modules for c in g if c != min(g)]
    for k, (c, d) in enumerate(equal):
        W[(m + k, c)], W[(m + k, d)] = Fraction(1), Fraction(-1)
    zero, classes = proportional_classes(null_space_rows(exact(((i, j, v) for (i, j), v in W.items()),
                                                               (m + len(equal), n_columns))))
    assert not zero and n_columns == max(n, 1 + max(column.values())), "columns"
    of_column = labels(classes, n_columns)
    # the classes of ker(W), restricted to the complexes, are the kinetic modules
    skip = zero_complexes(Y) if detach else set()
    assert {frozenset(g) for g in it.modules if len(g) > 1} == \
        {frozenset(g) for g in ({c for c in cl if c < n and c not in skip} for cl in classes) if len(g) > 1}, "modules"
    YtL = complex_directions(Y, free_directions(Y, it.modules, it.acr))
    fixed = {c for c in range(n) if c not in skip and not YtL[c]}
    return it.modules, of_column, fixed, it.rounds


def tangent(A, Y, N, detach):
    """Kinetic modules by the tangent route (as kimonu.tangent_route), the signature of every column of W,
    and the fixed complexes."""
    n, m = A.shape[0], Y.shape[0]
    skip = zero_complexes(Y) if detach else set()
    flux = interior_flux(N, 0)
    W, (_, n_columns) = steady_state_map(A, Y, N, flux, detach)
    Ycsc = sparse.csc_matrix(Y)
    J = {}
    for (i, c), w in W.items():
        if c >= n:
            J[(i, m + c - n)] = w
            continue
        for k, y in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]]):
            J[(i, int(k))] = J.get((i, int(k)), Fraction(0)) + w * Fraction(str(y))
    rows = null_space_rows(exact(((i, k, v) for (i, k), v in J.items()), (m, m + n_columns - n)))
    rows = [tuple(sorted(r.items())) for r in rows]
    modules = readout(Y, rows[:m], skip).modules + [{c} for c in sorted(skip)]
    signature = complex_directions(Y, rows[:m]) + rows[m:]
    keys = {}
    of_column = [keys.setdefault(s, len(keys)) for s in signature]
    for c in skip:                                 # a detached zero complex is no reaction's substrate
        of_column[c] = -1 - c
    fixed = {c for c in range(n) if c not in skip and not signature[c]}
    return modules, of_column, fixed


def reactions(fb, same, ma, fixed, column, n):
    """Reaction pairs with a fixed rate ratio under mass action, against flux coupling."""
    q = len(fb); lab = [ma[column[j]] for j in range(q)]
    assert refines(fb, lab), "flux coupling must refine the mass-action classes"
    free = join(fb, same)                          # what follows from flux coupling and shared substrates alone
    assert refines(free, lab)
    out = summary(lab)
    out.update(beyond_flux_coupling=out["pairs"] - n_pairs(Counter(fb).values()),
               beyond_flux_coupling_and_shared_substrate=out["pairs"] - n_pairs(Counter(free).values()),
               fixed_rate=sum(column[j] < n and column[j] in fixed for j in range(q)))
    return out, lab


def check_src(A, Y, N, detach, route, lab):
    """kimonu.kinetic_flux_coupling must give the classes of this script; returns its fixed-rate count."""
    src = kimonu.kinetic_flux_coupling(A, Y, N, detach, route=route)
    assert len(set(zip(src.class_of, lab))) == len(set(src.class_of)) == len(set(lab)), (route, "src classes")
    return len(src.fixed_rate)


def complexes(conc, balanced_too, modules, n):
    """Concordant pairs against kinetically coupled pairs. The two relations are not nested."""
    km = labels([g for g in modules if len(g) > 1], n)
    in_module = {c for g in modules if len(g) > 1 for c in g}
    in_conc = {c for g in conc.modules if len(g) > 1 for c in g}
    out = {}
    for name, fb in (("concordant", labels(conc.modules, n)), ("concordant_or_balanced", balanced_too)):
        s = summary(fb); b = both(fb, km)
        out[name] = dict(flux_balance=s["pairs"], mass_action=n_pairs(Counter(km).values()), both=b,
                         only_flux_balance=s["pairs"] - b, only_mass_action=n_pairs(Counter(km).values()) - b,
                         flux_balance_refines_mass_action=refines(fb, km))
    out.update(members_concordance=len(in_conc), members_modules=len(in_module),
               members_both=len(in_conc & in_module), balanced_in_modules=len(set(conc.balanced) & in_module))
    return out


if __name__ == "__main__":
    M = sys.argv[1]; name = "attached" if "attached" in sys.argv[2:] else "detached"; detach = name == "detached"
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    f = f"{B}/kinetic_variants/uptake_modes/{M}.json"; one = f.replace(".json", f".{name}.json")
    stored = json.load(open(one if os.path.exists(one) else f))[name]
    substrate, zero = substrate_complexes(A), zero_complexes(Y)
    column, sources = {}, {}                      # the column of W that carries each reaction, as in steady_state_map
    for j, c in substrate.items():
        column[j] = sources.setdefault(j, nc + len(sources)) if detach and c in zero else c
    assert len(column) == nr
    out = dict(model=M, uptakes=name, n_reactions=nr, n_complexes=nc, n_sources=len(sources), gap=stored["gap"])

    # flux balance: full coupling from K, concordance from M = A K
    t0 = time.time()
    blocked, classes = proportional_classes(null_space_rows(N))
    assert not blocked, "blocked reactions"
    fb = labels(classes, nr); same = [column[j] for j in range(nr)]
    conc = kimonu.concordance_from_matrices(A, N)
    balanced_too = labels(conc.modules + [conc.balanced], nc)
    out["flux_coupling"] = summary(fb)
    out["shared_substrate"] = dict(summary(same), also_flux_coupled=both(fb, same))
    out["flux_coupling_and_shared_substrate"] = summary(join(fb, same))
    out["concordance"] = dict(summary(labels(conc.modules, nc)), balanced=len(conc.balanced))
    out["seconds_flux_balance"] = time.time() - t0

    # mass action, iterated route
    t0 = time.time()
    modules, of_column, fixed, out["rounds"] = iterated(A, Y, N, detach, column)
    assert module_record(modules) == {k: stored["it"][k] for k in ("modules", "giant", "pairs")}, (M, "iterated")
    out["rates_iterated"], lab_it = reactions(fb, same, of_column, fixed, column, nc)
    out["rates_iterated"]["fixed_rate_src"] = check_src(A, Y, N, detach, "iterated", lab_it)
    out["complexes_iterated"] = complexes(conc, balanced_too, modules, nc)
    out["fixed_complexes_iterated"] = len(fixed)
    out["seconds_iterated"] = time.time() - t0

    # mass action, tangent route
    out["rates_tangent"] = out["complexes_tangent"] = None      # None: the same lists as the iterated route
    if detach or stored["gap"]:
        t0 = time.time()
        modules_t, of_column, fixed_t = tangent(A, Y, N, detach)
        assert module_record(modules_t) == {k: stored["tan"][k] for k in ("modules", "giant", "pairs")}, (M, "tangent")
        out["rates_tangent"], lab_tan = reactions(fb, same, of_column, fixed_t, column, nc)
        out["rates_tangent"]["fixed_rate_src"] = check_src(A, Y, N, detach, "tangent", lab_tan)
        assert refines(lab_it, lab_tan), "the iterated classes must refine the tangent classes"
        out["complexes_tangent"] = complexes(conc, balanced_too, modules_t, nc)
        out["fixed_complexes_tangent"] = len(fixed_t)
        out["rates_routes_agree"] = len(set(lab_it)) == len(set(lab_tan)) and \
            out["rates_iterated"]["fixed_rate"] == out["rates_tangent"]["fixed_rate"]
        out["seconds_tangent"] = time.time() - t0
    else:
        assert all(stored["it"][k] == stored["tan"][k] for k in ("modules", "giant", "pairs")), M
        out["rates_routes_agree"] = True          # no sources and equal modules: equal rate classes
    os.makedirs(f"{B}/kinetic_variants/fb_vs_ma", exist_ok=True)
    json.dump(out, open(f"{B}/kinetic_variants/fb_vs_ma/{M}{'' if detach else '.attached'}.json", "w"), indent=1)
    print(M, name, json.dumps(out), flush=True)
