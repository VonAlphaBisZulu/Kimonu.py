"""Search for a natural replacement of the completeness toy network (now in Figure 2).

The toy must show both failure modes at once:
  P1  the iterated route misses a true ACR (or ACRR): a sign certificate proves it, the iterated route
      does not report it;
  P2  the tangent route reports an ACR (or ACRR) that is false: two positive steady states with the same
      rate constants differ in it (local ACR in the sense of Pascual-Escudero & Feliu, but not ACR).
Naturalness: few species, connected, one sink per species, at most bimolecular if possible.

The grammar is a supply-and-switch template. A protein X is made at a constant rate, degraded, and
assembled (2X or X) into an inactive kinase C; the active form P activates C by a positive feedback
(bimolecular C+P, cooperative C+2P, or through its dimer D); C, P, D may be degraded, P may revert to C.
Every network in the grammar is screened exactly (iterated vs tangent route), and the ones that differ
are classified per species and per pair: certify (proof) or search (counterexample), both two-point.

    python natural_example_search.py            # from scripts/kinetic_variants
"""
import itertools, json, os, sys, time
import numpy as np, scipy.sparse as sp

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src"), os.path.join(os.path.dirname(B), "tests")]
import records                                         # paths of the records (KIMONU_OUT)
import kimonu
from networks import matrices
from kinetic_variants.tangent_one_shot import one_shot_exact
from kinetic_variants.two_point import substrate_composition, certify, search, iterated_constraints

SPECIES = ["X", "C", "P", "D"]
REQUIRED = [("0", "X")]                                   # constant synthesis of X
POOL = [("X", "0"), ("X+X", "C"), ("X", "C"), ("C", "0"), ("P", "0"), ("D", "0"),
        ("C+P", "P+P"), ("C+P+P", "P+P+P"), ("P+P", "D"), ("D", "P+P"), ("C+D", "P+D"), ("P", "C")]


def used_species(R):
    return [s for s in SPECIES if any(s in c.split("+") for r in R for c in r if c != "0")]


def connected(R, S):
    """species graph: two species are linked when they meet in one reaction"""
    adj = {s: set() for s in S}
    for r in R:
        sp_ = {t for c in r if c != "0" for t in c.split("+")}
        for a in sp_:
            adj[a] |= sp_
    seen, todo = set(), [S[0]]
    while todo:
        s = todo.pop()
        if s not in seen:
            seen.add(s); todo += list(adj[s] - seen)
    return seen == set(S)


def routes(R, S):
    """iterated route (attached uptakes) and tangent route, both exact. None if the network has blocked
    reactions (no strictly positive steady-state flux)."""
    A, Y, N, cx = matrices(R, S)
    try:
        it = kimonu.iterated_route(kimonu.crn_from_matrices(A, Y), detach_uptakes=False)
    except ValueError:
        return None
    _, tacr, tgroups, dim = one_shot_exact(A, Y, N, detach=False, want_groups=True)
    pairs = lambda groups: {frozenset(p) for g in groups for p in itertools.combinations(sorted(g), 2)}
    return dict(A=A, Y=Y, N=N, it_acr=set(it.acr), it_acrr={frozenset(p) for p in it.acrr},
                tan_acr=set(tacr), tan_acrr=pairs(tgroups), tan_dim=dim)


def classify(R, S, rel, time_limit=60):
    """'true' (certificate infeasible, with the iterated equalities), 'false' (an explicit pair of steady
    states), or 'open'"""
    A, Y, N = rel["A"], rel["Y"], rel["N"]
    G = substrate_composition(A, Y); Nf = sp.csr_matrix(N, dtype=float)
    V, classes, _ = iterated_constraints(A, Y, N)
    out = {}
    m = len(S)
    cands = [("ACR", (i,)) for i in sorted(rel["tan_acr"])] + [("ACRR", tuple(sorted(p))) for p in rel["tan_acrr"]]
    for kind, idx in cands:
        c = np.zeros(m)
        if kind == "ACR": c[idx[0]] = 1
        else: c[idx[0]], c[idx[1]] = 1, -1
        name = S[idx[0]] if kind == "ACR" else f"{S[idx[0]]}/{S[idx[1]]}"
        if "infeasible" in certify(Nf, G, c, V=V, classes=classes, time_limit=time_limit)[0]:
            out[name] = "true"
        elif search(Nf, G, c, starts=8, seed=1, time_limit=time_limit) is not None:
            out[name] = "false"
        else:
            out[name] = "open"
    return out


def screen():
    hits = []
    t0 = time.time()
    for k in range(2, 9):
        for extra in itertools.combinations(POOL, k):
            R = REQUIRED + list(extra)
            S = used_species(R)
            if len(S) < 3 or not connected(R, S):
                continue
            rel = routes(R, S)
            if rel is None:
                continue
            missed = (rel["tan_acr"] - rel["it_acr"]) | {p for p in rel["tan_acrr"] - rel["it_acrr"]}
            if not missed:
                continue
            verdict = classify(R, S, rel, time_limit=20)
            truths = [n for n, v in verdict.items() if v == "true"]
            falses = [n for n, v in verdict.items() if v == "false"]
            it_names = {S[i] for i in rel["it_acr"]} | {f"{S[min(p)]}/{S[max(p)]}" for p in rel["it_acrr"]}
            p1 = [n for n in truths if n not in it_names]
            if p1 and falses:
                termolecular = any(len(c.split("+")) > 2 for r in R for c in r)
                sinks = sum(1 for r in R if r[1] == "0")
                hits.append(dict(reactions=R, species=S, missed_true=p1, false_local=falses,
                                 iterated=sorted(it_names), n_reactions=len(R), termolecular=termolecular,
                                 sinks=sinks, verdict=verdict))
    hits.sort(key=lambda h: (h["termolecular"], h["n_reactions"], len(h["species"])))
    return hits, time.time() - t0


if __name__ == "__main__":
    hits, sec = screen()
    out = records.write(records.EXPERIMENTS, "natural_example_search", "natural_example_hits.json")
    json.dump(dict(seconds=sec, hits=hits), open(out, "w"), indent=1)
    print(f"{len(hits)} networks with both failure modes ({sec:.0f} s); first ten:")
    for h in hits[:10]:
        print(h["n_reactions"], "termolecular" if h["termolecular"] else "bimolecular", h["reactions"],
              "| missed true:", h["missed_true"], "| false local:", h["false_local"])
