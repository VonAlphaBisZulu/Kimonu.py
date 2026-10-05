"""Complexes that are no reaction's substrate, joined to the kinetic module they are coupled to.

Langary et al. (2025, Supplement S.1) call two complexes coupled when psi_c = lambda psi_d at every positive
steady state; nothing restricts this to substrate complexes. Kimonu's iterated route leaves non-substrate
(terminal) complexes as singletons, since they appear in no rate law and so never in W. Their psi is still a
monomial of x, and the log closure decides their coupling like anyone else's: c joins the module of d when
y_c - y_d lies in im(YΔ), the span of the module's fixed composition differences plus the ACR unit vectors
(Proposition S3-4 with Remark S3-6).

Joining changes module membership only. The test is a span membership, so the column y_d - y_c it adds to
YΔ is already in im(YΔ): no ACR or ACRR changes, and no two existing modules can merge (if c matched d and e,
y_d - y_e would be in im(YΔ) and d, e would have merged already). Two terminal complexes can still join each
other. Feeding u_c = u_d back into W changes nothing either: W has a zero column for c, so the new row only
fixes u_c.

    python terminal_join.py <model> [--detach]     -> json line: modules, giant, ACR before and after
"""
import json, os, sys, time

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.acr import coupling_signatures
from kimonu.io import substrate_complexes, zero_complexes


def join_terminals(A, Y, result, detach_uptakes=False):
    """Modules of ``result`` (a KineticResult) with every non-substrate complex joined where it is coupled.

    With detached uptakes the zero complex stays out, as kinetic_modules keeps it out of every module."""
    substrate = set(substrate_complexes(A).values())
    skip = zero_complexes(Y) if detach_uptakes else set()
    multi = [set(m) for m in result.modules if len(m) > 1]
    lone = [next(iter(m)) for m in result.modules if len(m) == 1 and not (set(m) & skip)]
    sets = multi + [{c} for c in lone]
    sig = coupling_signatures(Y, sets, result.known_acr)
    if sig is None:
        return [set(m) for m in result.modules]
    merged = {}
    for s, key in zip(sets, sig):
        merged.setdefault(key, set()).update(s)
    out = list(merged.values()) + [{c} for c in skip]
    # only non-substrate complexes may move: every merge must involve at most one set with a substrate
    for g in out:
        parts = [s for s in sets if s <= g]
        assert sum(1 for s in parts if s & substrate) <= 1, "a join merged two substrate modules"
    return out


def summarise(modules):
    return dict(n_kinetic=sum(1 for m in modules if len(m) > 1), giant=max(len(m) for m in modules))


if __name__ == "__main__":
    import numpy as np
    from run_kimonu_crn import triplets_to_sparse, to_exact
    M = sys.argv[1]; detach = "--detach" in sys.argv
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))      # blocked reactions already removed
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    t = time.time()
    res = kimonu.kinetic_modules(A, Y, N, detach_uptakes=detach)
    before = kimonu.acr_acrr(Y, res.modules, known_acr=res.known_acr)
    joined = join_terminals(A, Y, res, detach_uptakes=detach)
    after = kimonu.acr_acrr(Y, joined, known_acr=res.known_acr)
    substrate = set(substrate_complexes(A).values())
    moved = sum(1 for m in joined if len(m) > 1 for c in m if c not in substrate)
    print(json.dumps(dict(model=M, detach=detach, seconds=round(time.time() - t, 1),
                          before=summarise([set(m) for m in res.modules]), after=summarise(joined),
                          terminal_joined=moved,
                          acr_same=sorted(before.acr) == sorted(after.acr),
                          acrr_same=sorted(map(tuple, before.acrr)) == sorted(map(tuple, after.acrr)))), flush=True)
