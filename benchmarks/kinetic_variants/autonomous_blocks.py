"""Autonomous blocks: ACR from uniqueness, the step no linear or log-linear closure can take.

A species set S' is autonomous when every reaction that changes a species of S' (N_sj != 0) has all of its
substrate species in S'. Then dx_S'/dt depends on x_S' alone,

    f'(x_S') = N_S',R' (k_R' o x_S'^G_R',S'),    R' = reactions that change S',

and the S'-part of every positive steady state of the network is a positive steady state of f'. If N_S',R'
has rank |S'| (no conservation law inside S', so f' has a single compatibility class) and f' is injective
for every k > 0, f' has at most one positive zero for each k. All positive steady states of the network then
share x_S': every species of S' has ACR, for every choice of rate constants (vacuously if there is none).

Species already known to have ACR may stay outside S' even when S' depends on them: for fixed k their value
is the same at every positive steady state, so x_t^G_jt is a positive constant absorbed into k_j, and
injectivity for every k covers it. That is how proofs from this step and from the iterated route feed each
other.

Uptakes. Attached (zero complex, psi = 1): an uptake has no substrate species and runs at its rate constant;
it never breaks autonomy. Detached: every uptake has its own source of free level, a substrate species that
is never in S'; a species fed by an uptake is in no autonomous set, and the step proves nothing about it.
A conservation law inside S' leaves the totals free, so x_S' is not unique; the step does not apply then.

Injectivity is checked with the sign MILP of injectivity.py ("infeasible" proves it, up to the dynamic range
of its sign thresholds); blocks above ``max_block`` species are not tested.

    python autonomous_blocks.py <model> [--detach] [--max-block N] [--time-limit S]
"""
import json, os, sys, time
from collections import deque

import numpy as np
import scipy.sparse as sp

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]


def dependencies(Nf, G, detach_uptakes=False):
    """deps[s]: substrate species of the reactions that change s; open: species fed by a detached uptake."""
    Nc = sp.csc_matrix(Nf); Gr = sp.csr_matrix(G); m, q = Nc.shape
    subs = [set(Gr.indices[Gr.indptr[j]:Gr.indptr[j + 1]].tolist()) for j in range(q)]
    deps, changes, open_ = [set() for _ in range(m)], [set() for _ in range(m)], set()
    for j in range(q):
        rows = Nc.indices[Nc.indptr[j]:Nc.indptr[j + 1]][Nc.data[Nc.indptr[j]:Nc.indptr[j + 1]] != 0]
        for s in rows:
            deps[s] |= subs[j]; changes[s].add(j)
            if detach_uptakes and not subs[j]:
                open_.add(int(s))
    return deps, changes, open_


def closure(s, deps, known, open_, cap):
    """Smallest autonomous set containing s, with known-ACR species treated as constants; None if too large or open."""
    seen, todo = {s}, deque([s])
    while todo:
        t = todo.popleft()
        if t in open_:
            return None
        for u in deps[t]:
            if u not in seen and u not in known:
                seen.add(u); todo.append(u)
                if len(seen) > cap:
                    return None
    return frozenset(seen)


def block_system(Nf, G, block, changes):
    """N_S',R' and G_R',S' of an autonomous block."""
    S = sorted(block); R = sorted(set().union(*(changes[s] for s in S)))
    Nb = sp.csr_matrix(Nf)[S][:, R]; Gb = sp.csr_matrix(G)[R][:, S]
    return Nb, Gb


def autonomous_acr(Nf, G, known=(), detach_uptakes=False, max_block=60, time_limit=600, threads=8, log=print):
    """Species proven ACR by autonomous injective blocks, iterating as proofs remove dependencies."""
    from kinetic_variants.injectivity import certify_injective
    deps, changes, open_ = dependencies(Nf, G, detach_uptakes)
    known, proven, tested = set(known), set(), {}
    while True:
        new = set()
        blocks = {closure(s, deps, known | proven, open_, max_block) for s in range(Nf.shape[0]) if s not in known | proven}
        for b in sorted((b for b in blocks if b), key=len):
            if b in tested or b <= known | proven or b & new:
                continue
            Nb, Gb = block_system(Nf, G, b, changes)
            rank = np.linalg.matrix_rank(Nb.toarray()) if Nb.shape[1] else 0
            if rank < len(b):
                tested[b] = "conservation"; continue
            st, sec, wit = certify_injective(Nb, Gb, time_limit=time_limit, threads=threads)
            tested[b] = st
            log(f"block of {len(b)} species, {Nb.shape[1]} reactions: {st} ({sec:.1f}s)")
            if "infeasible" in st:
                new |= b - known - proven
        if not new:
            return proven, tested
        proven |= new


def augment(Y, modules, known):
    """Log closure (Proposition S3-4) and ACR (Remark S3-6) restarted with extra known-ACR species."""
    import kimonu
    from kimonu.acr import coupling_signatures
    sets, known = [set(m) for m in modules], set(known)
    while True:
        sig = coupling_signatures(Y, sets, known)
        merged = {}
        for s, key in zip(sets, sig):
            merged.setdefault(key, set()).update(s)
        merged = list(merged.values())
        found = known | set(kimonu.acr_acrr(Y, merged, known_acr=known).acr)
        if len(merged) == len(sets) and found == known:
            return merged, known
        sets, known = merged, found


def run(A, Y, N, Nf, detach_uptakes=False, max_block=60, time_limit=600, log=print):
    """Iterated route, then alternate autonomous blocks and the log closure until neither proves more."""
    import kimonu
    from kinetic_variants.two_point import substrate_composition
    G = substrate_composition(A, Y)
    it = kimonu.kinetic_modules(A, Y, N, detach_uptakes=detach_uptakes)
    base = set(kimonu.acr_acrr(Y, it.modules, known_acr=it.known_acr).acr)
    modules, known, rounds, tested = it.modules, set(base), 0, {}
    while True:
        rounds += 1
        proven, t = autonomous_acr(Nf, G, known, detach_uptakes, max_block, time_limit, log=log)
        tested.update(t)
        if not proven:
            break
        modules, known2 = augment(Y, modules, known | proven)
        log(f"round {rounds}: blocks prove {len(proven)}, closure brings ACR to {len(known2)} (iterated {len(base)})")
        if known2 == known:
            break
        known = known2
    acr = kimonu.acr_acrr(Y, modules, known_acr=known)
    return dict(iterated_acr=sorted(base), acr=sorted(acr.acr), acrr=acr.acrr, modules=modules,
                rounds=rounds,
                tested=sorted(((len(b), s) for b, s in tested.items()), key=lambda t: t[0]))


if __name__ == "__main__":
    import argparse
    from run_kimonu_crn import triplets_to_sparse, to_exact
    ap = argparse.ArgumentParser()
    ap.add_argument("model"); ap.add_argument("--detach", action="store_true")
    ap.add_argument("--max-block", type=int, default=60); ap.add_argument("--time-limit", type=float, default=600)
    a = ap.parse_args()
    crn = json.load(open(f"{B}/unblocked/{a.model}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    Nf = triplets_to_sparse(crn["S"], nm, nr); mets = crn["metabolites"]
    t = time.time()
    out = run(A, Y, to_exact(Nf), Nf, a.detach, a.max_block, a.time_limit, log=lambda s: print(f"[{a.model}] {s}", file=sys.stderr, flush=True))
    gained = sorted(set(out["acr"]) - set(out["iterated_acr"]))
    relf = f"/scratch/sadreyev_lab/pschneider/kimonu_twopoint/{a.model}_relations.json"
    rel_status = None
    if os.path.exists(relf):
        idx = {mm: i for i, mm in enumerate(mets)}
        acr, pairs = set(out["acr"]), {tuple(sorted(p)) for p in out["acrr"]}
        def holds(c):
            ks = [idx[k] for k in c]
            return all(k in acr for k in ks) if len(ks) == 1 else (tuple(sorted(ks)) in pairs or all(k in acr for k in ks))
        rel_status = [holds(r["c"]) for r in json.load(open(relf))["relations"]]
    sizes = {}
    for n, s in out["tested"]:
        sizes.setdefault(s if "infeasible" not in s else "injective", []).append(n)
    print(json.dumps(dict(model=a.model, detach=a.detach, seconds=round(time.time() - t, 1),
                          iterated_acr=len(out["iterated_acr"]), acr=len(out["acr"]), gained=[mets[i] for i in gained],
                          rounds=out["rounds"], blocks={k: [len(v), max(v)] for k, v in sizes.items()},
                          tangent_only_relations_proven=None if rel_status is None else [sum(rel_status), len(rel_status)])),
          flush=True)
