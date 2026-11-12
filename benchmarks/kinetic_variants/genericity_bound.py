"""The Schwartz–Zippel bound of the uniformly drawn reference flux, per network, from exact ranks.

kimonu.interior_flux(net, seed, draw="uniform", bits) draws r* = K (T c + s), with s uniform over the integers
[0, 2^bits) in every free coordinate (the LP point, draw="lp", is the default and has no such bound). Both routes read
their relations off matrices whose entries are linear in r*: W(r*), with the feedback rows of each round, for the
iterated route, and J(r*) = W(r*) Yᵀ for the tangent route. The result at r* differs from the generic one only if
r* is a zero of a nonzero polynomial among the minors that decide the rank of these matrices and whether a candidate
vector (e_i for ACR, e_i - e_j for ACRR, y_c - y_d for a complex pair) lies in the row space. Each has degree at most
rank + 1, and a shift by the fixed T c keeps the degree, so each vanishes at r* with probability at most
(rank + 1) / 2^bits. With the accounting of cert_modp.py (one test per candidate relation and one for the rank),

    tests          = m + m (m - 1) / 2 + n' (n' - 1) / 2 + 1     (n' complexes without the zero complex)
    tangent bound  = tests (rank J + 1) / 2^bits
    iterated bound = rounds tests (max rank W + 1) / 2^bits

bounds the probability that the reported result is not the generic one. The ranks are exact and taken at the
uniform r* itself (seed 0, detached uptakes). rank W comes from the rounds of the iterated route, repeated here
(seconds per network), and checked against the record of sz_flux.py at the same point. rank J is the number of
coordinates of J minus dim ker J, the tangent dimension of that record, so the slow tangent route is not run again.
bits comes from the sz_flux.py record, or the default of kimonu.interior_flux. A later change of bits needs a new
sz_flux.py record only for the check and the tangent dimension, which does not depend on the point.

    python kinetic_variants/genericity_bound.py [model ...]   -> $GB_OUT/<model>.json, $GB_OUT/summary.json

SZ_DIR is the sz_flux.py output read (default the rerun of 2026-11-11), GB_OUT the output directory.
"""
import inspect, json, os, statistics, subprocess, sys, time
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu import free_directions, readout, steady_state_map
from kimonu.network import zero_complexes
from kimonu.linalg import nullspace
from kimonu.readout import row_readout
from run_kimonu_crn import triplets_to_sparse

RERUN = "/scratch/sadreyev_lab/pschneider/kimonu_rerun_1111"
SZ = os.environ.get("SZ_DIR", f"{RERUN}/sz_flux_33")
OUT = os.environ.get("GB_OUT", f"{RERUN}/genericity_bound")
SEED = 0
BITS_DEFAULT = inspect.signature(kimonu.interior_flux).parameters["bits"].default


def iterated_ranks(net, flux):
    """The rounds of kimonu.iterated_route (detached uptakes) with the exact rank of W in each round."""
    n = net.A.shape[0]
    equal, ranks = [], []
    while True:
        W, n_columns = steady_state_map(net, flux, True)
        m = len(W)
        W += [tuple(sorted(((c, 1), (d, -1)))) for c, d in equal]          # the feedback rows u_c - u_d = 0
        rows = nullspace(W, n_columns)
        ranks.append(n_columns - len({j for r in rows for j, _ in r}))
        _, classes, _ = row_readout(rows)
        coupled = [s for s in ({c for c in g if c < n} for g in classes) if s]
        found = readout(net, free_directions(net, coupled))
        if {frozenset(s) for s in found.modules} == {frozenset(s) for s in coupled}:      # as kimonu.iterated_route
            break
        equal = [(min(s), c) for s in found.modules for c in s if c != min(s)]
    pairs = sorted(sorted((int(a), int(b))) for g in found.acrr_groups for i, a in enumerate(sorted(g)) for b in sorted(g)[i + 1:])
    return dict(ranks=ranks, rounds=len(ranks), n_columns=n_columns, rows=m,
                modules=sorted(sorted(int(c) for c in s) for s in found.modules if len(s) > 1),
                acr=sorted(map(int, found.acr)), acrr=pairs)


def network(M):
    crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    net = kimonu.crn_from_matrices(A, Y)                 # N = Y·A, the S of the crn.json
    rec = json.load(open(f"{SZ}/{M}.json"))
    bits = rec.get("bits", BITS_DEFAULT)
    t0 = time.time()
    flux = kimonu.interior_flux(net, SEED, "uniform", bits=bits)
    it = iterated_ranks(net, flux)
    seconds = time.time() - t0
    u = rec["uniform"]
    check = dict(modules=it["modules"] == u["iterated"]["modules"], acr=it["acr"] == u["iterated"]["acr"],
                 acrr=it["acrr"] == sorted(sorted(p) for p in u["iterated"]["acrr"]), rounds=it["rounds"] == u["iterated"]["rounds"])
    assert all(check.values()), (M, check)
    skip = zero_complexes(net)
    ncx = nc - len(skip)
    tests = nm + nm * (nm - 1) // 2 + ncx * (ncx - 1) // 2 + 1
    coords_J = nm + it["n_columns"] - nc                  # metabolites, then the detached sources
    dim = u["tangent"]["dimension"]
    rank_J = coords_J - dim
    return dict(model=M, bits=bits, seed=SEED, draw="uniform", detached_uptakes=True, m=nm, n=nc, zero_complexes=len(skip),
                tests=tests, W_columns=it["n_columns"], W_rows=it["rows"], ranks_W=it["ranks"], rounds=it["rounds"],
                J_coordinates=coords_J, tangent_dimension=dim, tangent_dimension_lp=rec["lp"]["tangent"]["dimension"],
                rank_J=rank_J, bound_tangent=tests * (rank_J + 1) / 2 ** bits,
                bound_iterated=it["rounds"] * tests * (max(it["ranks"]) + 1) / 2 ** bits,
                check_against_sz_flux=check, seconds_iterated_ranks=seconds, sz_record=f"{SZ}/{M}.json")


if __name__ == "__main__":
    models = sys.argv[1:] or [l.strip() for l in open(f"{B}/kinetic_variants/models33.txt") if l.strip()]
    os.makedirs(OUT, exist_ok=True)
    for M in models:
        d = network(M)
        json.dump(d, open(f"{OUT}/{M}.json", "w"), indent=1)
        print(f"{M}: tests {d['tests']}, rank J {d['rank_J']}, ranks W {d['ranks_W']}, bound tangent {d['bound_tangent']:.2e}, "
              f"iterated {d['bound_iterated']:.2e} ({d['seconds_iterated_ranks']:.1f} s)", flush=True)
    rows = [json.load(open(f"{OUT}/{l.strip()}.json")) for l in open(f"{B}/kinetic_variants/models33.txt")
            if l.strip() and os.path.exists(f"{OUT}/{l.strip()}.json")]
    head = subprocess.run(["git", "-C", B, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    summary = dict(networks=len(rows), bits=sorted({r["bits"] for r in rows}), kimonu=head)
    for key in ("bound_tangent", "bound_iterated"):
        v = [r[key] for r in rows]
        summary[key] = dict(min=min(v), median=statistics.median(v), max=max(v),
                            argmax=max(rows, key=lambda r: r[key])["model"]) if v else None
    summary["bound_either"] = dict(max=max(max(r["bound_tangent"], r["bound_iterated"]) for r in rows)) if rows else None
    json.dump(summary, open(f"{OUT}/summary.json", "w"), indent=1)
    print(json.dumps(summary))
