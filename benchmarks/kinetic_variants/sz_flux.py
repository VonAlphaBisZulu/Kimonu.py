"""A uniformly drawn interior flux, for a Schwartz-Zippel bound on genericity: cost and agreement.

Kimonu's r* is an LP point snapped exactly onto ker(N) (kimonu.linalg.positive_flux). An LP optimum is not a
uniform random draw, so no probability bound applies to it. Here r* = K b with b = T c + s: c is the
default draw (integer free coordinates at scale 10^6), s is uniform over the integers [0, 2^BITS) in every
free coordinate, and T is large enough that r* stays strictly positive whatever s is. A fixed shift keeps
the degree of a polynomial, so a nonzero polynomial of degree D in b vanishes at r* with probability at most
D / 2^BITS (Schwartz-Zippel). The script runs both routes with the default r* and with this one, and
reports wall time and whether kinetic modules, ACR and ACRR agree (detached uptakes).

    python sz_flux.py <model> [bits]   -> $SZ_OUT/<model>.json
"""
import json, os, sys, time
from fractions import Fraction
from math import ceil, lcm
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu, kimonu.kinetic, kimonu.tangent
from kimonu.linalg import null_space_rows, positive_flux
from run_kimonu_crn import triplets_to_sparse, to_exact

OUT = os.environ.get("SZ_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_sz_flux")
M = sys.argv[1]; BITS = int(sys.argv[2]) if len(sys.argv) > 2 else 64
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
rows = null_space_rows(N)


def uniform_flux(N, seed=0, scale=None):
    """r* = K (T c + s), s uniform over [0, 2^BITS) per free coordinate, integer and strictly positive."""
    c_flux = positive_flux(N, seed, scale=10**6)                    # integer, strictly positive, in ker N
    free = {}
    for i, row in enumerate(rows):
        if len(row) == 1:
            (j, v), = row.items()
            free.setdefault(j, c_flux[i] / v)                       # c: the free coordinates of the default r*
    rng = np.random.default_rng(10**6 + seed)
    s = {j: Fraction(int.from_bytes(rng.bytes((BITS + 7) // 8), "big") % 2**BITS) for j in free}
    # |K s|_i <= 2^BITS * sum_j |K_ij|, and K c >= min r_c > 0: T = 2 * 2^BITS * max_i sum_j |K_ij| / min r_c
    reach = max(sum(abs(v) for v in row.values()) for row in rows)
    T = ceil(2 * 2**BITS * reach / min(c_flux))
    flux = [sum((v * (T * free[j] + s[j]) for j, v in row.items()), Fraction(0)) for row in rows]
    assert all(x > 0 for x in flux)
    d = lcm(*(x.denominator for x in flux))
    return [x * d for x in flux]


def run(draw):
    kimonu.kinetic.positive_flux = draw if draw else positive_flux
    kimonu.tangent.positive_flux = draw if draw else positive_flux
    t0 = time.time(); km = kimonu.kinetic_modules(A, Y, N, detach_uptakes=True); ti = time.time() - t0
    R = kimonu.acr.readout(Y, kimonu.acr.free_directions(Y, km.modules), kimonu.io.zero_complexes(Y))
    t0 = time.time(); tr = kimonu.tangent_route(A, Y, N, detach_uptakes=True); tt = time.time() - t0
    mods = lambda ms: sorted(sorted(int(c) for c in s) for s in ms if len(s) > 1)
    return dict(iterated=dict(seconds=ti, modules=mods(km.modules), acr=sorted(map(int, R.acr)),
                              acrr=sorted(sorted(map(int, g)) for g in R.acrr), rounds=km.rounds,
                              flux_bits=max(int(x).bit_length() if Fraction(x).denominator == 1 else
                                            max(Fraction(x).numerator.bit_length(), Fraction(x).denominator.bit_length())
                                            for x in km.flux)),
                tangent=dict(seconds=tt, modules=mods(tr.modules), acr=sorted(map(int, tr.acr)),
                             acrr=sorted(map(list, tr.acrr)), dimension=tr.dimension,
                             flux_bits=max(int(x).bit_length() for x in tr.flux)))


out = dict(model=M, bits=BITS, default=run(None), uniform=run(uniform_flux))
for route in ("iterated", "tangent"):
    a, b = out["default"][route], out["uniform"][route]
    out[f"{route}_agree"] = all(a[k] == b[k] for k in ("modules", "acr", "acrr"))
    print(f"{M} {route}: default {a['seconds']:.1f} s ({a['flux_bits']} bit), uniform {b['seconds']:.1f} s "
          f"({b['flux_bits']} bit), agree {out[f'{route}_agree']}", flush=True)
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"))
