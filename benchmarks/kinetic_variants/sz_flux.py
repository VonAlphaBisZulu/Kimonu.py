"""The uniformly drawn interior flux against the LP point: cost and agreement of both routes.

kimonu.interior_flux draws r* in two ways. draw="lp" is an LP point snapped exactly onto ker(N), its free coordinates
rounded at scale 10^6. An LP optimum is not a uniform random draw, so no probability bound applies to it.
draw="uniform", opt-in (the LP point is the default), is r* = K b with b = T c + s. c holds the free coordinates of the LP
point, s is uniform over the integers [0, 2^bits) in every free coordinate (bits = 64 by default), and T is large
enough that r* stays strictly positive whatever s is. A fixed shift keeps the degree of a polynomial, so a nonzero
polynomial of degree D in b vanishes at r* with probability at most D / 2^bits (Schwartz-Zippel). The script runs both routes at both draws and
reports wall time and whether kinetic modules, ACR and ACRR agree (detached uptakes).

Records written before 2026-11-11 call the LP arm "default", and their iterated route ran at the LP point without
rounding. figures/paper_numbers.py reads both.

    python sz_flux.py <model> [bits]   -> $SZ_OUT/<model>.json
"""
import json, os, sys, time
from fractions import Fraction
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from run_kimonu_crn import triplets_to_sparse, to_exact

OUT = os.environ.get("SZ_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_sz_flux")
M = sys.argv[1]; BITS = int(sys.argv[2]) if len(sys.argv) > 2 else 64
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = to_exact(triplets_to_sparse(crn["S"], nm, nr))


def run(draw):
    flux = kimonu.interior_flux(N, 0, draw, bits=BITS)
    t0 = time.time(); it = kimonu.iterated_route(A, Y, N, detach_uptakes=True, flux=flux); ti = time.time() - t0
    t0 = time.time(); tr = kimonu.tangent_route(A, Y, N, detach_uptakes=True, flux=flux); tt = time.time() - t0
    mods = lambda ms: sorted(sorted(int(c) for c in s) for s in ms if len(s) > 1)
    bits = max(int(x).bit_length() for x in flux)
    return dict(iterated=dict(seconds=ti, modules=mods(it.modules), acr=sorted(map(int, it.acr)),
                              acrr=sorted(map(list, it.acrr)), rounds=it.rounds, flux_bits=bits),
                tangent=dict(seconds=tt, modules=mods(tr.modules), acr=sorted(map(int, tr.acr)),
                             acrr=sorted(map(list, tr.acrr)), dimension=tr.dimension, flux_bits=bits))


out = dict(model=M, bits=BITS, lp=run("lp"), uniform=run("uniform"))
for route in ("iterated", "tangent"):
    a, b = out["lp"][route], out["uniform"][route]
    out[f"{route}_agree"] = all(a[k] == b[k] for k in ("modules", "acr", "acrr"))
    print(f"{M} {route}: lp {a['seconds']:.1f} s ({a['flux_bits']} bit), uniform {b['seconds']:.1f} s "
          f"({b['flux_bits']} bit), agree {out[f'{route}_agree']}", flush=True)
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"))
