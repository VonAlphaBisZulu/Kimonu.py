"""Exact values of the iterated route's relations at given rate constants (kimonu.relation_values), 33 models.

Per model and uptake mode (attached, detached):
1. Round trip. x is a random small rational on every metabolite whose substrate coefficients are all
   integers and 1 on the rest, so ψ(x) stays rational; k = r* / ψ(x) (kimonu.rate_constants) with the
   integer flux r* of the tangent route. relation_values(k) must return every ACR value as x_i, every ACRR
   ratio as x_j / x_i and every module ratio as ψ_d(x) / ψ_c(x), exactly, and as many relations as the
   iterated route found in uptake_modes/<model>.json.
2. Binomial values. With k_j the j-th prime under a random permutation, a value that is a monomial in k
   factors completely over those primes; two permutations giving the same exponent per reaction confirm
   it. Every other value is a rational function of k that is no monomial (or a root of one).

    python relation_values.py <model>   -> $RV_OUT/<model>.json
"""
import json, os, random, sys, time
from fractions import Fraction
import numpy as np
from scipy import sparse
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.io import substrate_complexes
from kimonu.linalg import positive_flux
from run_kimonu_crn import triplets_to_sparse, to_exact

OUT = os.environ.get("RV_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_values")
M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
Ycsc = sparse.csc_matrix(Y)
ref = json.load(open(f"{B}/kinetic_variants/uptake_modes/{M}.json"))
if os.path.exists(f"{B}/kinetic_variants/uptake_modes/{M}.detached.json"):
    ref["detached"] = json.load(open(f"{B}/kinetic_variants/uptake_modes/{M}.detached.json"))["detached"]


def column(c):
    return {int(i): Fraction(str(v)) for i, v in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]],
                                                    Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]])}


def psi(c, x):
    out = Fraction(1)
    for i, y in column(c).items():
        if x[i] != 1:
            out *= x[i] ** int(y)
    return out


def primes(n):
    out, p = [], 2
    while len(out) < n:
        if all(p % d for d in out if d * d <= p):
            out.append(p)
        p += 1
    return out


def exponents(value, prime_of):
    """{reaction: exponent} when the Value factors over the assigned primes, else None"""
    out = {}
    for part, scale in ((value.r, 1), (value.b, Fraction(1, value.q))):
        for n, sign in ((part.numerator, 1), (part.denominator, -1)):
            for p, j in prime_of.items():
                while n % p == 0:
                    n //= p
                    out[j] = out.get(j, 0) + sign * scale
            if n != 1:
                return None
    return {j: e for j, e in out.items() if e}


def all_values(v):
    """(kind, key, Value) for every reported value"""
    yield from (("acr", i, x) for i, x in v.acr.items())
    yield from (("acrr", (r, j), x) for r, ratios in v.acrr for j, x in ratios.items())
    yield from (("module", (c, d), x) for c, ratios in v.modules for d, x in ratios.items())


t0 = time.time()
flux = positive_flux(N, 0, scale=10**6)
out = dict(model=M, n_reactions=nr, n_complexes=nc, n_metabolites=nm, seconds_flux=time.time() - t0)
fractional = {i for c in set(substrate_complexes(A).values()) for i, y in column(c).items() if y.denominator != 1}
rng = random.Random(0)
x = [Fraction(1) if i in fractional else Fraction(rng.randint(1, 5), rng.randint(1, 5)) for i in range(nm)]
P = primes(nr)
for mode, detach in (("attached", False), ("detached", True)):
    res = out[mode] = {}
    # 1) round trip
    k = kimonu.rate_constants(A, Y, flux, x)
    t0 = time.time()
    try:
        v = kimonu.relation_values(A, Y, N, k, detach_uptakes=detach)
    except (ValueError, NotImplementedError) as e:
        res["error"] = repr(e); print(M, mode, e, flush=True); continue
    res["seconds"], res["rounds"] = time.time() - t0, v.rounds
    wrong = [kind for kind, key, val in all_values(v)
             if val != (x[key] if kind == "acr" else x[key[1]] / x[key[0]] if kind == "acrr" else psi(key[1], x) / psi(key[0], x))]
    acrr_pairs = sum(len(r) * (len(r) + 1) // 2 for _, r in v.acrr)
    pairs = sum(len(r) * (len(r) + 1) // 2 for _, r in v.modules)
    it = ref[mode]["it"]
    res["round_trip"] = dict(acr=len(v.acr), acrr_pairs=acrr_pairs, module_pairs=pairs, wrong=len(wrong),
                             same_as_iterated=dict(acr=len(v.acr) == it["acr"], acrr_pairs=acrr_pairs == it["acrr_pairs"],
                                                   module_pairs=pairs == it["pairs"]))
    # 2) monomial values at prime rate constants, two permutations
    found = []
    for seed in (1, 2):
        perm = list(range(nr)); random.Random(seed).shuffle(perm)
        prime_of = {P[perm[j]]: j for j in range(nr)}
        try:
            vp = kimonu.relation_values(A, Y, N, [P[perm[j]] for j in range(nr)], detach_uptakes=detach)
        except (ValueError, NotImplementedError) as e:
            res["primes_error"] = repr(e); break
        found.append({(kind, key): exponents(val, prime_of) for kind, key, val in all_values(vp)})
    if len(found) == 2:
        tally = {}
        for key, e in found[0].items():
            kind = key[0]
            mono = e is not None and found[1].get(key) == e
            t = tally.setdefault(kind, dict(values=0, monomial=0, reactions_max=0))
            t["values"] += 1; t["monomial"] += mono
            if mono:
                t["reactions_max"] = max(t["reactions_max"], len(e))
        res["binomial"] = tally
    print(M, mode, json.dumps(res), flush=True)
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"), indent=1)
