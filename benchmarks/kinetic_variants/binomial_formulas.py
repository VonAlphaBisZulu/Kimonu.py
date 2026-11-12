"""Monomial formulas from flux coupling (kimonu.binomial_formulas) on the 33 models, checked against exact values.

Per model and uptake mode: k = r* / ψ(x) as in relation_values.py, the exact values there
(kimonu.relation_values), and every monomial formula evaluated at the same k, which must give the same
Value. Reported: how many ACR values, ACRR ratios and module ratios have a monomial formula, the largest
number of rate constants in one, and the constraints on k.

    python binomial_formulas.py <model>   -> $BF_OUT/<model>.json
"""
import json, os, random, sys, time
from fractions import Fraction
import numpy as np
from scipy import sparse
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.network import substrate_complexes
from kimonu import interior_flux
from run_kimonu_crn import triplets_to_sparse

OUT = os.environ.get("BF_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_binomial")
M = sys.argv[1]
crn = json.load(open(f"{B}/unblocked/{M}_fixed.crn.json"))
nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
net = kimonu.crn_from_matrices(A, Y)                 # N = Y·A, the S of the crn.json
Ycsc = sparse.csc_matrix(Y)
fractional = {int(i) for c in set(substrate_complexes(net).values())
              for i, y in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]], Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]])
              if Fraction(str(y)).denominator != 1}
rng = random.Random(0)
x = [Fraction(1) if i in fractional else Fraction(rng.randint(1, 5), rng.randint(1, 5)) for i in range(nm)]
k = kimonu.rate_constants(net, interior_flux(net, 0, draw="lp"), x)
out = dict(model=M)
for mode, detach in (("attached", False), ("detached", True)):
    v = kimonu.relation_values(net, k, detach_uptakes=detach)
    t0 = time.time()
    b = kimonu.binomial_formulas(net, detach_uptakes=detach, draw="lp")
    res = out[mode] = dict(seconds=time.time() - t0, constraints=len(b.constraints),
                           constraints_hold=all(c(k) == 1 for c in b.constraints))
    triples = [("acr", v.acr[i], b.acr[i]) for i in v.acr]
    for kind, values, formulas in (("acrr", v.acrr, b.acrr), ("module", v.modules, b.modules)):
        formulas = dict(formulas)
        triples += [(kind, val, formulas[r][j]) for r, ratios in values for j, val in ratios.items()]
    for kind, val, f in triples:
        t = res.setdefault(kind, dict(values=0, monomial=0, wrong=0, rate_constants_max=0))
        t["values"] += 1
        if f is not None:
            t["monomial"] += 1; t["wrong"] += f(k) != val
            t["rate_constants_max"] = max(t["rate_constants_max"], len(f.exponents))
    print(M, mode, json.dumps(res), flush=True)
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/{M}.json", "w"), indent=1)
