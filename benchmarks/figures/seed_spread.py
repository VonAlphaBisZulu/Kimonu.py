"""Seed spread for fig. S4: per model and arm, the relations that change with the seed.

For each stochastic arm (COCOA.jl both settings, three seeds; sampled reference, three seeds) and each model:
coupled complex pairs found by some seed but not by all (from the module partitions, by inclusion-exclusion over
partition meets), and ACR metabolites plus ACRR pairs among free metabolites found by some seed but not by all.
Kimonu: seeds 0, 1, 2 of r* (kinetic_variants/seed_invariance.py) gave identical modules, ACR and ACRR on all 33
models, so its spread is zero. Output: figures/seed_spread.json.

    python figures/seed_spread.py     # from benchmarks/
"""
import collections, glob, itertools, json, os, re, sys
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(B); sys.path.insert(0, "figures")
from gather_figdata import free_set, read_pairs
D = json.load(open("figdata.json"))
SEEDS = os.environ.get("SEED_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_seeds")

def P(l):
    return sum(v * (v - 1) // 2 for v in collections.Counter(x for x in l if x is not None).values())
def meet(*ls):
    return [tuple(x) if all(v is not None for v in x) else None for x in zip(*ls)]
def pair_spread(labs):
    """(pairs in some but not all, pairs in some) for module labellings over the same complexes"""
    inter = P(meet(*labs))
    union = sum(P(l) for l in labs) - sum(P(meet(a, b)) for a, b in itertools.combinations(labs, 2)) + inter
    return union - inter, union
def set_spread(sets):
    u, i = set.union(*sets), set.intersection(*sets); return len(u - i), len(u)
def rel(acr, acrr, free):
    return {("acr", a) for a in acr if a in free} | {("acrr",) + tuple(sorted(p)) for p in acrr if all(x in free for x in p)}

out = {"cocoa_false": {}, "cocoa_true": {}, "ref_sampled": {}, "kimonu_det": {}}
for arm, sub in (("cocoa_false", "cocoa_v110"), ("cocoa_true", "cocoa_eff")):
    for m in D["models"]:
        fs = sorted(glob.glob(f"sweep/{sub}/{sub}_{m}_HiGHS_p64_seed*.json"))
        if len(fs) < 2: continue
        js = [json.load(open(f)) for f in fs]; free = free_set(m)
        labs = [[x if x >= 0 else None for x in j["kinetic_module"]] for j in js]
        out[arm][m] = dict(seeds=len(js), pairs=pair_spread(labs), fixed=set_spread([rel(j["acr"], j["acrr"], free) for j in js]))
root = "reference_final_sampled/acr/Results"
for m in D["ref_sampled"]:
    labs, sets, free = [], [], free_set(m)
    for rec in D["ref_sampled"][m]:
        st = f"{m}_seed{rec['seed']}"
        mods = [x for x in json.load(open(f"{root}/modules_{st}.json"))["modules"] if x]
        l = [None] * (max(max(x) for x in mods) + 1)
        for i, x in enumerate(mods):
            for c in x: l[c] = i
        labs.append(l)
        import csv
        single = [r[1] for r in list(csv.reader(open(f"{root}/MetSingle/MetSingle_{st}.csv"), delimiter=" "))[1:]]
        sets.append(rel(single, read_pairs(f"{root}/MetDouble/MetDouble_{st}.csv"), free))
    L = max(len(l) for l in labs); labs = [l + [None] * (L - len(l)) for l in labs]
    out["ref_sampled"][m] = dict(seeds=len(labs), pairs=pair_spread(labs), fixed=set_spread(sets))
for f in sorted(glob.glob(f"{SEEDS}/*.json")):
    d = json.load(open(f)); assert all(d["identical"].values()), d["model"]
    out["kimonu_det"][d["model"]] = dict(seeds=len(d["runs"]), pairs=(0, None), fixed=(0, None))
json.dump(out, open("figures/seed_spread.json", "w"), indent=1)
for arm, v in out.items():
    nz_p = sum(x["pairs"][0] > 0 for x in v.values()); nz_f = sum(x["fixed"][0] > 0 for x in v.values())
    print(f"{arm:12s} {len(v):2d} models | pairs seed-dependent in {nz_p:2d} | ACR+ACRR seed-dependent in {nz_f:2d}")
