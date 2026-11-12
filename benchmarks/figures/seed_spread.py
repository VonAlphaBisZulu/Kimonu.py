"""Seed spread for fig. S4: per model and arm, the relations that change with the seed.

For each stochastic arm (COCOA.jl both settings, three seeds; sampled reference, three seeds) and each model:
coupled complex pairs found by some seed but not by all (from the module partitions, by inclusion-exclusion over
partition meets), and ACR metabolites plus ACRR pairs among free metabolites found by some seed but not by all
("fixed", as listed). "fixed_ratios" counts the same with ACRR as independent ratios (k - 1 per group of k), so it
does not depend on how a tool lists its pairs. Each seed's ACR and ACRR (all metabolites) are closed under what they
imply, as in pair_sets.py, and restricted to the free metabolites, which gives a partition of the free metabolites and
a block of fixed ones. The independent relations of a partition are its rank, and fixed_ratios is the rank of what
some seed finds (the join of the partitions) minus the rank of what every seed finds (their meet).
Kimonu: the three runs of the timing arm det (seeds 0, 1, 2 of r*, KIMONU_DET_DIR) must give identical modules, ACR
and ACRR on every model, which is checked here, so its spread is zero. Output: figures/seed_spread.json.

    python figures/seed_spread.py     # from benchmarks/
"""
import collections, glob, itertools, json, os, re, sys
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(B); sys.path.insert(0, "figures")
from gather_figdata import free_set, read_pairs
D = json.load(open("figdata.json"))
DET = os.environ.get("KIMONU_DET_DIR", "/scratch/sadreyev_lab/pschneider/kimonu_timing_seeds/det")
FIXED = ("fixed",)

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
def closed(runs, free):
    """block of each free metabolite (and of FIXED) under the ACR and ACRR of all runs together, closed under implication"""
    parent = {FIXED: FIXED}
    def find(a):
        while parent.setdefault(a, a) != a:
            a = parent[a]
        return a
    for acr, acrr in runs:
        for a in acr: parent[find(a)] = find(FIXED)
        for p in acrr:
            p = list(p)
            if len(set(p)) == 2: parent[find(p[0])] = find(p[1])
    return {a: find(a) for a in list(parent) if a == FIXED or a in free}
def rank(labels):
    return len(labels) - len(set(labels.values()))
def ratio_spread(runs, free):
    """(rank of the join minus rank of the meet, rank of the join) of the runs' partitions of the free metabolites"""
    parts = [closed([r], free) for r in runs]
    join = closed(runs, free)
    meet = {a: tuple(p.get(a, ("alone", a)) for p in parts) for a in join}
    return rank(join) - rank(meet), rank(join)

out = {"cocoa_false": {}, "cocoa_true": {}, "ref_sampled": {}, "kimonu_det": {}}
for arm, sub in (("cocoa_false", "cocoa_v110"), ("cocoa_true", "cocoa_eff")):
    for m in D["models"]:
        fs = sorted(glob.glob(f"sweep/{sub}/{sub}_{m}_HiGHS_p64_seed*.json"))
        if len(fs) < 2: continue
        js = [json.load(open(f)) for f in fs]; free = free_set(m)
        labs = [[x if x >= 0 else None for x in j["kinetic_module"]] for j in js]
        out[arm][m] = dict(seeds=len(js), pairs=pair_spread(labs), fixed=set_spread([rel(j["acr"], j["acrr"], free) for j in js]),
                           fixed_ratios=ratio_spread([(j["acr"], j["acrr"]) for j in js], free))
root = "reference_final_sampled/acr/Results"
for m in D["ref_sampled"]:
    labs, sets, runs, free = [], [], [], free_set(m)
    for rec in D["ref_sampled"][m]:
        st = f"{m}_seed{rec['seed']}"
        mods = [x for x in json.load(open(f"{root}/modules_{st}.json"))["modules"] if x]
        l = [None] * (max(max(x) for x in mods) + 1)
        for i, x in enumerate(mods):
            for c in x: l[c] = i
        labs.append(l)
        import csv
        single = [r[1] for r in list(csv.reader(open(f"{root}/MetSingle/MetSingle_{st}.csv"), delimiter=" "))[1:]]
        double = read_pairs(f"{root}/MetDouble/MetDouble_{st}.csv")
        sets.append(rel(single, double, free)); runs.append((single, double))
    L = max(len(l) for l in labs); labs = [l + [None] * (L - len(l)) for l in labs]
    out["ref_sampled"][m] = dict(seeds=len(labs), pairs=pair_spread(labs), fixed=set_spread(sets), fixed_ratios=ratio_spread(runs, free))
def kimonu_run(f):
    d = json.load(open(f)); g = {}
    for c in d["complexes"]:
        g.setdefault(c["kinetic_module"], set()).add(c["id"])
    return {frozenset(x) for x in g.values()}, set(d["acr"]), {frozenset(p) for p in d["acrr"]}
for m in D["models"]:
    fs = sorted(glob.glob(f"{DET}/kimonu_{m}_rep*.json"))
    runs = [kimonu_run(f) for f in fs]
    assert len(runs) == 3 and all(r == runs[0] for r in runs), (m, DET, len(runs))
    out["kimonu_det"][m] = dict(seeds=len(runs), pairs=(0, None), fixed=(0, None), fixed_ratios=(0, None))
json.dump(out, open("figures/seed_spread.json", "w"), indent=1)
for arm, v in out.items():
    nz_p = sum(x["pairs"][0] > 0 for x in v.values()); nz_f = sum(x["fixed"][0] > 0 for x in v.values())
    nz_r = sum(x["fixed_ratios"][0] > 0 for x in v.values())
    print(f"{arm:12s} {len(v):2d} models | pairs seed-dependent in {nz_p:2d} | ACR+ACRR seed-dependent in {nz_f:2d} "
          f"(as ratios {nz_r:2d})")
