"""Figure 3B/C on one counting rule, and set containment against Kimonu, on the models where the exhaustive
reference finished.

Figure 3C counts independent fixed relations: one per ACR metabolite and k - 1 per ACRR group of k metabolites,
which does not depend on how a tool lists its pairs ("fixed"; "fixed_pairs" is ACR plus ACRR pairs). Before
counting, each tool's ACR and ACRR output is closed under what it implies, then counted among free
metabolites. A metabolite in a fixed ratio with an ACR metabolite has ACR itself, and fixed ratios are transitive,
so ACRR groups are closed. ACRR pairs are counted within groups of metabolites without ACR: a pair of two ACR
metabolites is fixed only because each concentration is. Kimonu's output is closed already; COCOA.jl and the
reference list some implied relations (pairs of ACR metabolites) and leave others implicit.

Containment: complexes are matched across tools by composition (their column of Y), metabolites by name. For
every arm, the coupled complex pairs (pairs within one kinetic module), the free ACR metabolites and the free
ACRR pairs (closed under the rule above) are compared with Kimonu's default (iterated route, detached uptakes),
and reported as: contained in Kimonu's, and how many are not.

    python figures/pair_sets.py     # from benchmarks/  -> figures/pair_sets.json; updates "fixed" in figures/pair_counts.json
"""
import glob, json, os, re, sys
import numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(B)
sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
import kimonu
from kimonu.acr import free_directions, readout
from kimonu.io import zero_complexes
from run_kimonu_crn import triplets_to_sparse, to_exact
from kinetic_variants.s210_criterion import s210
from figures.gather_figdata import read_pairs
import csv

REFMOD = os.environ.get("KIMONU_REFMOD_DIR", "/scratch/sadreyev_lab/pschneider/kimonu_trip_2026-11-02/refmode_modules")
D = json.load(open("figdata.json")); REF = sorted(D["ref_exh"], key=lambda m: D["n_complexes"][m])
key = lambda comp: frozenset((k, float(v)) for k, v in comp.items() if v)


def cpairs(modules):
    """coupled complex pairs, as frozensets of two composition keys"""
    return {frozenset((a, b)) for mod in modules for i, a in enumerate(mod) for b in list(mod)[i + 1:] if a != b}


def rule(acr, acrr, free):
    """Close ACR and ACRR (all metabolites) under implication; return free ACR names, free ACRR pairs within
    non-ACR groups, and how many listed free pairs were pairs of two ACR metabolites."""
    parent = {}
    find = lambda a: a if parent.setdefault(a, a) == a else (parent.__setitem__(a, find(parent[a])) or parent[a])
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    FIX = ("fixed",)
    for a in acr: union(a, FIX)
    for p in acrr:
        p = list(p)
        if len(set(p)) == 2: union(p[0], p[1])
    groups = {}
    for a in list(parent):
        if a != FIX: groups.setdefault(find(a), set()).add(a)
    fix = find(FIX)
    acr_c = {a for a in groups.get(fix, ()) if a in free}
    pairs = {frozenset((x, y)) for r, g in groups.items() if r != fix
             for x in g for y in g if x < y and x in free and y in free}
    listed = {frozenset(p) for p in acrr if len(set(p)) == 2 and all(x in free for x in p)}
    return acr_c, pairs, len({p for p in listed if all(x in acr_c for x in p)})


def independent(acr, pairs):
    """Independent fixed relations: one per ACR metabolite and k - 1 per ACRR group of k metabolites."""
    parent = {}
    def find(a):
        while parent.setdefault(a, a) != a:
            a = parent[a]
        return a
    for p in pairs:
        x, y = sorted(p)
        rx, ry = find(x), find(y)
        if rx != ry: parent[rx] = ry
    nodes = set(parent)
    return len(acr) + len(nodes) - len({find(a) for a in nodes})


def cocoa_key(s, mets):
    if s in ("", "0", "∅"): return frozenset()
    comp = {}
    for t in s.split("+"):
        mm = re.match(r"^(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)_(.+)$", t)
        c, name = (float(mm.group(1)), mm.group(2)) if mm and mm.group(2) in mets and t not in mets else (1.0, t)
        assert name in mets, (s, t)
        comp[name] = comp.get(name, 0) + c
    return key(comp)


out, counts = {}, json.load(open("figures/pair_counts.json"))
for m in REF:
    crn = json.load(open(f"unblocked/{m}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64); Y = triplets_to_sparse(crn["Y"], nm, nc)
    N = to_exact(triplets_to_sparse(crn["S"], nm, nr))
    mets = crn["metabolites"]; cut = mets.index("E_1") if "E_1" in mets else nm; free = set(mets[:cut]); mset = set(mets)
    Yc = Y.tocsc()
    ckey = [key({mets[i]: v for i, v in zip(Yc.indices[Yc.indptr[c]:Yc.indptr[c + 1]], Yc.data[Yc.indptr[c]:Yc.indptr[c + 1]])}) for c in range(nc)]
    skip = zero_complexes(Y)
    K = lambda mods: [[ckey[c] for c in mod if c not in skip] for mod in mods if len(set(mod) - skip) > 1]
    arms = {}
    # Kimonu default
    km = kimonu.kinetic_modules(A, Y, N, detach_uptakes=True)
    R = readout(Y, free_directions(Y, km.modules), skip)
    acrr = [(mets[a], mets[b]) for g in R.acrr for i, a in enumerate(g) for b in g[i + 1:]]
    arms["kimonu_det"] = [(cpairs(K(km.modules)),) + rule([mets[i] for i in R.acr], acrr, free)]
    # Kimonu reference-equivalent mode: its modules, ACR/ACRR by the reference's pairwise rule
    k = json.load(open(f"{REFMOD}/kimonu_{m}.json"))
    a, r = s210(Y, k["modules"])
    arms["kimonu_ref"] = [(cpairs(K(k["modules"])),) + rule([mets[i] for i in a], [(mets[i], mets[j]) for i, j in r], free)]
    # reference, exhaustive and sampled: complexes in the reference's own list, matched by composition
    def ref_arm(rec_json, root, stem):
        comps = [key(c["composition"]) for c in json.load(open(rec_json))["complexes"]]
        mods = json.load(open(f"{root}/Results/modules_{stem}.json"))["modules"]
        # module indices are 0-based into the reference's complex list (checked against the reference-equivalent mode)
        km_ = [[comps[i] for i in mod] for mod in mods if len(mod) > 1]
        single = [row[1] for row in list(csv.reader(open(f"{root}/Results/MetSingle/MetSingle_{stem}.csv"), delimiter=" "))[1:]]
        double = [tuple(p) for p in read_pairs(f"{root}/Results/MetDouble/MetDouble_{stem}.csv")]
        return (cpairs(km_),) + rule(single, double, free)
    stem = next(os.path.basename(p)[8:-5] for p in glob.glob(f"reference_final/acr/Results/modules_{m}_concordant_fixed*.json"))
    rt = f"reference_final/retest/reference_{m}_retest.json"
    arms["ref_exh"] = [ref_arm(rt if "retest" in stem else f"reference_final/reference_{m}.json", "reference_final/acr", stem)]
    arms["ref_sampled"] = [ref_arm(f, "reference_final_sampled/acr", f"{m}_seed{re.search(r'seed(\d+)', f).group(1)}")
                           for f in sorted(glob.glob(f"reference_final_sampled/reference_{m}_seed*.json"), key=lambda f: int(re.search(r"seed(\d+)", f).group(1)))]
    for arm, sub in (("cocoa_false", "cocoa_v110"), ("cocoa_true", "cocoa_eff")):
        arms[arm] = []
        for f in sorted(glob.glob(f"sweep/{sub}/{sub}_{m}_HiGHS_p64_seed*.json"), key=lambda f: int(re.search(r"seed(\d+)", f).group(1))):
            d = json.load(open(f)); keys = [cocoa_key(s, mset) for s in d["complex_ids"]]
            g = {}
            for c, lab in enumerate(d["kinetic_module"]):
                if lab >= 0 and keys[c]: g.setdefault(lab, []).append(keys[c])
            arms[arm].append((cpairs([v for v in g.values() if len(v) > 1]),) + rule(d["acr"], [tuple(p) for p in d["acrr"]], free))
    ref_pairs, ref_acr, ref_acrr, _ = arms["kimonu_det"][0]
    assert arms["kimonu_ref"][0][0] == arms["ref_exh"][0][0], f"{m}: reference-equivalent pairs differ from the reference's (index offset?)"
    out[m] = {}
    for arm, recs in arms.items():
        out[m][arm] = [dict(pairs=len(p), fixed=independent(a, r), fixed_pairs=len(a) + len(r), acr=len(a), acrr=len(r), acr_acr_pairs_dropped=n,
                            pairs_not_in_kimonu=len(p - ref_pairs), acr_not_in_kimonu=sorted(a - ref_acr),
                            acrr_not_in_kimonu=[sorted(x) for x in r - ref_acrr if not all(y in ref_acr for y in x)],
                            kimonu_pairs_missing=len(ref_pairs - p))
                       for p, a, r, n in recs]
    print(m, {arm: [(v["pairs"], v["fixed"], v["pairs_not_in_kimonu"], len(v["acr_not_in_kimonu"]) + len(v["acrr_not_in_kimonu"]), v["acr_acr_pairs_dropped"]) for v in recs]
              for arm, recs in out[m].items()}, flush=True)
    # one rule for Figure 3C: replace "fixed" in pair_counts.json, keep the old value
    for arm, v in counts[m].items():
        new = out[m][arm]
        # match by seed order: pair_counts lists are sorted by seed, as are the files read here
        if isinstance(v, list):
            for x, y in zip(sorted(v, key=lambda r: r["seed"]), new):
                x.setdefault("fixed_as_reported", x["fixed"]); x["fixed"] = y["fixed"]
        else:
            v.setdefault("fixed_as_reported", v["fixed"]); v["fixed"] = new[0]["fixed"]
json.dump(out, open("figures/pair_sets.json", "w"), indent=1)
json.dump(counts, open("figures/pair_counts.json", "w"), indent=1)
