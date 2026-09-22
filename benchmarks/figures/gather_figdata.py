"""Build figdata.json (schema in make_prelim.py) from the benchmark result directories.

Arms: kimonu from sweep_shared/ (giant_overlap from the sweep_tentative/ logs), COCOA from sweep/cocoa_v110 (efficient=false) and sweep/cocoa_eff
(efficient=true), reference from reference_final/ (exhaustive, fixed algorithm and bridge; retested partitions from
reference_final/retest/ where present) and reference_final_sampled/. Earlier reference runs are not read: they predate the concordance fixes.
Free metabolites are those listed before "E_1" in the model's own metabolite list.

    python figures/gather_figdata.py [--check]   # --check diffs kimonu/COCOA arms against the committed file
"""
import csv, glob, json, os, re, sys
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(B)

def free_set(m):
    mets = json.load(open(f"unblocked/{m}_fixed.crn.json"))["metabolites"]
    cut = mets.index("E_1") if "E_1" in mets else len(mets)
    return set(mets[:cut])

def n_free(items, free):
    """items: metabolite names (ACR) or name pairs (ACRR); a pair is free when both members are."""
    return sum(all(x in free for x in (it if isinstance(it, list) else [it])) for it in items)

def kimonu_rec(m):
    d = json.load(open(f"sweep_shared/kimonu_{m}.json")); free = free_set(m)
    return dict(balanced=d["n_balanced"], conc=d["n_concordance_modules"], giant=d["giant_kinetic_module"],
                acr_free=n_free(d["acr_free"], free), acrr_free=n_free(d["acrr_free"], free),
                acr=len(d["acr"]), acrr=len(d["acrr"]), seconds=d["seconds"], sec_conc=d["seconds_concordance"],
                sec_kin=d["seconds_kinetic"], sec_acr=d["seconds_acr"], giant_overlap=giant_overlap(m),
                n_kinetic=d["n_kinetic_modules"], **s210(m)), d["n_complexes"]

def s210(m):
    """ACR/ACRR by the reference's combinatorial rule on Kimonu's modules (kinetic_variants/s210_criterion.py)."""
    f = f"kinetic_variants/s210_{m}.json"
    if not os.path.isfile(f): return {}
    d = json.load(open(f))
    return dict(acr_s210=len(d["acr"]), acrr_s210=len(d["acrr"]), acr_s210_free=d["acr_free"], acrr_s210_free=d["acrr_free"])

def giant_overlap(m):
    """Largest module under overlap-only merging (the reference's rule), from the sweep_tentative/ logs; the four-phase and two-phase upstream sets are identical, see kinetic_variants/."""
    for line in open(f"sweep_tentative/{m}.log"):
        if "after overlap merge" in line:
            return int(re.search(r"largest\s+(\d+)", line).group(1))
    raise ValueError(m)

def cocoa_recs(m, sub, tag):
    out = []
    for f in sorted(glob.glob(f"sweep/{sub}/{tag}_{m}_HiGHS_p64_seed*.json")):
        d = json.load(open(f)); free = free_set(m)
        cm = d["concordance_module"]; balanced = sum(1 for x in cm if x == 0)
        sizes = {}
        for x in cm:
            if x > 0: sizes[x] = sizes.get(x, 0) + 1
        ksizes = {}
        for x in d["kinetic_module"]:
            if x >= 0: ksizes[x] = ksizes.get(x, 0) + 1
        out.append(dict(seed=d["seed"], balanced=balanced, conc=sum(1 for v in sizes.values() if v > 1),
                        n_kinetic=sum(1 for v in ksizes.values() if v > 1),
                        giant=d["giant_module"], acr=len(d["acr"]), acrr=len(d["acrr"]),
                        acr_free=n_free(d["acr"], free), acrr_free=n_free(d["acrr"], free), seconds=d["seconds"]))
    return out

def ref_counts(d):
    mods = {}
    for c in d["complexes"]:
        if c["module"] > 0: mods[c["module"]] = mods.get(c["module"], 0) + 1
    return sum(1 for c in d["complexes"] if c["module"] == 0), sum(1 for v in mods.values() if v > 1)

def read_pairs(path):
    rows = list(csv.reader(open(path), delimiter=" "))[1:]
    vals = [r[1:] for r in rows]
    if vals and len(vals[0]) == 1:          # a single pair is written as a two-row column
        return [[vals[0][0], vals[1][0]]] if len(vals) == 2 else []
    return vals

def ref_exh(m):
    """Exhaustive reference; where the rejected pairs were re-tested with the scaling-bound fix
    (reference_final/retest/), the retested partition and its R stage are used and `retest` says so."""
    f = f"reference_final/reference_{m}.json"
    if not os.path.isfile(f): return None
    stem, rec = f"{m}_concordant_fixed", {}
    rt = f"reference_final/retest/reference_{m}_retest.json"
    if os.path.isfile(rt):
        f, stem = rt, f"{m}_concordant_fixed_retest"
        w = open(f"reference_final/retest/retest_{m}.txt").readline().split()
        rec["retest"] = dict(retested=int(w[1]), accepted=int(w[3]), residual=int(w[5]))
    d = json.load(open(f)); b, c = ref_counts(d)
    rec.update(balanced=b, conc=c, seconds=d["seconds"], sec_conc=d["seconds_concordant"])
    ov = f"reference_final/acr/Results/Overview/Overview_{stem}.csv"
    if os.path.isfile(ov):
        v = {int(r[0]): r[1] for r in list(csv.reader(open(ov), delimiter=" "))[1:]}
        free = free_set(m)
        single = [r[1] for r in list(csv.reader(open(f"reference_final/acr/Results/MetSingle/MetSingle_{stem}.csv"), delimiter=" "))[1:]]
        double = read_pairs(f"reference_final/acr/Results/MetDouble/MetDouble_{stem}.csv")
        mj = f"reference_final/acr/Results/modules_{stem}.json"
        if os.path.isfile(mj):
            md = json.load(open(mj))
            rec["n_kinetic"] = sum(1 for x in md["modules"] if len(x) > 1)     # Overview[4]-[5] is off by the late singleton append
        rec.update(giant=int(v[6]), acr=int(v[14]), acrr=int(v[15]),
                   acr_free=n_free(single, free), acrr_free=n_free(double, free))
        assert len(single) == rec["acr"] and len(double) == rec["acrr"], (m, len(single), len(double), rec)
    return rec

def ref_sampled(m):
    out = []
    for f in sorted(glob.glob(f"reference_final_sampled/reference_{m}_seed*.json")):
        d = json.load(open(f)); b, c = ref_counts(d)
        out.append(dict(seed=int(re.search(r"seed(\d+)", f).group(1)), balanced=b, conc=c, seconds=d["seconds"]))
    return out

if __name__ == "__main__":
    old = json.load(open("figdata.json"))
    models = old["models"]
    D = dict(models=models, n_complexes={}, kimonu={}, cocoa_false={}, cocoa_true={}, ref_exh={}, ref_sampled={})
    for m in models:
        D["kimonu"][m], D["n_complexes"][m] = kimonu_rec(m)
        cf, ct = cocoa_recs(m, "cocoa_v110", "cocoa_v110"), cocoa_recs(m, "cocoa_eff", "cocoa_eff")
        if cf: D["cocoa_false"][m] = cf
        if ct: D["cocoa_true"][m] = ct
        r = ref_exh(m)
        if r: D["ref_exh"][m] = r
        s = ref_sampled(m)
        if s: D["ref_sampled"][m] = s

    if "--check" in sys.argv:
        for arm in ("n_complexes", "kimonu", "cocoa_false", "cocoa_true"):
            for m in models:
                a, b = old[arm].get(m), D[arm].get(m)
                if isinstance(a, list): a = sorted(a, key=lambda r: r["seed"]); b = sorted(b or [], key=lambda r: r["seed"])
                if a != b: print("DIFF", arm, m, "\n  old", a, "\n  new", b)
        print("check done")
    else:
        json.dump(D, open("figdata.json", "w"), indent=1)
        print({k: len(v) for k, v in D.items() if isinstance(v, dict)})
