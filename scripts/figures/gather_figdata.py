"""Build figdata.json from the benchmark result directories. Every figure of the benchmark and Data S1 read it.

All paths are paths of the records (records.py): read from KIMONU_RECORDS, figdata.json written to
KIMONU_OUT/figures/source_data/. Arms: kimonu (iterated route, attached uptakes; seconds the median of three seeds),
kimonu_det (iterated route, detached uptakes, Kimonu's default) and kimonu_ref (reference-equivalent mode,
run_kimonu_refmode.py) from the timing arms att, det and ref (results/kimonu/timing/, the rerun of 2026-11-11 by
timing_exclusive.sbatch), kimonu_up from results/kimonu/upstream_set_route_2026-09/kimonu_runs/ (the upstream-set
route of September, kept in the file but drawn by no figure of the paper; compat/upstream.py; giant_overlap from its
overlap_merge_logs/, the pairwise-rule fields *_s210 from its shinar_feinberg_rule/), COCOA from
results/cocoa_jl/efficient_false and efficient_true, reference from results/reference/exhaustive/ (fixed algorithm and
bridge; retested partitions from its retest/ where present) and results/reference/sampled/. Earlier reference runs are
not read: they predate the concordance fixes. Free metabolites are those listed before "E_1" in the model's own
metabolite list.

    python figures/gather_figdata.py [--check]   # --check diffs kimonu/COCOA arms against the figdata.json of the records

Complexes are counted as in the deficiency, the zero complex included and the complexes that occur in no reaction left
out. The crn.json lists 51 such complexes for iJB785 and 29 for iSynCJ816. Kimonu and the reference count them as
balanced, COCOA.jl leaves them out (and the zero complex as well), so every balanced count here leaves them out, checked
against the per-complex records of each tool.

figdata.json schema:
  models, n_complexes[m]   (complexes as above), complexes_in_no_reaction[m]   (left out of n_complexes and balanced)
  kimonu[m], kimonu_det[m], kimonu_ref[m]: balanced, conc, n_kinetic, giant, acr, acrr, acr_free, acrr_free,
      acrr_ratios, acrr_ratios_free, seconds, sec_conc, sec_kin, sec_acr
  cocoa_false[m] / cocoa_true[m]: [{seed, balanced, conc, n_kinetic, giant, acr, acrr, acr_free, acrr_free, acrr_ratios,
      acrr_ratios_free, seconds}, ...]
  ref_exh[m]: balanced, conc, seconds, sec_conc (+ n_kinetic, giant, acr, acrr, acr_free, acrr_free, acrr_ratios,
      acrr_ratios_free where the R stage ran; retest where rejected pairs were retested)
  ref_sampled[m]: [{seed, balanced, conc, seconds, ...}, ...]
acrr and acrr_free count the ACRR pairs as each tool lists them. acrr_ratios and acrr_ratios_free count independent
ratios, k - 1 per ACRR group of k metabolites (free: k free metabolites of a group), with the groups joined from the
listed pairs, so they do not depend on how a tool lists its pairs. The paper counts ACRR this way.
A network a tool did not finish within 7 days is missing from that arm.
"""
import csv, glob, json, os, re, sys
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B)
import records
TIMING = records.read(records.TIMING)
DET_DIR, REF_DIR, ATT_DIR = f"{TIMING}/det", f"{TIMING}/ref", f"{TIMING}/att"
UP = records.read(records.UPSTREAM_ROUTE)
COCOA = {"cocoa_v110": records.read(records.COCOA_FALSE), "cocoa_eff": records.read(records.COCOA_TRUE)}
REF_EXH, REF_SAMPLED = records.read(records.REF_EXHAUSTIVE), records.read(records.REF_SAMPLED)
NO_REACTION = {}


def complexes(m):
    """The number of complexes as the paper counts them (the zero complex included, the complexes that occur in no
    reaction left out) and the ids of those left out, complex_<row of A, from 1> as every tool's record names them."""
    crn = json.load(open(records.network(m)))
    used = {i for i, _, _ in crn["A"]}
    NO_REACTION[m] = {f"complex_{c}" for c in range(1, crn["n_complexes"] + 1) if c not in used}
    return len(used), crn["n_complexes"]

def free_set(m):
    mets = json.load(open(records.network(m)))["metabolites"]
    cut = mets.index("E_1") if "E_1" in mets else len(mets)
    return set(mets[:cut])

def n_free(items, free):
    """items: metabolite names (ACR) or name pairs (ACRR); a pair is free when both members are."""
    return sum(all(x in free for x in (it if isinstance(it, list) else [it])) for it in items)

def ratios(pairs, free=None):
    """Independent ACRR ratios, k - 1 per group of k metabolites, the groups joined from the listed pairs. With free,
    each group counts only its free metabolites (k free metabolites of a group give k - 1 free ratios)."""
    parent = {}
    def find(a):
        while parent.setdefault(a, a) != a:
            a = parent[a]
        return a
    for p in pairs:
        a, b = list(p)
        parent[find(a)] = find(b)
    groups = {}
    for a in list(parent):
        groups.setdefault(find(a), set()).add(a)
    return sum(max(len(g if free is None else g & free) - 1, 0) for g in groups.values())

def kimonu_rec(m, root=ATT_DIR):
    """psi route; the replicates agree on every count, so only the times vary."""
    reps = [json.load(open(f)) for f in sorted(glob.glob(f"{root}/kimonu_{m}_rep*.json"))]
    d, free = reps[0], free_set(m)
    med = lambda k: sorted(r[k] for r in reps)[len(reps) // 2]
    n, n_crn = complexes(m)
    assert d["n_complexes"] == len(d["complexes"]) == n_crn - 1, m          # Kimonu lists all but the zero complex
    assert all(c["module"] == 0 for c in d["complexes"] if c["id"] in NO_REACTION[m]), m   # and calls these balanced
    return dict(balanced=d["n_balanced"] - len(NO_REACTION[m]), conc=d["n_concordance_modules"], giant=d["giant_kinetic_module"],
                acr_free=n_free(d["acr_free"], free), acrr_free=n_free(d["acrr_free"], free),
                acr=len(d["acr"]), acrr=len(d["acrr"]), acrr_ratios=ratios(d["acrr"]), acrr_ratios_free=ratios(d["acrr"], free),
                seconds=med("seconds"), sec_conc=med("seconds_concordance"),
                sec_kin=med("seconds_kinetic"), sec_acr=med("seconds_acr"), n_kinetic=d["n_kinetic_modules"]), n


def kimonu_ref_rec(m):
    """reference-equivalent mode (run_kimonu_refmode.py): the reference's definitions from null spaces. Its records
    hold counts and modules, so the ACRR pairs for the ratios come from the Shinar-Feinberg rule on those modules."""
    reps = [json.load(open(f)) for f in sorted(glob.glob(f"{REF_DIR}/kimonu_{m}_rep*.json"))]
    if not reps: return None
    d = reps[0]; med = lambda k: sorted(r[k] for r in reps)[len(reps) // 2]
    sys.path[:0] = [B, os.path.join(os.path.dirname(B), "src")]
    import kimonu
    from kimonu.reference import shinar_feinberg
    from run_kimonu_crn import triplets_to_sparse
    crn = json.load(open(records.network(m)))
    mets = crn["metabolites"]; free = free_set(m)
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    net = kimonu.crn_from_matrices(triplets_to_sparse(crn["A"], nc, nr, dtype=int), triplets_to_sparse(crn["Y"], nm, nc))
    acr, acrr = shinar_feinberg(net, d["modules"])
    pairs = [(mets[i], mets[j]) for i, j in acrr]
    assert (len(acr), len(pairs), n_free(list(map(list, pairs)), free)) == (d["acr"], d["acrr"], d["acrr_free"]), m
    return dict(balanced=d["balanced"] - len(NO_REACTION[m]), conc=d["conc"], n_kinetic=d["n_kinetic"], giant=d["giant"],
                acr=d["acr"], acrr=d["acrr"], acr_free=d["acr_free"], acrr_free=d["acrr_free"],
                acrr_ratios=ratios(pairs), acrr_ratios_free=ratios(pairs, free),
                seconds=med("seconds"), sec_conc=med("seconds_concordance"), sec_kin=med("seconds_kinetic"), sec_acr=med("seconds_acr"))

def kimonu_up_rec(m):
    d = json.load(open(f"{UP}/kimonu_runs/kimonu_{m}.json")); free = free_set(m)
    assert all(c["module"] == 0 for c in d["complexes"] if c["id"] in NO_REACTION[m]), m
    return dict(balanced=d["n_balanced"] - len(NO_REACTION[m]), conc=d["n_concordance_modules"], giant=d["giant_kinetic_module"],
                acr_free=n_free(d["acr_free"], free), acrr_free=n_free(d["acrr_free"], free),
                acr=len(d["acr"]), acrr=len(d["acrr"]), seconds=d["seconds"], sec_conc=d["seconds_concordance"],
                sec_kin=d["seconds_kinetic"], sec_acr=d["seconds_acr"], giant_overlap=giant_overlap(m),
                n_kinetic=d["n_kinetic_modules"], **s210(m))

def s210(m):
    """ACR/ACRR by the reference's pairwise rule on the September modules, from the records of the removed
    kinetic_variants/s210_criterion.py (HISTORY.md); the rule is now kimonu.reference.shinar_feinberg."""
    f = f"{UP}/shinar_feinberg_rule/s210_{m}.json"
    if not os.path.isfile(f): return {}
    d = json.load(open(f))
    return dict(acr_s210=len(d["acr"]), acrr_s210=len(d["acrr"]), acr_s210_free=d["acr_free"], acrr_s210_free=d["acrr_free"])

def giant_overlap(m):
    """Largest module under overlap-only merging (the reference's rule), from the overlap_merge_logs/; the four-phase and two-phase upstream sets are identical, see kinetic_variants/."""
    for line in open(f"{UP}/overlap_merge_logs/{m}.log"):
        if "after overlap merge" in line:
            return int(re.search(r"largest\s+(\d+)", line).group(1))
    raise ValueError(m)

def cocoa_recs(m, sub, tag):
    out = []
    for f in sorted(glob.glob(f"{COCOA[sub]}/{tag}_{m}_HiGHS_p64_seed*.json")):
        d = json.load(open(f)); free = free_set(m)
        assert d["n_complexes"] == len(d["complex_ids"]) == complexes(m)[0] - 1, m   # neither the zero complex nor those of no reaction
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
                        acr_free=n_free(d["acr"], free), acrr_free=n_free(d["acrr"], free),
                        acrr_ratios=ratios(d["acrr"]), acrr_ratios_free=ratios(d["acrr"], free), seconds=d["seconds"]))
    return out

def ref_counts(d, m):
    """balanced complexes and concordance modules of two or more; the reference lists every complex of the crn.json,
    the zero complex (module -1) and those of no reaction (module 0, balanced) included"""
    assert len(d["complexes"]) == complexes(m)[1], m
    mods = {}
    for c in d["complexes"]:
        if c["module"] > 0: mods[c["module"]] = mods.get(c["module"], 0) + 1
    return sum(1 for c in d["complexes"] if c["module"] == 0 and c["id"] not in NO_REACTION[m]), sum(1 for v in mods.values() if v > 1)

def read_pairs(path):
    rows = list(csv.reader(open(path), delimiter=" "))[1:]
    vals = [r[1:] for r in rows]
    if vals and len(vals[0]) == 1:          # a single pair is written as a two-row column
        return [[vals[0][0], vals[1][0]]] if len(vals) == 2 else []
    return vals

def ref_exh(m):
    """Exhaustive reference; where the rejected pairs were re-tested with the scaling-bound fix
    (retest/), the retested partition and its R stage are used and `retest` says so."""
    f = f"{REF_EXH}/reference_{m}.json"
    if not os.path.isfile(f): return None
    stem, rec = f"{m}_concordant_fixed", {}
    rt = f"{REF_EXH}/retest/reference_{m}_retest.json"
    if os.path.isfile(rt):
        f, stem = rt, f"{m}_concordant_fixed_retest"
        w = open(f"{REF_EXH}/retest/retest_{m}.txt").readline().split()
        rec["retest"] = dict(retested=int(w[1]), accepted=int(w[3]), residual=int(w[5]))
    d = json.load(open(f)); b, c = ref_counts(d, m)
    rec.update(balanced=b, conc=c, seconds=d["seconds"], sec_conc=d["seconds_concordant"])
    return r_stage(rec, m, f"{REF_EXH}/acr", stem)

def r_stage(rec, m, root, stem):
    """Add giant, kinetic-module count and ACR/ACRR from a reference R-stage run when present."""
    ov = f"{root}/Results/Overview/Overview_{stem}.csv"
    if not os.path.isfile(ov): return rec
    v = {int(r[0]): r[1] for r in list(csv.reader(open(ov), delimiter=" "))[1:]}
    free = free_set(m)
    single = [r[1] for r in list(csv.reader(open(f"{root}/Results/MetSingle/MetSingle_{stem}.csv"), delimiter=" "))[1:]]
    double = read_pairs(f"{root}/Results/MetDouble/MetDouble_{stem}.csv")
    mj = f"{root}/Results/modules_{stem}.json"
    if os.path.isfile(mj):
        rec["n_kinetic"] = sum(1 for x in json.load(open(mj))["modules"] if len(x) > 1)
    rec.update(giant=int(v[6]), acr=int(v[14]), acrr=int(v[15]), acr_free=n_free(single, free), acrr_free=n_free(double, free),
               acrr_ratios=ratios(double), acrr_ratios_free=ratios(double, free))
    assert len(single) == rec["acr"] and len(double) == rec["acrr"], (m, stem, len(single), len(double), rec)
    return rec

def ref_sampled(m):
    out = []
    for f in sorted(glob.glob(f"{REF_SAMPLED}/reference_{m}_seed*.json")):
        d = json.load(open(f)); b, c = ref_counts(d, m); seed = int(re.search(r"seed(\d+)", f).group(1))
        out.append(r_stage(dict(seed=seed, balanced=b, conc=c, seconds=d["seconds"]), m, f"{REF_SAMPLED}/acr", f"{m}_seed{seed}"))
    return out

if __name__ == "__main__":
    old = json.load(open(records.read(records.SOURCE_DATA, "figdata.json")))
    models = old["models"]
    D = dict(models=models, n_complexes={}, complexes_in_no_reaction={}, kimonu={}, kimonu_det={}, kimonu_up={}, kimonu_ref={},
             cocoa_false={}, cocoa_true={}, ref_exh={}, ref_sampled={})
    for m in models:
        D["kimonu"][m], D["n_complexes"][m] = kimonu_rec(m)
        D["complexes_in_no_reaction"][m] = len(NO_REACTION[m])
        D["kimonu_det"][m] = kimonu_rec(m, DET_DIR)[0]
        D["kimonu_up"][m] = kimonu_up_rec(m)
        r = kimonu_ref_rec(m)
        if r: D["kimonu_ref"][m] = r
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
        json.dump(D, open(records.write(records.SOURCE_DATA, "figdata.json"), "w"), indent=1)
        print({k: len(v) for k, v in D.items() if isinstance(v, dict)})
