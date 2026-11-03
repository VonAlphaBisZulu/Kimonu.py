"""Pair counts for Figure 3B on the models where the exhaustive reference finished: coupled complex pairs
(sum of C(size, 2) over the kinetic modules) and ACR metabolites plus ACRR pairs among free metabolites, per arm.

Sources: Kimonu default (iterated route, detached uptakes) from kinetic_variants/uptake_modes/ (detached, it); the
reference-equivalent mode from a rerun of run_kimonu_refmode.py that also writes its modules (KIMONU_REFMOD_DIR);
the reference from reference_final/acr/Results/modules_<stem>.json (exhaustive) and reference_final_sampled/acr/
Results/modules_<model>_seed<seed>.json (sampled); COCOA.jl from the kinetic_module labels in sweep/cocoa_v110 and
sweep/cocoa_eff. ACR and ACRR counts are those of figdata.json. Output: figures/pair_counts.json.

    python figures/pair_counts.py     # from benchmarks/
"""
import glob, json, os, re
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(B)
REFMOD = os.environ.get("KIMONU_REFMOD_DIR", "/scratch/sadreyev_lab/pschneider/kimonu_trip_2026-11-02/refmode_modules")
D = json.load(open("figdata.json")); REF = sorted(D["ref_exh"], key=lambda m: D["n_complexes"][m])
pairs = lambda sizes: sum(s * (s - 1) // 2 for s in sizes)

def labels(lab):
    g = {}
    for x in lab:
        if x >= 0: g[x] = g.get(x, 0) + 1
    return g.values()

def fixed(rec):
    return rec["acr_free"] + rec["acrr_free"]

out = {}
for m in REF:
    r = out[m] = {}
    f = f"kinetic_variants/uptake_modes/{m}.json"; det = f.replace(".json", ".detached.json")
    it = (json.load(open(det)) if os.path.exists(det) else json.load(open(f)))["detached"]["it"]
    r["kimonu_det"] = dict(pairs=it["pairs"], fixed=it["acr_free"] + it["acrr_pairs_free"])
    assert (it["acr_free"], it["acrr_pairs_free"]) == (D["kimonu_det"][m]["acr_free"], D["kimonu_det"][m]["acrr_free"]), m
    k = json.load(open(f"{REFMOD}/kimonu_{m}.json"))
    assert all(k[x] == D["kimonu_ref"][m][x] for x in ("n_kinetic", "giant", "acr_free", "acrr_free")), m
    r["kimonu_ref"] = dict(pairs=pairs(len(x) for x in k["modules"]), fixed=fixed(k))
    stem = next(os.path.basename(p)[8:-5] for p in glob.glob(f"reference_final/acr/Results/modules_{m}_concordant_fixed*.json"))
    mods = json.load(open(f"reference_final/acr/Results/modules_{stem}.json"))["modules"]
    assert sum(len(x) > 1 for x in mods) == D["ref_exh"][m]["n_kinetic"], m
    r["ref_exh"] = dict(pairs=pairs(len(x) for x in mods), fixed=fixed(D["ref_exh"][m]))
    r["ref_sampled"] = []
    for rec in sorted(D["ref_sampled"].get(m, []), key=lambda x: x["seed"]):
        mods = json.load(open(f"reference_final_sampled/acr/Results/modules_{m}_seed{rec['seed']}.json"))["modules"]
        assert sum(len(x) > 1 for x in mods) == rec["n_kinetic"], (m, rec["seed"])
        r["ref_sampled"].append(dict(seed=rec["seed"], pairs=pairs(len(x) for x in mods), fixed=fixed(rec)))
    for arm, sub in (("cocoa_false", "cocoa_v110"), ("cocoa_true", "cocoa_eff")):
        r[arm] = []
        for rec in sorted(D[arm].get(m, []), key=lambda x: x["seed"]):
            d = json.load(open(f"sweep/{sub}/{sub}_{m}_HiGHS_p64_seed{rec['seed']}.json"))
            g = labels(d["kinetic_module"]); assert sum(s > 1 for s in g) == rec["n_kinetic"], (m, arm)
            r[arm].append(dict(seed=rec["seed"], pairs=pairs(g), fixed=fixed(rec)))
    # the reference-equivalent mode reproduces the exhaustive reference, pairs included
    assert r["kimonu_ref"] == r["ref_exh"], (m, r["kimonu_ref"], r["ref_exh"])
json.dump(out, open("figures/pair_counts.json", "w"), indent=1)
for m in REF:
    med = lambda v: sorted(x["pairs"] for x in v)[len(v) // 2] if isinstance(v, list) else v["pairs"]
    medf = lambda v: sorted(x["fixed"] for x in v)[len(v) // 2] if isinstance(v, list) else v["fixed"]
    print(f"{m[:16]:16s}", " ".join(f"{a}:{med(v)}/{medf(v)}" for a, v in out[m].items()))
