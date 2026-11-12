"""Data S1 and Data S2 of the Supplement, as CSV files.

Data S1 (DataS1_values.csv) has one row per network, tool and mode, and for the stochastic tools one row per seed.
Every file has the model file id (model) and the name the paper uses (network, figures/names.py).
The values are those of figdata.json (gather_figdata.py). Kimonu's coupled complex pairs are given on all 33 networks:
the iterated route's from kinetic_variants/uptake_modes/ (it.pairs, detached and attached), those of the
reference-equivalent mode from the modules of its runs ($KIMONU_REF_DIR). Those of the other tools come from
figures/pair_counts.json (pair_sets.py), which has them for the exhaustive reference on the 11 networks that it
finished, and must agree with the above where it has Kimonu's. The tangent route (detached uptakes, Kimonu's default)
has its own rows. Its kinetic modules, ACR, ACRR and pairs are those of uptake_modes/ (tan) and
kinetic_variants/member_counts/ (the ACRR group sizes, tan, or it where the routes agree), checked against the timed
runs in $KIMONU_TIMING/tan, and its wall time is the median of those three runs (seeds 0, 1, 2, which agree on every
count). That route computes no balanced complexes or concordance modules, so these two cells stay empty.
A tool that did not finish within 7 days gets "not finished in 7 days".
ACRR is counted as independent ratios, k - 1 per group of k metabolites (acrr_ratios, and acrr_ratios_free for k free
metabolites of a group), as in the paper. The last two columns keep the ACRR pairs as each tool lists them.

Data S2 (DataS2_uptake_acr.csv) lists every free metabolite with ACR under attached uptakes, with its compartment and
what becomes of it under detached uptakes (ACR in both, in an ACRR group, or free). It reads the lists that
uptake_lists.py wrote (iterated route, both uptake modes, metabolites as indices into the network's own list).
Free metabolites are those before "E_1" in that list. The script prints the totals that Supplementary Text S5 quotes.

    python figures/dataS1_S2.py <out_dir>        # from scripts/

UPTAKE_LISTS replaces the lists of the paper, /scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/uptake_lists (e.g. a rerun),
KIMONU_TIMING and KIMONU_REF_DIR the timing runs (default the rerun of 2026-11-11, timing/ and timing/ref).
Complexes are counted as in table S2 (with the zero complex, without the complexes of no reaction), and so are the
balanced complexes of every tool (figdata.json, gather_figdata.py).
"""
import csv, glob, json, os, statistics, sys

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from names import network
LISTS = os.environ.get("UPTAKE_LISTS", "/scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/uptake_lists")
TIMING = os.environ.get("KIMONU_TIMING", "/scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/timing")
REF_DIR = os.environ.get("KIMONU_REF_DIR", f"{TIMING}/ref")

ARMS = [("kimonu_det", "Kimonu, iterated route, detached uptakes"),
        ("kimonu_tan", "Kimonu, tangent route, detached uptakes"),
        ("kimonu", "Kimonu, iterated route, attached uptakes"),
        ("kimonu_ref", "Kimonu, reference-equivalent mode"),
        ("ref_exh", "Reference, exhaustive (Langary et al., 2025)"),
        ("ref_sampled", "Reference, sampled (Langary et al., 2025)"),
        ("cocoa_false", "COCOA.jl v1.1.0, efficient=false"),
        ("cocoa_true", "COCOA.jl v1.1.0, efficient=true")]
KEYS = ["balanced", "conc", "n_kinetic", "giant", "acr", "acrr_ratios", "acr_free", "acrr_ratios_free", "seconds"]
LISTED = [("acrr", "acrr_pairs_as_listed"), ("acrr_free", "acrr_pairs_as_listed_free")]
NOT_FINISHED = "not finished in 7 days"


pairs_of = lambda sizes: sum(k * (k - 1) // 2 for k in sizes)
ratios = lambda sizes: sum(k - 1 for k in sizes)                 # k - 1 independent ratios per ACRR group of k


def kimonu_rows(m, D, PC):
    """Kimonu's arms of network m with the coupled complex pairs on every network, and the tangent route (see above)"""
    um = json.load(open(f"{B}/kinetic_variants/uptake_modes/{m}.json"))
    mc = json.load(open(f"{B}/kinetic_variants/member_counts/{m}.json"))
    out = {}
    for arm, mode in (("kimonu_det", "detached"), ("kimonu", "attached")):
        it, v = um[mode]["it"], D[arm][m]
        assert (v["n_kinetic"], v["giant"], v["acr"], v["acrr"], v["acr_free"], v["acrr_free"]) == \
            (it["modules"], it["giant"], it["acr"], it["acrr_pairs"], it["acr_free"], it["acrr_pairs_free"]), (m, arm)
        out[arm] = dict(v, pairs=it["pairs"])
    assert out["kimonu_det"]["pairs"] == pairs_of(mc["it"]["module_sizes"]), m
    crn = json.load(open(f"{B}/unblocked/{m}_fixed.crn.json"))
    zero = set(range(crn["n_complexes"])) - {j - 1 for _, j, _ in crn["Y"]}          # the zero complex (Y triplets from 1)
    ref = json.load(open(sorted(glob.glob(f"{REF_DIR}/kimonu_{m}_rep*.json"))[0]))  # the mode has no random point
    assert not zero & {c for s in ref["modules"] for c in s} and len(ref["modules"]) == D["kimonu_ref"][m]["n_kinetic"], m
    out["kimonu_ref"] = dict(D["kimonu_ref"][m], pairs=pairs_of(len(s) for s in ref["modules"]))
    for arm in ("kimonu_det", "kimonu_ref"):
        if arm in PC.get(m, {}):
            assert PC[m][arm]["pairs"] == out[arm]["pairs"], (m, arm)
    # tangent route
    tan, g = um["detached"]["tan"], mc["tan"] or mc["it"]
    reps = [json.load(open(f)) for f in sorted(glob.glob(f"{TIMING}/tan/tangent_{m}_rep*.json"))]
    assert len(reps) == 3 and all(r["detached"] for r in reps), m
    mets = crn["metabolites"]; free = set(range(mets.index("E_1") if "E_1" in mets else len(mets)))
    for r in reps:
        assert (r["n_modules"], r["giant_module"], pairs_of(len(s) for s in r["modules"]), len(r["acr"]), len(r["acrr"]),
                sum(i in free for i in r["acr"]), sum(a in free and b in free for a, b in r["acrr"])) == \
            (tan["modules"], tan["giant"], tan["pairs"], tan["acr"], tan["acrr_pairs"], tan["acr_free"], tan["acrr_pairs_free"]), (m, r["seed"])
    assert (tan["pairs"], tan["acr"], tan["acr_free"]) == (pairs_of(g["module_sizes"]), g["acr"], g["acr_free"]), m
    out["kimonu_tan"] = dict(n_kinetic=tan["modules"], giant=tan["giant"], acr=tan["acr"], acrr_ratios=ratios(g["acrr_group_sizes"]),
                             acr_free=tan["acr_free"], acrr_ratios_free=ratios(g["acrr_free_group_sizes"]),
                             seconds=statistics.median(r["seconds"] for r in reps), pairs=tan["pairs"],
                             acrr=tan["acrr_pairs"], acrr_free=tan["acrr_pairs_free"])
    return out


def data_s1(path):
    D = json.load(open(f"{B}/figdata.json"))
    PC = json.load(open(f"{B}/figures/pair_counts.json"))
    nC, n = D["n_complexes"], 0
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "network", "complexes", "arm", "seed"] + KEYS + ["coupled_complex_pairs"] + [h for _, h in LISTED])
        for m in sorted(nC, key=lambda m: nC[m]):
            K = kimonu_rows(m, D, PC)
            for arm, label in ARMS:
                v = K.get(arm) or D[arm].get(m)
                if v is None:      # only the reference and COCOA.jl have networks they did not finish
                    w.writerow([m, network(m), nC[m], label, ""] + [NOT_FINISHED if arm.startswith(("ref", "cocoa")) else ""] * len(KEYS) + [""] * (1 + len(LISTED)))
                    n += 1
                    continue
                pc = PC.get(m, {}).get(arm)
                pairs = v["pairs"] if arm in K else pc.get("pairs", "") if isinstance(pc, dict) else ""     # per-seed pair counts are not listed
                for r in (v if isinstance(v, list) else [dict(v, seed="")]):
                    w.writerow([m, network(m), nC[m], label, r.get("seed", "")] + [r.get(k, "") for k in KEYS] + [pairs] + [r.get(k, "") for k, _ in LISTED])
                    n += 1
    return n


def data_s2(path):
    rows, tot = [], dict(attached=0, both=0, lost=0, extracellular=0, periplasmic=0, lost_in_acrr=0, gained=0)
    files = sorted(glob.glob(f"{LISTS}/*.json"))
    assert len(files) == 33, (LISTS, len(files))
    for f in files:
        d = json.load(open(f))
        m = d["model"]
        mets = [x if isinstance(x, str) else x.get("id") for x in json.load(open(f"{B}/unblocked/{m}_fixed.crn.json"))["metabolites"]]
        nfree = mets.index("E_1") if "E_1" in mets else len(mets)
        att = [i for i in d["attached"]["acr"] if i < nfree]
        det = {i for i in d["detached"]["acr"] if i < nfree}
        grp = {i for g in d["detached"]["acrr_groups"] for i in g}
        tot["gained"] += len(det - set(att))
        for i in att:
            comp = mets[i].rsplit("[", 1)[-1].rstrip("]") if "[" in mets[i] else mets[i].rsplit("_", 1)[-1]
            fate = "ACR in both" if i in det else ("in an ACRR group" if i in grp else "free")
            rows.append((m, network(m), mets[i], comp, fate))
            tot["attached"] += 1
            if i in det:
                tot["both"] += 1
            else:
                tot["lost"] += 1
                tot["extracellular"] += comp == "e"
                tot["periplasmic"] += comp == "p"
                tot["lost_in_acrr"] += i in grp
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model", "network", "metabolite", "compartment", "with detached uptakes"])
        w.writerows(rows)
    return len(rows), tot


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    print("Data S1 rows", data_s1(f"{out}/DataS1_values.csv"))
    n, tot = data_s2(f"{out}/DataS2_uptake_acr.csv")
    print("Data S2 rows", n, tot)
