"""Data S1 and Data S2 of the Supplement, as CSV files.

Data S1 (DataS1_values.csv) has one row per network, tool and mode, and for the stochastic tools one row per seed.
The values are those of figdata.json (gather_figdata.py). The coupled complex pairs come from figures/pair_counts.json
(pair_sets.py), which has them for Kimonu's default and reference-equivalent mode and for the exhaustive reference on
the 11 networks that it finished. A tool that did not finish within 7 days gets "not finished in 7 days".

Data S2 (DataS2_uptake_acr.csv) lists every free metabolite with ACR under attached uptakes, with its compartment and
what becomes of it under detached uptakes (ACR in both, in an ACRR group, or free). It reads the lists that
uptake_lists.py wrote (iterated route, both uptake modes, metabolites as indices into the network's own list).
Free metabolites are those before "E_1" in that list. The script prints the totals that Supplementary Text S5 quotes.

    python figures/make_data_s1_s2.py <out_dir>        # from benchmarks/

UPTAKE_LISTS replaces /scratch/sadreyev_lab/pschneider/kimonu_uptake_modes/lists (e.g. a rerun).
"""
import csv, glob, json, os, sys

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LISTS = os.environ.get("UPTAKE_LISTS", "/scratch/sadreyev_lab/pschneider/kimonu_uptake_modes/lists")

ARMS = [("kimonu_det", "Kimonu, iterated route, detached uptakes"),
        ("kimonu", "Kimonu, iterated route, attached uptakes"),
        ("kimonu_ref", "Kimonu, reference-equivalent mode"),
        ("ref_exh", "Reference, exhaustive (Langary et al., 2025)"),
        ("ref_sampled", "Reference, sampled (Langary et al., 2025)"),
        ("cocoa_false", "COCOA.jl v1.1.0, efficient=false"),
        ("cocoa_true", "COCOA.jl v1.1.0, efficient=true")]
KEYS = ["balanced", "conc", "n_kinetic", "giant", "acr", "acrr", "acr_free", "acrr_free", "seconds"]
NOT_FINISHED = "not finished in 7 days"


def data_s1(path):
    D = json.load(open(f"{B}/figdata.json"))
    PC = json.load(open(f"{B}/figures/pair_counts.json"))
    nC, n = D["n_complexes"], 0
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "complexes", "arm", "seed"] + KEYS + ["coupled_complex_pairs"])
        for m in sorted(nC, key=lambda m: nC[m]):
            for arm, label in ARMS:
                v = D[arm].get(m)
                if v is None:      # only the reference and COCOA.jl have networks they did not finish
                    w.writerow([m, nC[m], label, ""] + [NOT_FINISHED if arm.startswith(("ref", "cocoa")) else ""] * len(KEYS) + [""])
                    n += 1
                    continue
                pc = PC.get(m, {}).get(arm)
                pairs = pc.get("pairs", "") if isinstance(pc, dict) else ""     # per-seed pair counts are not listed
                for r in (v if isinstance(v, list) else [dict(v, seed="")]):
                    w.writerow([m, nC[m], label, r.get("seed", "")] + [r.get(k, "") for k in KEYS] + [pairs])
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
            fate = "ACR in both" if i in det else ("detached: in an ACRR group" if i in grp else "detached: free")
            rows.append((m, mets[i], comp, fate))
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
        w.writerow(["model", "metabolite", "compartment", "with detached uptakes"])
        w.writerows(rows)
    return len(rows), tot


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    print("Data S1 rows", data_s1(f"{out}/DataS1_values.csv"))
    n, tot = data_s2(f"{out}/DataS2_uptake_acr.csv")
    print("Data S2 rows", n, tot)
