"""Detached-uptake data for fig3AB_hierarchy.py: the detached block of uptake_modes/<model>[.detached].json in the
schema of hierarchy/<model>.json, with the ACRR counted as independent ratios (k - 1 per group of k metabolites, fields
acrr_ratios and acrr_ratios_free) from the group sizes of member_counts/<model>.json. The pair counts stay as
acrr_pairs and acrr_pairs_free, and the group sizes are checked against them.

    python kinetic_variants/hierarchy_detached.py   -> results/kimonu/hierarchy_detached/<model>.json

Inputs results/kimonu/uptake_modes/ and member_counts/ of KIMONU_RECORDS, output to KIMONU_OUT (records.py).
"""
import glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import records
UM, OUT = records.read(records.UPTAKE_MODES), records.write(records.HIERARCHY_DETACHED, folder=True)
MC = records.read(records.MEMBER_COUNTS)


def with_ratios(rec, members):
    """rec (counts of uptake_modes.py) with the k - 1 counts of the ACRR groups in members (member_counts.py)"""
    g, gf = members["acrr_group_sizes"], members["acrr_free_group_sizes"]
    assert sum(k * (k - 1) // 2 for k in g) == rec["acrr_pairs"] and sum(k * (k - 1) // 2 for k in gf) == rec["acrr_pairs_free"]
    return dict(rec, acrr_ratios=sum(k - 1 for k in g), acrr_ratios_free=sum(k - 1 for k in gf))

os.makedirs(OUT, exist_ok=True)
done = set()
for f in sorted(glob.glob(f"{UM}/*.json")):
    d = json.load(open(f))
    if "detached" not in d or d["model"] in done:
        continue
    x = d["detached"]
    mc = json.load(open(f"{MC}/{d['model']}.json"))
    assert mc["uptakes"] == "detached" and mc["gap"] == x["gap"], d["model"]
    out = dict(model=d["model"], n_complexes=d["n_complexes"], n_metabolites=d["n_metabolites"], tangent_dim=x["tangent_dim"],
               seconds_iterated=x["seconds_iterated"], seconds_tangent=x["seconds_tangent"],
               it=with_ratios(x["it"], mc["it"]), tan=with_ratios(x["tan"], mc["tan"] or mc["it"]),
               violations=x["violations"], gap=x["gap"])
    json.dump(out, open(f"{OUT}/{d['model']}.json", "w"), indent=1)
    done.add(d["model"])
print(len(done), "models")
