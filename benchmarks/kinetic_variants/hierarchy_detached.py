"""Detached-uptake data for make_hierarchy.py: the detached block of uptake_modes/<model>[.detached].json in the
schema of hierarchy/<model>.json.

    python kinetic_variants/hierarchy_detached.py   -> kinetic_variants/hierarchy_detached/<model>.json

UM_DIR replaces the input directory kinetic_variants/uptake_modes and HD_OUT the output directory (a rerun).
"""
import glob, json, os
B = os.path.dirname(os.path.abspath(__file__))
UM = os.environ.get("UM_DIR", f"{B}/uptake_modes"); OUT = os.environ.get("HD_OUT", f"{B}/hierarchy_detached")
os.makedirs(OUT, exist_ok=True)
done = set()
for f in sorted(glob.glob(f"{UM}/*.json")):
    d = json.load(open(f))
    if "detached" not in d or d["model"] in done:
        continue
    x = d["detached"]
    out = dict(model=d["model"], n_complexes=d["n_complexes"], n_metabolites=d["n_metabolites"], tangent_dim=x["tangent_dim"],
               seconds_iterated=x["seconds_iterated"], seconds_tangent=x["seconds_tangent"],
               it=x["it"], tan=x["tan"], violations=x["violations"], gap=x["gap"])
    json.dump(out, open(f"{OUT}/{d['model']}.json", "w"), indent=1)
    done.add(d["model"])
print(len(done), "models")
