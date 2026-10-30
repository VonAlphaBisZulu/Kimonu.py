"""Detached-uptake data for make_hierarchy.py: the detached block of uptake_modes/<model>[.detached].json in the
schema of hierarchy/<model>.json.

    python kinetic_variants/hierarchy_detached.py   -> kinetic_variants/hierarchy_detached/<model>.json
"""
import glob, json, os
B = os.path.dirname(os.path.abspath(__file__))
os.makedirs(f"{B}/hierarchy_detached", exist_ok=True)
done = set()
for f in sorted(glob.glob(f"{B}/uptake_modes/*.json")):
    d = json.load(open(f))
    if "detached" not in d or d["model"] in done:
        continue
    x = d["detached"]
    out = dict(model=d["model"], n_complexes=d["n_complexes"], n_metabolites=d["n_metabolites"], tangent_dim=x["tangent_dim"],
               seconds_iterated=x["seconds_iterated"], seconds_tangent=x["seconds_tangent"],
               it=x["it"], tan=x["tan"], violations=x["violations"], gap=x["gap"])
    json.dump(out, open(f"{B}/hierarchy_detached/{d['model']}.json", "w"), indent=1)
    done.add(d["model"])
print(len(done), "models")
