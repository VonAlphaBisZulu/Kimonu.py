"""Preprocessing of a network for Kimonu: load the SBML, then kimonu.preprocess (split reversible reactions, find and
remove blocked reactions with the LP, N = Y A). The step keeps the key "prepare" in the output, as in the records.
    python kinetic_variants/preprocess_time.py <model.xml>   Prints the seconds per step (preprocess_time.sbatch)."""
import sys, time
t0 = time.perf_counter()
import cobra
from cobra.io import read_sbml_model
import kimonu
cobra.Configuration().solver = "glpk"
t = {"import": time.perf_counter() - t0}; a = time.perf_counter()
m = read_sbml_model(sys.argv[1]); t["load"] = time.perf_counter() - a; a = time.perf_counter()
crn = kimonu.preprocess(m); t["prepare"] = time.perf_counter() - a
t["total_without_import"] = sum(v for k, v in t.items() if k != "import")
n_split = sum(2 if r.lower_bound < 0 < r.upper_bound else 1 for r in m.reactions)     # reactions after the split
print(sys.argv[1].split("/")[-1], len(crn.reaction_ids), "blocked", n_split - len(crn.reaction_ids),
      {k: round(v, 2) for k, v in t.items()}, flush=True)
