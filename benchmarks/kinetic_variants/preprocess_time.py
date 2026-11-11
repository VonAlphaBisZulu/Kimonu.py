"""Preprocessing of a network for Kimonu: load the SBML, split reversible reactions, find and remove blocked reactions with the LP, N = Y A.
    python kinetic_variants/preprocess_time.py <model.xml>   Prints the seconds per step (preprocess_time.sbatch)."""
import sys, time
t0 = time.perf_counter()
import cobra
from cobra.io import read_sbml_model
from kimonu.io import build_crn, irreversible
from kimonu.concordance import _preprocess
cobra.Configuration().solver = "glpk"
t = {"import": time.perf_counter() - t0}; a = time.perf_counter()
m = read_sbml_model(sys.argv[1]); t["load"] = time.perf_counter() - a; a = time.perf_counter()
m = irreversible(m).copy(); t["split"] = time.perf_counter() - a; a = time.perf_counter()
blocked, pins = _preprocess(m); t["lp"] = time.perf_counter() - a; a = time.perf_counter()
if blocked: m.remove_reactions(blocked, remove_orphans=True)
t["remove"] = time.perf_counter() - a; a = time.perf_counter()
crn = build_crn(m); t["factorize"] = time.perf_counter() - a
t["total_without_import"] = sum(v for k, v in t.items() if k != "import")
print(sys.argv[1].split("/")[-1], len(m.reactions), "blocked", len(blocked), {k: round(v, 2) for k, v in t.items()}, flush=True)
