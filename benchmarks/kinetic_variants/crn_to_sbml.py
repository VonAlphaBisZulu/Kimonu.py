"""Write a crn.json network (sweep/ or unblocked/) as SBML, every metabolite in compartment c, for preprocess_time.py.

    python kinetic_variants/crn_to_sbml.py <model>_fixed.crn.json <out.xml>
"""
import json, sys, cobra
from cobra.io import write_sbml_model
cobra.Configuration().solver = "glpk"
crn = json.load(open(sys.argv[1]))
m = cobra.Model("split"); mets = [cobra.Metabolite(x, compartment="c") for x in crn["metabolites"]]; m.add_metabolites(mets)
terms = {}
for i, j, v in crn["S"]: terms.setdefault(int(j) - 1, {})[mets[int(i) - 1]] = float(v)
rx = [cobra.Reaction(crn["reactions"][j], lower_bound=crn["lb"][j], upper_bound=crn["ub"][j]) for j in range(crn["n_reactions"])]
m.add_reactions(rx)
for j, r in enumerate(rx):
    if j in terms: r.add_metabolites(terms[j])
write_sbml_model(m, sys.argv[2]); print("written", len(rx))
