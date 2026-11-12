"""The example of the README: every step of Kimonu on e_coli_core, which ships with COBRApy as "textbook".

    python examples/quickstart.py

The code below the docstring is the code block of README.md, and tests/test_quickstart.py checks that both agree and
that the script prints the output the README shows.
"""
import cobra, kimonu
model = cobra.io.load_model("textbook")     # e_coli_core, or cobra.io.read_sbml_model("network.xml")
crn = kimonu.preprocess(model)       # reversible reactions split, blocked ones removed, exact N = Y·A
conc = kimonu.concordance(crn)       # K, M, balanced complexes, concordance modules and activity ratios
it = kimonu.iterated_route(crn)      # proven kinetic modules, ACR and ACRR
tan = kimonu.tangent_route(crn)      # every true relation, false couplings possible
cert = kimonu.certify(crn)           # both routes at one r*: certified, proven lists, open relations
k = kimonu.rate_constants(crn, it.flux)     # rate constants at which x = 1 is a steady state, or your own
vals = kimonu.relation_values(crn, k)       # exact values of the relations at these rate constants

modules = [m for m in cert.modules if len(m) > 1]                           # a kinetic module has two or more complexes
gap = len(cert.open_pairs) + len(cert.open_acr) + len(cert.open_acrr)       # relations only the tangent route reports
share = sum(len(m) for m in modules) / len(crn.complexes)                   # the zero complex counts as a complex
print(f"certified {cert.certified}, {len(modules)} kinetic modules, gap {gap}, "
      f"{share:.1%} of {len(crn.complexes)} complexes in kinetic modules")
