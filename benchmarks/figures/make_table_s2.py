"""Table S2 of the Supplement (the 33 networks) as a Markdown table.

Sizes are those of the networks every tool received (unblocked/<model>_fixed.crn.json). Complexes are counted
without the zero complex (figdata.json). The deficiency counts the zero complex and leaves out complex rows of the
crn.json that occur in no reaction, which only iJB785 and iSynCJ816 have. Linkage classes, rank N, deficiency and dim ker N
come from deficiency/deficiency_33.json (deficiency_split.py), balanced complexes and concordance modules with two or
more complexes from Kimonu's default run (figdata.json, kimonu_det). Rows are sorted by the number of complexes.

    python figures/make_table_s2.py [out.md]        # from benchmarks/; prints the table when no file is given
"""
import json, os, sys

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHORT = {"iIS312_Trypomastigote": "iIS312 Trypomastigote", "iIS312_Amastigote": "iIS312 Amastigote",
         "iIS312_Epimastigote": "iIS312 Epimastigote", "Saha2011iRS1597": "iRS1597",
         "ArabidopsisCoreModel": "Arabidopsis core model", "Chlorella_variabilis_iAJ526": "iAJ526"}


def rows():
    D = json.load(open(f"{B}/figdata.json"))
    nC, det = D["n_complexes"], D["kimonu_det"]
    F = {r["model"]: r for r in json.load(open(f"{B}/deficiency/deficiency_33.json"))}
    out = []
    for m in sorted(nC, key=lambda m: nC[m]):
        crn, f = json.load(open(f"{B}/unblocked/{m}_fixed.crn.json")), F[m]
        # the deficiency counts the zero complex but not complexes that occur in no reaction (iJB785 51, iSynCJ816 29),
        # which the complex column (Kimonu's count) includes
        unused = nC[m] + f["zero_cx"] - f["n"]
        assert unused >= 0, m
        if unused:
            print(f"{m}: {unused} complexes occur in no reaction", file=sys.stderr)
        assert (f["balanced"], f["conc_nonsingle"]) == (det[m]["balanced"], det[m]["conc"]), m
        out.append((SHORT.get(m, m), len(crn["metabolites"]), crn["n_reactions"], nC[m], f["ell"], f["rankN"],
                    f["delta"], f["dimK"], det[m]["balanced"], det[m]["conc"]))
    return out


def markdown(rs):
    md = ["| network | metabolites | reactions | complexes | linkage classes | rank N | deficiency | dim ker N | balanced | concordance modules |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    md += ["| " + " | ".join(str(x) for x in r) + " |" for r in rs]
    return "\n".join(md) + "\n"


if __name__ == "__main__":
    text = markdown(rows())
    if len(sys.argv) > 1:
        open(sys.argv[1], "w").write(text)
        print(f"{sys.argv[1]}: {text.count(chr(10)) - 2} networks")
    else:
        print(text, end="")
