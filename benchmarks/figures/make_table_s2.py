"""Tables S2 (the 33 networks) and S3 (kinetic modules, ACR and ACRR per network) of the Supplement as Markdown tables.

Table S2. Sizes are those of the networks every tool received (unblocked/<model>_fixed.crn.json). Complexes are counted
as in the deficiency, the zero complex included and the complexes that occur in no reaction left out (figdata.json
n_complexes, gather_figdata.py). The crn.json lists such complexes only for iJB785 (51) and iSynCJ816 (29). They do not
change the deficiency, since each would add one complex and one linkage class. Linkage classes, rank N, deficiency and
dim ker N come from deficiency/deficiency_33.json (deficiency_split.py), balanced complexes and concordance modules with
two or more complexes from Kimonu's default run (figdata.json, kimonu_det, balanced without the complexes of no
reaction). Rows are sorted by the number of complexes.

Table S3. Iterated route with detached uptakes (kinetic_variants/member_counts/, member_counts.py), in the order of
table S2. A kinetic module has two or more complexes, and the share of complexes in a module is taken of the complexes
of table S2 (Fig. 3C). ACRR is counted as k - 1 independent ratios per group of k metabolites without ACR, and among
free metabolites as k - 1 per group of k free metabolites. The share of free metabolites with ACR or in an ACRR group is
that of Fig. 3D. Gap networks are marked with a dagger. The deficiency is given as a share of dim ker N.

    python figures/make_table_s2.py [s2.md [s3.md]]     # from benchmarks/; prints both tables when no file is given
"""
import json, os, sys

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHORT = {"iIS312_Trypomastigote": "iIS312 Trypomastigote", "iIS312_Amastigote": "iIS312 Amastigote",
         "iIS312_Epimastigote": "iIS312 Epimastigote", "Saha2011iRS1597": "iRS1597",
         "ArabidopsisCoreModel": "Arabidopsis core model", "Chlorella_variabilis_iAJ526": "iAJ526"}
D = json.load(open(f"{B}/figdata.json"))
nC, det = D["n_complexes"], D["kimonu_det"]
F = {r["model"]: r for r in json.load(open(f"{B}/deficiency/deficiency_33.json"))}
MC = {m: json.load(open(f"{B}/kinetic_variants/member_counts/{m}.json")) for m in nC}
ORDER = sorted(nC, key=lambda m: nC[m])


def rows_s2():
    out = []
    for m in ORDER:
        crn, f, gone = json.load(open(f"{B}/unblocked/{m}_fixed.crn.json")), F[m], D["complexes_in_no_reaction"][m]
        # the deficiency's n counts the zero complex and leaves out the complexes of no reaction, as n_complexes does;
        # its balanced count (Kimonu's, from M) still holds them
        assert f["n"] == nC[m] == crn["n_complexes"] - gone and f["delta"] == f["n"] - f["ell"] - f["rankN"], m
        assert (f["balanced"] - gone, f["conc_nonsingle"]) == (det[m]["balanced"], det[m]["conc"]), m
        out.append((SHORT.get(m, m), len(crn["metabolites"]), crn["n_reactions"], nC[m], f["ell"], f["rankN"],
                    f["delta"], f["dimK"], det[m]["balanced"], det[m]["conc"]))
    return out


def rows_s3():
    out, tot = [], [0, 0, 0, 0]
    for m in ORDER:
        d, f = MC[m], F[m]
        it = d["it"]
        assert d["n_complexes"] - D["complexes_in_no_reaction"][m] == nC[m], m     # member_counts counts every row of the crn.json
        acrr = sum(k - 1 for k in it["acrr_group_sizes"])
        acrr_free = sum(k - 1 for k in it["acrr_free_group_sizes"])
        v = [len(it["module_sizes"]), f"{100 * it['complexes_in_modules'] / nC[m]:.1f}", it["acr"], acrr, it["acr_free"], acrr_free,
             f"{100 * it['fixed_members_free'] / d['n_free']:.1f}", f"{100 * f['delta'] / f['dimK']:.1f}"]
        assert min(it["module_sizes"]) >= 2, m
        tot = [a + b for a, b in zip(tot, (it["acr"], acrr, it["acr_free"], acrr_free))]
        out.append((SHORT.get(m, m) + ("†" if d["gap"] else ""), *v))
    out.append(("all 33", "", "", *tot, "", ""))
    return out


def markdown(head, rs):
    md = ["| " + " | ".join(head) + " |", "|---|" + "---:|" * (len(head) - 1)]
    md += ["| " + " | ".join(str(x) for x in r) + " |" for r in rs]
    return "\n".join(md) + "\n"


S2_HEAD = ["network", "metabolites", "reactions", "complexes", "linkage classes", "rank N", "deficiency", "dim ker N",
           "balanced", "concordance modules"]
S3_HEAD = ["network", "kinetic modules", "complexes in a module (%)", "ACR, all", "ACRR, all", "ACR, free", "ACRR, free",
           "free metabolites with ACR or ACRR (%)", "deficiency (% of dim ker N)"]

if __name__ == "__main__":
    s2, s3 = markdown(S2_HEAD, rows_s2()), markdown(S3_HEAD, rows_s3())
    if len(sys.argv) > 1:
        open(sys.argv[1], "w").write(s2)
        print(f"{sys.argv[1]}: {s2.count(chr(10)) - 2} networks")
        if len(sys.argv) > 2:
            open(sys.argv[2], "w").write(s3)
            print(f"{sys.argv[2]}: {s3.count(chr(10)) - 3} networks and the sum")
    else:
        print(s2 + "\n" + s3, end="")
