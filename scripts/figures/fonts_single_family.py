"""Latin Modern Roman 10 with one family name, so that Windows, Chrome and Illustrator all see family "Latin Modern Roman".

The originals carry two family names: "LM Roman 10" (name ID 1) and "Latin Modern Roman" with style "10 Bold" (IDs 16, 17).
Windows lists the first, Illustrator groups by the second, and the SVG names the second, which Chrome on Windows does not
match. Here IDs 1 and 2 become "Latin Modern Roman" and the plain style, 16 and 17 go, and the PostScript names stay
(LMRoman10-Bold and so on), so existing figures and svg_finish keep finding the faces.
    python rename.py <font dir> <out dir>
"""
import os, sys
from fontTools.ttLib import TTFont

src, out = sys.argv[1:3]
FAMILY = "Latin Modern Roman"
for f, style in [("lmroman10-regular.otf", "Regular"), ("lmroman10-bold.otf", "Bold"), ("lmroman10-italic.otf", "Italic"), ("lmroman10-bolditalic.otf", "Bold Italic")]:
    t = TTFont(os.path.join(src, f))
    name = t["name"]
    ps = name.getDebugName(6)
    full = f"{FAMILY} {style}"
    for rec in list(name.names):
        if rec.nameID in (16, 17, 18, 21, 22):
            name.removeNames(nameID=rec.nameID)
    for pid, eid, lid in {(r.platformID, r.platEncID, r.langID) for r in name.names}:
        for nid, val in [(1, FAMILY), (2, style), (3, f"2.004;UKWN;{ps};single-family"), (4, full)]:
            name.setName(val, nid, pid, eid, lid)
    td = t["CFF "].cff.topDictIndex[0]
    td.FamilyName, td.FullName = FAMILY, full
    dst = os.path.join(out, f"LatinModernRoman-{style.replace(' ', '')}.otf")
    t.save(dst)
    print(os.path.basename(dst), "| PostScript", ps)
