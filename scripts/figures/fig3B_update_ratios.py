"""Figure 3B on the count of independent ACRR ratios, edited in place on the vault Figure 3.

Figure 3B showed per network the ACR metabolites plus the ACRR pairs among free metabolites. The paper counts ACRR as
independent ratios, k - 1 per group of k free metabolites, as in Fig. 4D and Data S1. This script moves the marks of
panel B to the new counts and leaves everything else as it was. Each dot moves vertically to its new value, the bar of
each network with a gap gets new ends, its "+n" the new difference, and the title says "ratios" instead of "pairs"
(centred as before). The log axis is kept, and the new values are checked to lie inside it.

The positions are calibrated on the old counts (acr_free + acrr_pairs_free). The dot centres must lie on one line
y = a - b log10(max(v, 0.7)) to within 0.05 pt, and the bar ends must lie on it as well. Each network is found by the x
of its tick label. With --pairs the old counts are drawn again, and the output must render as the input (the check
of this script).

    python figures/fig3B_update_ratios.py <vault Figure3.svg> <out.svg> [--pairs]     # from scripts/

results/kimonu/hierarchy_detached/ of the records (records.py, written by hierarchy_detached.py) gives the counts.
"""
import glob, json, math, os, re, sys, tempfile
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish, svg_panels

src, dest = sys.argv[1], sys.argv[2]
PAIRS = "--pairs" in sys.argv[3:]
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import records
HIERARCHY_DIR = records.read(records.HIERARCHY_DETACHED)
from names import NAME as SHORT                    # the tick labels of the vault file (figures/names.py)
OLD_TITLE, NEW_TITLE = "ACR metabolites + ACRR pairs (free)", "ACR metabolites + ACRR ratios (free)"
GREY, DARK = "#6b7480", "#3a4250"
NS = {"s": svg_finish.SVG}
ET.register_namespace("", svg_finish.SVG); ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")

D = {d["model"]: d for d in (json.load(open(f)) for f in glob.glob(f"{HIERARCHY_DIR}/*.json"))}
assert len(D) == 33, len(D)
old_v = lambda r: r["acr_free"] + r["acrr_pairs_free"]
new_v = old_v if PAIRS else (lambda r: r["acr_free"] + r["acrr_ratios_free"])
floor = lambda v: max(v, 0.7)                             # zero counts sit on the floor of the log axis, as in fig3AB_hierarchy.py

text = open(src, encoding="utf-8").read()
root = ET.fromstring(text)
panel = next(g for g in root.iter(f"{{{svg_finish.SVG}}}g") if g.get("id") == "panel_B")
num = r"-?\d*\.?\d+(?:e-?\d+)?"


def translate(el):
    m = re.fullmatch(rf"translate\(({num})[ ,]+({num})\)", el.get("transform", "translate(0 0)"))
    return float(m.group(1)), float(m.group(2))


# tick labels: x of each network
label_x, texts = {}, []
for el in panel.iter(f"{{{svg_finish.SVG}}}text"):
    s = "".join(el.itertext()).strip()
    texts.append((el, s))
    m = re.fullmatch(rf"matrix\(\S+ -1 1 \S+ ({num}) ({num})\)", el.get("transform", ""))
    if m and el.get("font-size") == "6":                  # the network names (the y label "count" is 7 pt)
        label_x[s] = float(m.group(1))
name_of = {SHORT.get(m, m): m for m in D}
assert set(label_x) == set(name_of), set(label_x) ^ set(name_of)

# dots (grey filled circles) and bars (dark strokes)
dots, bars = [], []
for el in panel.iter(f"{{{svg_finish.SVG}}}path"):
    if el.get("fill") == GREY and not el.get("stroke"):
        dots.append(el)
    elif el.get("stroke") == DARK:
        bars.append(el)


def nearest(x):
    s = min(label_x, key=lambda k: abs(label_x[k] - x))
    return name_of[s], label_x[s] - x


def dot_geometry(el):
    d = el.get("d")
    assert re.match(r"M[^a-zA-Z]+c[^a-zA-Z]+s", d), d
    t = [float(v) for v in re.findall(num, d)]          # M x,y  c (3 points)  s (2 points) ...
    x, ybot = t[0], t[1]
    r = -(t[7] + t[11])                                   # the first quarter of the circle climbs by its radius
    return x, ybot - r + translate(el)[1], r


cal, at = [], {}
for el in dots:
    x, yc, r = dot_geometry(el)
    m, off = nearest(x)
    assert not D[m]["gap"] and m not in at, m
    at[m] = ("dot", el, x, yc)
    cal.append((math.log10(floor(old_v(D[m]["it"]))), yc, off))
offs = [o for *_, o in cal]
assert max(offs) - min(offs) < 0.1, (min(offs), max(offs))
A = np.array([[1, -l] for l, _, _ in cal]); a, b = np.linalg.lstsq(A, np.array([y for _, y, _ in cal]), rcond=None)[0]
Y = lambda v: a - b * math.log10(floor(v))
res = max(abs(Y(10 ** l) - y) for l, y, _ in cal)
assert res < 0.05, res
print(f"calibration: y = {a:.3f} - {b:.3f} log10(v), largest residual {res:.3f} pt over {len(cal)} dots")

# bars: one vertical line and two caps per gap network
vert = {}
for el in bars:
    d = el.get("d")
    mv = re.fullmatch(rf"M({num}),({num})v({num})", d)
    mh = re.fullmatch(rf"M({num}),({num})h({num})", d)
    ty = translate(el)[1]
    if mv:
        x = float(mv.group(1)); m, _ = nearest(x)
        vert[m] = dict(line=el, x=x, lo=float(mv.group(2)) + ty, hi=float(mv.group(2)) + float(mv.group(3)) + ty, ty=ty, caps=[])
for el in bars:
    mh = re.fullmatch(rf"M({num}),({num})h({num})", el.get("d"))
    if mh:
        xc = float(mh.group(1)) + float(mh.group(3)) / 2
        m = min(vert, key=lambda k: abs(vert[k]["x"] - xc))
        vert[m]["caps"].append(el)
gaps = sorted(m for m in D if D[m]["gap"])
assert sorted(vert) == gaps and all(len(v["caps"]) == 2 for v in vert.values()), (sorted(vert), gaps)
for m, v in vert.items():
    for end, rec in (("lo", D[m]["it"]), ("hi", D[m]["tan"])):
        assert abs(Y(old_v(rec)) - v[end]) < 0.05, (m, end, Y(old_v(rec)), v[end])

# the axis: every new value inside the drawn range (between the floor and the top tick label)
ticks = [(float(s.replace(",", "")), translate(el)[1]) for el, s in texts if re.fullmatch(r"[\d,]+", s) and "translate" in el.get("transform", "")]
top = max(v for v, _ in ticks)
assert all(0.7 <= floor(new_v(D[m][k])) <= top * 3 for m in D for k in ("it", "tan")), "a new value leaves the axis"

# move the marks
changed = {}
for m, (_, el, x, yc) in at.items():
    v0, v1 = old_v(D[m]["it"]), new_v(D[m]["it"])
    assert D[m]["it"] == D[m]["tan"] or new_v(D[m]["it"]) == new_v(D[m]["tan"]), m
    if v1 != v0:
        tx, ty = translate(el)
        el.set("transform", f"translate({tx:.3f} {ty + Y(v1) - Y(v0):.3f})")
        changed[m] = (v0, v1)
plus = [(el, s) for el, s in texts if re.fullmatch(r"\+[\d,]+", s)]
assert len(plus) == len(gaps), plus
for m, v in vert.items():
    lo0, hi0, lo1, hi1 = old_v(D[m]["it"]), old_v(D[m]["tan"]), new_v(D[m]["it"]), new_v(D[m]["tan"])
    assert hi1 > lo1, (m, "the routes agree on the new count; draw a dot instead of a bar")
    ylo, yhi = Y(lo1), Y(hi1)
    v["line"].set("d", f"M{v['x']:.2f},{ylo - v['ty']:.2f}v{yhi - ylo:.2f}")
    lower_first = sorted(v["caps"], key=lambda e: -float(re.match(rf"M{num},({num})", e.get("d")).group(1)))
    for cap, y_new in zip(lower_first, (ylo, yhi)):           # the lower cap marks the iterated route
        mh = re.fullmatch(rf"M({num}),({num})h({num})", cap.get("d"))
        cap.set("d", f"M{mh.group(1)},{y_new - v['ty']:.2f}h{mh.group(3)}")
    el, s = min(plus, key=lambda p: abs(translate(p[0])[0] - v["x"]))
    assert s == f"+{hi0 - lo0:,}", (m, s, hi0 - lo0)
    tx, ty = translate(el)
    el.set("transform", f"translate({tx:.3f} {ty + Y(hi1) - Y(hi0):.3f})")
    tsp = el.find("s:tspan", NS)
    (tsp if tsp is not None else el).text = f"+{hi1 - lo1:,}"
    changed[m] = ((lo0, hi0), (lo1, hi1))

# the title, centred as before
title = [el for el, s in texts if s == OLD_TITLE]
assert len(title) == 1, title
if not PAIRS:
    el = title[0]
    face = svg_panels.face_of(el, {})
    tx, ty = translate(el)
    w0, w1 = (svg_panels.text_width(s, face, float(el.get("font-size"))) for s in (OLD_TITLE, NEW_TITLE))
    el.set("transform", f"translate({tx + (w0 - w1) / 2:.3f} {ty:.3f})")
    tsp = el.find("s:tspan", NS)
    (tsp if tsp is not None else el).text = NEW_TITLE

out, used = svg_finish.finish(ET.tostring(root, encoding="unicode"))
open(dest, "w", encoding="utf-8").write(out)

# outside panel B nothing may change; with --pairs nothing at all
tmp = tempfile.mkdtemp()
pa, pb = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(text, pa); svg_finish._render(out, pb)
P, Q = (np.asarray(Image.open(p).convert("L"), dtype=int) for p in (pa, pb))
k = P.shape[1] / 518.4
diff = np.abs(P - Q) > 40
inside = diff[:int(151.5 * k), int(259 * k):].sum()
outside = diff.sum() - inside
print(f"pixels changed: panel B {inside}, elsewhere {outside} | fonts {sorted(used)}")
assert outside == 0 and (inside == 0 or not PAIRS)
print(json.dumps({m: c for m, c in sorted(changed.items())}))
