"""Figure 4 on the final results, edited in place on the vault Figure 4, which carries Phil's Illustrator edits.

make_figure3.py draws Figure 4 from scratch, but the vault file was edited after that (bold titles, legend layout,
the labels of the reference lines in B, live formulas in E, "generic k"). This script keeps all of it and moves only
what the data decide.

A  every marker to its new position (Kimonu's times of the final runs, every tool's network size as table S2 counts
   it), every fit line, and the slopes in the legend.
B  the markers of the exact nullspace (tangent route, timing arms small and tan of the final runs), the slowest time in
   the legend, and the label of the 1 h line without a colon (its right edge kept).
F  Kimonu's times on iLJ478, the size of iLJ478 (3440 complexes with the zero complex) and the cells without colons.
C, D, their legend and E stay as they are (Kimonu's pairs and relations did not change), which is checked by
rendering both files.

The positions are calibrated on the tick marks of each axis, and the markers are matched to networks through the data
the vault figure was drawn from: figdata.json of commit 451a748 (A) and the tangent runs in OLD_TIMING (B). Every
matched marker must lie within 0.3 pt of its old position, and every old fit line must start where the old data put it.

    python figures/figure4_update.py <vault Figure4.svg> <out.svg>     # from benchmarks/, after gather_figdata.py

KIMONU_TIMING (default the rerun of 2026-11-11) gives the new tangent times, OLD_TIMING (default
/scratch/sadreyev_lab/pschneider/kimonu_timing_seeds) the old ones, OLD_REV (default 451a748) the old figdata.json.
"""
import glob, json, os, re, statistics, subprocess, sys, tempfile
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish, svg_panels

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(B)
src, dest = sys.argv[1], sys.argv[2]
NEW_T = os.environ.get("KIMONU_TIMING", "/scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/timing")
OLD_T = os.environ.get("OLD_TIMING", "/scratch/sadreyev_lab/pschneider/kimonu_timing_seeds")
OLD_REV = os.environ.get("OLD_REV", "451a748")
NEW = json.load(open("figdata.json"))
OLD = json.loads(subprocess.run(["git", "show", f"{OLD_REV}:benchmarks/figdata.json"], capture_output=True, text=True, check=True).stdout)
SVG = svg_finish.SVG
ET.register_namespace("", SVG); ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
GREEN, BLUE, RED, WHITE = "#4daf4a", "#377eb8", "#e41a1c", "#fff"
ARMS = {"kimonu_det": (GREEN, True, False), "kimonu_ref": (GREEN, False, True), "ref_exh": (BLUE, True, False),
        "ref_sampled": (BLUE, False, True), "cocoa_false": (RED, True, False), "cocoa_true": (RED, False, True)}  # colour, filled, dashed
NUM = r"-?(?:\d+\.?\d*|\.\d+)"
med = statistics.median

text = open(src, encoding="utf-8").read()
root = ET.fromstring(text)
panel = {g.get("id"): g for g in root if g.tag == f"{{{SVG}}}g"}


def col(c):
    c = (c or "").lower()
    return {"#ffffff": WHITE, "white": WHITE}.get(c, c)


def offset(el, strict=True):
    """the translation of an element (strict: nothing but a translation)"""
    t = el.get("transform")
    if not t:
        return np.zeros(2)
    m = re.fullmatch(rf"\s*translate\(\s*({NUM})[ ,]+({NUM})\s*\)\s*" + ("" if strict else ".*"), t)
    assert m, t
    return np.array([float(m.group(1)), float(m.group(2))])


def center(el):
    p = np.array(svg_panels.path_points(el.get("d")))
    return (p.min(0) + p.max(0)) / 2 + offset(el)


def move(el, dx, dy):
    """translate a path whose first command is an absolute M and all others relative"""
    m = re.match(rf"\s*M\s*({NUM})[ ,]?({NUM})(.*)$", el.get("d"), re.S)
    assert m and not re.search(r"[A-Y]", m.group(3)), el.get("d")[:80]
    el.set("d", f"M{float(m.group(1)) + dx:.3f},{float(m.group(2)) + dy:.3f}{m.group(3)}")


def axis(g, xs_v, ys_v):
    """x(log10 v), y(log10 v) from the tick marks of a panel: vertical ticks below the x axis, horizontal ones left of y"""
    xt, yt = [], []
    for p in g.iter(f"{{{SVG}}}path"):
        d = p.get("d", "")
        if m := re.fullmatch(rf"M({NUM}),({NUM})v2\.5", d): xt.append(float(m.group(1)))
        if m := re.fullmatch(rf"M({NUM}),({NUM})h-2\.5", d): yt.append(float(m.group(2)))
    assert len(xt) == len(xs_v) and len(yt) == len(ys_v), (len(xt), len(yt))
    fx = np.polyfit(np.log10(xs_v), sorted(xt), 1); fy = np.polyfit(np.log10(ys_v), sorted(yt, reverse=True), 1)
    for f, v, t in ((fx, xs_v, sorted(xt)), (fy, ys_v, sorted(yt, reverse=True))):
        assert np.abs(np.polyval(f, np.log10(v)) - t).max() < 0.03, (f, t)
    return lambda x, y: np.array([np.polyval(fx, np.log10(x)), np.polyval(fy, np.log10(y))])


def groups(g):
    """marker paths of one scatter call (one marker id, as Illustrator keeps it in data-name), by base id"""
    out = {}
    for p in g.iter(f"{{{SVG}}}path"):
        i = p.get("data-name") or p.get("id")
        if i and re.fullmatch(r"m[0-9a-f]{10}", i):
            out.setdefault(i, []).append(p)
    return out


def match(paths, old_xy, tol=0.3):
    """paths matched to the old positions (one per network); returns {network: path}"""
    names = list(old_xy)
    P = np.array([old_xy[n] for n in names]); Q = np.array([center(p) for p in paths])
    assert len(P) == len(Q), (len(P), len(Q))
    cost = np.linalg.norm(P[:, None] - Q[None], axis=2)
    r, c = linear_sum_assignment(cost)
    worst = cost[r, c].max()
    assert worst < tol, f"a marker lies {worst:.2f} pt from where the old data put it"
    return {names[i]: paths[j] for i, j in zip(r, c)}, worst


sec = lambda v: med(r["seconds"] for r in v) if isinstance(v, list) else v["seconds"]
report = []

# ------------------------------------------------------------------ A
gA = panel["panel_A"]
XA = axis(gA, [1000, 2000, 5000, 10000, 20000], [0.1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6])
G = groups(gA)
slopes = {}
for arm, (c, filled, dashed) in ARMS.items():
    pick = [ps for ps in G.values() if len(ps) > 1 and col(ps[0].get("fill")) == (c if filled else WHITE)
            and col(ps[0].get("stroke")) == (WHITE if filled else c)]
    assert len(pick) == 1, (arm, len(pick))
    ms = [m for m in OLD["models"] if m in OLD[arm]]
    assert ms == [m for m in NEW["models"] if m in NEW[arm]], arm
    hit, worst = match(pick[0], {m: XA(OLD["n_complexes"][m], sec(OLD[arm][m])) for m in ms})
    for m, p in hit.items():
        dx, dy = XA(NEW["n_complexes"][m], sec(NEW[arm][m])) - XA(OLD["n_complexes"][m], sec(OLD[arm][m]))
        move(p, dx, dy)
    # the fit line, as make_figure3.py draws it
    fits = {}
    for key, D in (("old", OLD), ("new", NEW)):
        X = np.array([D["n_complexes"][m] for m in ms], float); Y = np.array([sec(D[arm][m]) for m in ms], float)
        b, a0 = np.polyfit(np.log10(X), np.log10(Y), 1); gx = np.linspace(np.log10(X.min()), np.log10(X.max()), 50)
        fits[key] = (b, [XA(10 ** x, 10 ** (a0 + b * x)) for x in gx])
    lines = [p for p in gA.iter(f"{{{SVG}}}path") if col(p.get("stroke")) == c and col(p.get("fill")) in ("none", "")
             and len(re.findall(r"[lL]", p.get("d", ""))) >= 1 and len(svg_panels.path_points(p.get("d"))) >= 20
             and (p.get("stroke-dasharray") is not None) == dashed]
    assert len(lines) == 1, (arm, len(lines))
    pts = np.array(svg_panels.path_points(lines[0].get("d"))) + offset(lines[0])
    assert np.abs(pts[0] - fits["old"][1][0]).max() < 0.3 and np.abs(pts[-1] - fits["old"][1][-1]).max() < 0.3, (arm, pts[0], fits["old"][1][0])
    lines[0].set("d", "M" + "L".join(f"{x - offset(lines[0])[0]:.3f},{y - offset(lines[0])[1]:.3f}" for x, y in fits["new"][1]))
    slopes[arm] = (fits["old"][0], fits["new"][0])
    report.append(f"A {arm}: {len(hit)} markers moved (matched within {worst:.2f} pt), slope {fits['old'][0]:.2f} -> {fits['new'][0]:.2f}")
done = set()
for ts in gA.iter(f"{{{SVG}}}tspan"):
    m = re.fullmatch(r"(\d+) models, slope (\d\.\d\d)", ts.text or "")
    if not m: continue
    arm = [a for a, (o, _) in slopes.items() if f"{o:.2f}" == m.group(2) and a not in done]
    assert len(arm) == 1, (ts.text, arm)
    ts.text = f"{m.group(1)} models, slope {slopes[arm[0]][1]:.2f}"; done.add(arm[0])
assert done == set(ARMS), done

# ------------------------------------------------------------------ B
gB = panel["panel_B"]
XB = axis(gB, [10, 100, 1000, 10000], [1e-3, 1e-2, 0.1, 1, 10, 100, 1e3])


def exact(T):
    """{network: (complexes, median seconds)} of the tangent runs drawn in B: arms small and tan"""
    runs = {}
    for arm in ("small", "tan"):
        for f in glob.glob(f"{T}/{arm}/tangent_*_rep*.json"):
            r = json.load(open(f)); runs.setdefault((arm, r["network"]), []).append(r)
    return {k: (v[0]["complexes"], float(np.median([x["seconds"] for x in v]))) for k, v in runs.items()}


rows = {r["network"] for r in json.load(open(f"{os.environ.get('KIMONU_LACR', '/scratch/sadreyev_lab/pschneider/kimonu_local_acr')}/small.json"))}
rows |= {os.path.basename(f)[7:-5] for f in glob.glob(f"{os.environ.get('KIMONU_LACR', '/scratch/sadreyev_lab/pschneider/kimonu_local_acr')}/ladder_*.json")}
eo, en = exact(OLD_T), exact(NEW_T)
keys = [k for k in eo if k[0] == "tan" or k[1] in rows]
assert sorted(keys) == sorted(k for k in en if k[0] == "tan" or k[1] in rows), "the networks of B differ"
pick = [ps for ps in groups(gB).values() if len(ps) == len(keys) and col(ps[0].get("fill")) == GREEN]
assert len(pick) == 1, len(pick)
hit, worst = match(pick[0], {k: XB(eo[k][0], max(eo[k][1], 1e-3)) for k in keys})
for k, p in hit.items():
    assert eo[k][0] == en[k][0], k
    dx, dy = XB(en[k][0], max(en[k][1], 1e-3)) - XB(eo[k][0], max(eo[k][1], 1e-3))
    move(p, dx, dy)
slow_o = max(v[1] for k, v in eo.items() if k[0] == "tan"); slow_n = max(v[1] for k, v in en.items() if k[0] == "tan")
n_gen = sum(1 for k in en if k[0] == "tan")
report.append(f"B: {len(hit)} markers moved (matched within {worst:.2f} pt), slowest of {n_gen} models {slow_o:.0f} s -> {slow_n:.0f} s")


def width(s, el):
    face = "Arial-BoldMT" if "Bold" in (el.get("font-family") or "") else "ArialMT"
    return svg_panels.text_width(s, face, float(el.get("font-size")))


def retext(g, old, new, keep="start", count=1):
    """replace the text of the tspans that read old; keep "end" moves the text so that its right edge stays"""
    hits = [(t, ts) for t in g.iter(f"{{{SVG}}}text") for ts in t.iter(f"{{{SVG}}}tspan") if (ts.text or "") == old]
    assert len(hits) == count, (old, len(hits))
    for t, ts in hits:
        if keep != "start":
            sh = (width(old, t) - width(new, t)) * (1 if keep == "end" else 0.5)
            x, y = offset(t); t.set("transform", f"translate({x + sh:.3f} {y:.3f})")
        ts.text = new


retext(gB, f"slowest of {n_gen} models: {slow_o:.0f} s, one core", f"at most {slow_n:.0f} s on {n_gen} models, one core")
retext(gB, "1 h limit (open: not finished)", "1 h limit (open markers did not finish)", keep="end")

# ------------------------------------------------------------------ F
gF = panel["panel_F"]
o, n = OLD["kimonu_det"]["iLJ478"], NEW["kimonu_det"]["iLJ478"]
for a, b in ((f"{o['sec_conc']:.1f} s", f"{n['sec_conc']:.1f} s"), (f"{o['sec_kin'] + o['sec_acr']:.1f} s", f"{n['sec_kin'] + n['sec_acr']:.1f} s"),
             (f"{OLD['kimonu_ref']['iLJ478']['seconds']:.1f} s", f"{NEW['kimonu_ref']['iLJ478']['seconds']:.1f} s")):
    retext(gF, a, b)
    report.append(f"F: {a} -> {b}")
nc = NEW["n_complexes"]["iLJ478"]
retext(gF, f"iLJ478: {OLD['n_complexes']['iLJ478']} complexes, 5.9 M complex pairs", f"iLJ478, {nc} complexes, {nc * (nc - 1) / 2e6:.1f} M complex pairs",
       keep="middle")
for a, b, k in (("skipped: 2.37 M balanced", "2.37 M balanced pairs skipped", 1), ("ACR, ACRR: Shinar–Feinberg rule", "ACR, ACRR by Shinar–Feinberg rule", 2),
                ("skipped: 5.82 M by three filters", "5.82 M skipped by three filters", 1), ("1 nullspace:", "1 nullspace,", 1),
                ("skipped: all, one readout", "no pair tested, one readout", 1), ("ACR, ACRR: zero and equal rows", "ACR, ACRR from zero and equal rows", 1)):
    retext(gF, a, b, count=k)
# the cells' texts must stay inside their boxes
boxes = [np.array(svg_panels.path_points(p.get("d"))) for p in gF.iter(f"{{{SVG}}}path") if len(svg_panels.path_points(p.get("d") or "M0,0")) >= 8
         and col(p.get("fill")) not in ("none", "")]
for t in gF.iter(f"{{{SVG}}}text"):
    s = "".join(ts.text or "" for ts in t.iter(f"{{{SVG}}}tspan"))
    if not s.strip() or "Arial" not in (t.get("font-family") or ""): continue
    x0, y0 = offset(t, False); x1 = x0 + width(s, t)
    inside = [b for b in boxes if b[:, 0].min() - 0.5 <= x0 <= b[:, 0].max() and b[:, 1].min() <= y0 <= b[:, 1].max()]
    for b in inside:
        assert x1 <= b[:, 0].max() - 1.0, f"{s!r} ends at {x1:.1f} pt, its box at {b[:, 0].max():.1f} pt"

svg, used = svg_finish.finish(ET.tostring(root, encoding="unicode"))
open(dest, "w", encoding="utf-8").write(svg)

# ------------------------------------------------------------------ check: what changed, per panel
tmp = tempfile.mkdtemp()
a, b = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(text, a); svg_finish._render(svg, b)
A = np.asarray(Image.open(a).convert("L"), dtype=int); Bm = np.asarray(Image.open(b).convert("L"), dtype=int)
assert A.shape == Bm.shape
k = A.shape[1] / 518.4
ch = np.abs(A - Bm) > 40
for name, g in panel.items():
    pts = [np.array(svg_panels.path_points(p.get("d"))) + offset(p, False) for p in g.iter(f"{{{SVG}}}path") if p.get("d")]
    pts += [offset(t, False)[None] for t in g.iter(f"{{{SVG}}}text")]
    P = np.vstack(pts); x0, y0 = P.min(0); x1, y1 = P.max(0)
    report.append(f"{name}: {int(ch[max(int(y0 * k), 0):int(y1 * k) + 1, max(int(x0 * k), 0):int(x1 * k) + 1].sum())} pixels changed")
print("\n".join(report))
print(f"fonts {sorted(used)}")
for name in ("panel_C", "panel_D", "panel_legend_CD", "panel_E"):
    g = panel[name]
    pts = [np.array(svg_panels.path_points(p.get("d"))) + offset(p, False) for p in g.iter(f"{{{SVG}}}path") if p.get("d")]
    P = np.vstack(pts); x0, y0 = P.min(0); x1, y1 = P.max(0)
    assert ch[int(y0 * k) + 2:int(y1 * k) - 1, int(x0 * k) + 2:int(x1 * k) - 1].sum() == 0, f"{name} changed"
