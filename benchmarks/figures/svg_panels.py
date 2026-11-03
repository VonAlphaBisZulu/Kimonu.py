"""Restructure a figure SVG for editing in Illustrator: one group per panel, sizes on a fixed scale, formulas as text.

1. Flatten. Every group is dissolved and its transform and inherited presentation attributes move to its
   children, so each drawn element (path, shape, text) sits directly in its panel group. A group with an
   opacity below one keeps its own group, since its children are composited together.
2. Panels. Each element goes to the panel whose rectangle holds the centre of its bounding box; panels are
   given as rectangles in pt, checked in order. One group per panel, id "panel_A", in the original stacking
   order within each panel.
3. Sizes. Font sizes are mapped to a fixed scale (whole points, half points only for scripts). A label keeps
   the edge it is aligned on: labels that share their right end with others stay right-aligned, labels whose
   centres line up stay centred, and text that continues a label on the same line (an exponent, "/R2" after
   "R1") moves with that label's end.
4. Formulas. Outlined formulas (paths of glyphs from pdflatex) inside a given box are replaced by live text
   typeset by matplotlib in Latin Modern, the same fonts as Figure 3, at the same left edge or centre and
   baseline.

    from svg_panels import Figure
"""
import math, re
import xml.etree.ElementTree as ET
import numpy as np
from fontTools.ttLib import TTFont

import svg_finish

SVG = svg_finish.SVG
NS = f"{{{SVG}}}"
DRAWN = {"path", "polygon", "polyline", "line", "rect", "circle", "ellipse", "text", "use", "image"}
INHERITED = ("fill", "fill-opacity", "fill-rule", "stroke", "stroke-width", "stroke-linecap", "stroke-linejoin",
             "stroke-miterlimit", "stroke-dasharray", "stroke-dashoffset", "stroke-opacity", "font-family", "font-size",
             "font-weight", "font-style", "text-anchor", "letter-spacing", "word-spacing", "visibility", "color")


# ---------------------------------------------------------------- geometry
def parse_transform(t):
    """3x3 matrix of an SVG transform list."""
    M = np.eye(3)
    for name, args in re.findall(r"(\w+)\s*\(([^)]*)\)", t or ""):
        a = [float(v) for v in re.split(r"[\s,]+", args.strip()) if v]
        if name == "matrix":
            m = np.array([[a[0], a[2], a[4]], [a[1], a[3], a[5]], [0, 0, 1]])
        elif name == "translate":
            m = np.array([[1, 0, a[0]], [0, 1, a[1] if len(a) > 1 else 0], [0, 0, 1]])
        elif name == "scale":
            sy = a[1] if len(a) > 1 else a[0]
            m = np.array([[a[0], 0, 0], [0, sy, 0], [0, 0, 1]])
        elif name == "rotate":
            r = math.radians(a[0]); c, s = math.cos(r), math.sin(r)
            m = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
            if len(a) == 3:
                m = np.array([[1, 0, a[1]], [0, 1, a[2]], [0, 0, 1]]) @ m @ np.array([[1, 0, -a[1]], [0, 1, -a[2]], [0, 0, 1]])
        elif name == "skewX":
            m = np.array([[1, math.tan(math.radians(a[0])), 0], [0, 1, 0], [0, 0, 1]])
        elif name == "skewY":
            m = np.array([[1, 0, 0], [math.tan(math.radians(a[0])), 1, 0], [0, 0, 1]])
        else:
            raise ValueError(name)
        M = M @ m
    return M


def matrix_str(M):
    a, c, e, b, d, f = M[0, 0], M[0, 1], M[0, 2], M[1, 0], M[1, 1], M[1, 2]
    if abs(a - 1) < 1e-9 and abs(d - 1) < 1e-9 and abs(b) < 1e-9 and abs(c) < 1e-9:
        return f"translate({e:.3f} {f:.3f})"
    return f"matrix({a:.6g} {b:.6g} {c:.6g} {d:.6g} {e:.3f} {f:.3f})"


_NUM = re.compile(r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?")


def path_points(d):
    """End and control points of a path, absolute (a superset of its outline's hull)."""
    pts, cur, start = [], np.zeros(2), np.zeros(2)
    toks = re.findall(r"[MmZzLlHhVvCcSsQqTtAa]|" + _NUM.pattern, d)
    i, cmd = 0, None
    def num():
        nonlocal i
        v = float(toks[i]); i += 1; return v
    while i < len(toks):
        if re.fullmatch(r"[A-Za-z]", toks[i]):
            cmd = toks[i]; i += 1
            if cmd in "Zz":
                cur = start.copy(); continue
        rel = cmd.islower(); C = cmd.upper()
        base = cur if rel else np.zeros(2)
        if C in "ML":
            p = base + [num(), num()]; pts.append(p); cur = p
            if C == "M":
                start = p.copy(); cmd = "l" if rel else "L"
        elif C == "H":
            cur = np.array([(cur[0] if rel else 0) + num(), cur[1]]); pts.append(cur)
        elif C == "V":
            cur = np.array([cur[0], (cur[1] if rel else 0) + num()]); pts.append(cur)
        elif C in "CSQT":
            n = {"C": 3, "S": 2, "Q": 2, "T": 1}[C]
            ps = [base + [num(), num()] for _ in range(n)]
            pts.extend(ps); cur = ps[-1]
        elif C == "A":
            rx, ry = num(), num(); num(); num(); num()
            p = base + [num(), num()]
            pts.extend([p + [rx, ry], p - [rx, ry]]); cur = p
    return np.array(pts) if pts else np.zeros((0, 2))


_ARIAL = {}


def text_width(s, face, size):
    """Advance width of s in pt, from the font's metrics (kerning ignored)."""
    if face not in _ARIAL:
        f = TTFont(svg_finish.FACES[face][0]); _ARIAL[face] = (f.getBestCmap(), f["hmtx"].metrics, f["head"].unitsPerEm)
    cmap, hmtx, upm = _ARIAL[face]
    return sum(hmtx[cmap[ord(c)]][0] if ord(c) in cmap else upm * 0.5 for c in s) * size / upm


def face_of(el, inherited):
    fam = el.get("font-family") or inherited.get("font-family") or "ArialMT, Arial"
    w = el.get("font-weight") or inherited.get("font-weight")
    st = el.get("font-style") or inherited.get("font-style")
    return svg_finish._face(fam, w in ("bold", "700"), st == "italic")


# ---------------------------------------------------------------- the figure
class Leaf:
    """A drawn element with its full transform (M) and the attributes it inherits."""
    def __init__(self, el, M, order):
        self.el, self.M, self.order = el, M, order

    def tag(self):
        return self.el.tag.replace(NS, "")

    def bbox(self):
        """(x0, y0, x1, y1) in figure pt, approximate for text."""
        el, t = self.el, self.tag()
        if t == "path":
            P = path_points(el.get("d", ""))
        elif t in ("polygon", "polyline"):
            v = [float(x) for x in _NUM.findall(el.get("points", ""))]
            P = np.array(v).reshape(-1, 2)
        elif t == "line":
            P = np.array([[float(el.get(k, 0)) for k in ("x1", "y1")], [float(el.get(k, 0)) for k in ("x2", "y2")]])
        elif t in ("rect", "image", "use"):
            x, y = float(el.get("x", 0)), float(el.get("y", 0))
            w, h = float(el.get("width", 0) or 0), float(el.get("height", 0) or 0)
            P = np.array([[x, y], [x + w, y + h]])
        elif t in ("circle", "ellipse"):
            cx, cy = float(el.get("cx", 0)), float(el.get("cy", 0))
            rx = float(el.get("r", el.get("rx", 0))); ry = float(el.get("r", el.get("ry", 0)))
            P = np.array([[cx - rx, cy - ry], [cx + rx, cy + ry]])
        elif t == "text":
            P = self.text_extent()
        else:
            P = np.zeros((0, 2))
        if not len(P):
            return None
        Q = (self.M @ np.c_[P, np.ones(len(P))].T).T[:, :2]
        return Q[:, 0].min(), Q[:, 1].min(), Q[:, 0].max(), Q[:, 1].max()

    def runs(self):
        """[(element owning the characters, xs or None, y, string, face, size)] in reading order, in the text's
        own coordinates; xs is None where the run continues from the previous one (no explicit x)."""
        out = []
        def walk(e, inh, y):
            inh = dict(inh)
            for k in ("font-family", "font-weight", "font-style", "font-size"):
                if e.get(k): inh[k] = e.get(k)
            xs = [float(v) for v in _NUM.findall(e.get("x", ""))] or None
            ys = [float(v) for v in _NUM.findall(e.get("y", ""))]
            y = ys[0] if ys else y
            face, size = face_of(e, inh), float(inh.get("font-size", "8").replace("px", ""))
            if e.text and e.text.strip("\n"):
                out.append((e, xs, y, e.text, face, size))
            for c in e:
                walk(c, inh, y)
                if c.tail and c.tail.strip("\n"):
                    out.append((e, None, y, c.tail, face, size))
        walk(self.el, self.inherited, 0.0)
        return out

    def anchor(self):
        return self.el.get("text-anchor") or _style(self.el).get("text-anchor", "start").strip()

    def text_extent(self):
        P, cur = [], 0.0
        for e, xs, y, s, face, size in self.runs():
            if xs and len(xs) > 1:                            # one x per character (matplotlib)
                for x, ch in zip(xs, s):
                    P += [[x, y - 0.72 * size], [x + text_width(ch, face, size), y + 0.2 * size]]
                cur = xs[-1] + text_width(s[-1], face, size)
            else:
                x0 = xs[0] if xs else cur
                cur = x0 + text_width(s, face, size)
                P += [[x0, y - 0.72 * size], [cur, y + 0.2 * size]]
        P = np.array(P)
        if len(P) and self.anchor() in ("end", "middle"):     # the runs are laid out from the anchor
            w = P[:, 0].max() - P[:, 0].min()
            P[:, 0] -= w if self.anchor() == "end" else w / 2
        return P


_FACE_FILE = {"latinmodern-math": "LatinModernMath-Regular", "lmroman10-regular": "LMRoman10-Regular",
              "lmroman10-bold": "LMRoman10-Bold", "lmroman10-italic": "LMRoman10-Italic",
              "lmroman10-bolditalic": "LMRoman10-BoldItalic"}
_RC = {}


def typeset(tex, size, colour):
    """A formula as live text: (<text> at the baseline origin, ink box (x0, y0, x1, y1) with y up), set by
    matplotlib in Latin Modern exactly as in Figure 3. Script sizes are rounded up to half points."""
    import matplotlib
    from matplotlib.font_manager import FontProperties
    from matplotlib.mathtext import MathTextParser
    from matplotlib.textpath import TextPath
    if not _RC: _RC.update(svg_finish.mathtext_rc())
    with matplotlib.rc_context(_RC):
        prop = FontProperties(size=size)
        _, _, _, glyphs, rects = MathTextParser("path").parse(f"${tex}$", 72, prop)
        V = TextPath((0, 0), f"${tex}$", prop=prop).vertices
    assert not rects, f"{tex}: rules (fraction bars, roots) are not handled"
    runs = []
    for font, fs, num, ox, oy in glyphs:
        name = __import__("os").path.basename(font.fname).rsplit(".", 1)[0]
        ch = chr(num)
        if name.lower().startswith("stix"):                    # lowercase Greek: Latin Modern Math italic
            face = "LatinModernMath-Regular"; ch = svg_finish._greek_math_italic(ch, "Bol" in name)
        else:
            face = _FACE_FILE[name]
        fs = math.ceil(round(fs, 3) * 2) / 2
        if runs and runs[-1][0] == face and runs[-1][1] == fs and abs(runs[-1][2] + oy) < 1e-6:
            runs[-1][3].append(ox); runs[-1][4].append(ch)
        else:
            runs.append([face, fs, -oy, [ox], [ch]])
    text = ET.Element(f"{NS}text", {"fill": colour})
    for face, fs, y, xs, chars in runs:
        _, fam, w, st = svg_finish.FACES[face]
        a = {"x": " ".join(f"{x:.3f}" for x in xs), "y": f"{y:.3f}", "font-family": fam, "font-size": f"{fs:g}"}
        if w: a["font-weight"] = w
        if st: a["font-style"] = st
        ET.SubElement(text, f"{NS}tspan", a).text = "".join(chars)
    return text, (V[:, 0].min(), V[:, 1].min(), V[:, 0].max(), V[:, 1].max())


def _style(el):
    return dict(kv.split(":", 1) for kv in (p.strip() for p in (el.get("style") or "").split(";")) if ":" in kv)


class Figure:
    def __init__(self, svg_text):
        self.root = ET.fromstring(svg_text)
        self.defs = [c for c in self.root if c.tag.replace(NS, "") in ("defs", "style", "metadata", "title")]
        self.leaves, self.nested_defs = [], []
        self._flatten()

    def _flatten(self):
        """Collect the drawn elements with composed transforms; dissolve every group."""
        order = [0]
        def walk(el, M, inh):
            for c in list(el):
                t = c.tag.replace(NS, "")
                if t == "defs":                               # definitions inside a group (matplotlib's markers)
                    self.nested_defs.extend(list(c)); continue
                if t in ("style", "metadata", "title", "clipPath", "mask", "symbol", "desc"):
                    continue
                Mc = M @ parse_transform(c.get("transform"))
                st = _style(c)
                own = dict(inh)
                for k in INHERITED:
                    v = c.get(k) or (st.get(k).strip() if st.get(k) else None)
                    if v is not None: own[k] = v
                if t == "g":
                    op = c.get("opacity") or st.get("opacity")
                    if op and float(op) < 1 and len(list(c.iter())) > 2:
                        raise ValueError("group with opacity: keep it as a unit (not handled yet)")
                    walk(c, Mc, own)
                elif t in DRAWN:
                    for k, v in inh.items():                  # inherited attributes the element does not set
                        if c.get(k) is None and k not in st:
                            c.set(k, v)
                    if "transform" in c.attrib: del c.attrib["transform"]
                    lf = Leaf(c, Mc, order[0]); lf.inherited = {k: c.get(k) for k in INHERITED if c.get(k)}
                    order[0] += 1; self.leaves.append(lf)
                else:
                    raise ValueError(f"unexpected element {t}")
        walk(self.root, np.eye(3), {})

    # ------------------------------------------------------------ panels
    def assign(self, panels):
        """panels: [(name, [(x0, y0, x1, y1), ...])]; each leaf goes to the first panel holding its bbox centre."""
        self.panel_of = {}
        missing = []
        for lf in self.leaves:
            b = lf.bbox()
            if b is None:
                missing.append(lf); continue
            cx, cy = (b[0] + b[2]) / 2, (b[1] + b[3]) / 2
            for name, rects in panels:
                if any(x0 <= cx <= x1 and y0 <= cy <= y1 for x0, y0, x1, y1 in rects):
                    self.panel_of[id(lf)] = name; break
            else:
                missing.append(lf)
        self.panels = [n for n, _ in panels]
        # a white background that reaches beyond its panel (matplotlib's figure patch) would cover the panels
        # drawn before it once each panel is one group; it carries nothing, so it goes
        rects = dict(panels)
        def white(lf):
            st = _style(lf.el)
            fill = (lf.el.get("fill") or st.get("fill", "")).strip().lower()
            stroke = (lf.el.get("stroke") or st.get("stroke", "none")).strip().lower()
            return fill in ("#fff", "#ffffff", "white") and stroke in ("none", "")
        drop = []
        for lf in self.leaves:
            n = self.panel_of.get(id(lf))
            if n is None or lf.tag() not in ("path", "rect") or not white(lf): continue
            b = lf.bbox()
            if not any(b[0] >= x0 - 2 and b[1] >= y0 - 2 and b[2] <= x1 + 2 and b[3] <= y1 + 2 for x0, y0, x1, y1 in rects[n]):
                drop.append(lf)
        gone = {id(l) for l in drop}
        self.leaves = [l for l in self.leaves if id(l) not in gone]
        self.dropped_backgrounds = len(drop)
        return missing

    # ------------------------------------------------------------ sizes
    def _frame(self, lf):
        """Origin, reading direction, up direction and the local-to-figure scale of a text."""
        O = (lf.M @ [0, 0, 1])[:2]
        X = (lf.M @ [1, 0, 1])[:2] - O; Y = (lf.M @ [0, 1, 1])[:2] - O
        k = float(np.hypot(*X))
        return O, X / k, -Y / np.hypot(*Y), k

    def resize(self, rule):
        """Map font sizes. rule(leaf, panel, text, size, face) -> new size (or the same). Labels keep their alignment
        edge (see module docstring); returns {id(leaf): shift of its end along its reading direction, in pt}."""
        texts = [lf for lf in self.leaves if lf.tag() == "text"]
        info = {}
        for lf in texts:
            runs = lf.runs()
            if not runs: continue
            O, d, up, k = self._frame(lf)
            P = lf.text_extent()
            s0, s1 = P[:, 0].min(), P[:, 0].max()               # local start and end (anchor applied)
            base = runs[0][5]
            new = {id(r[0]): rule(lf, self.panel_of[id(lf)], "".join(x[3] for x in runs), r[5], r[4]) for r in runs}
            info[id(lf)] = dict(lf=lf, O=O, d=d, up=up, k=k, s0=s0, s1=s1, base=base, runs=runs, new=new,
                                f=new[id(runs[0][0])] / base, start=O + d * s0 * k, end=O + d * s1 * k,
                                anchored=lf.anchor() in ("end", "middle"))
        # units: a text whose start follows another's end on the same line moves with it
        nxt = {}
        for a in info.values():
            if a["anchored"]: continue
            for b in info.values():
                if b is a or b["anchored"] or self.panel_of[id(a["lf"])] != self.panel_of[id(b["lf"])]: continue
                if float(a["d"] @ b["d"]) < 0.999: continue
                gap = (b["start"] - a["end"]) @ a["d"]; off = (b["start"] - a["start"]) @ a["up"]
                if -0.6 < gap < 1.2 and abs(off) < 0.75 * a["base"] and id(b["lf"]) not in nxt.values():
                    nxt[id(a["lf"])] = id(b["lf"]); break
        heads = [i for i in info if i not in nxt.values()]
        units = []
        for h in heads:
            u = [h]
            while u[-1] in nxt: u.append(nxt[u[-1]])
            units.append(u)
        # alignment of each unit, from the units it lines up with
        def span(u):
            a, z = info[u[0]], info[u[-1]]
            return (a["start"] @ a["d"], z["end"] @ a["d"], a)
        spans = {tuple(u): span(u) for u in units}
        shift = {}
        for u in units:
            U0, U1, a = spans[tuple(u)]
            if all(abs(info[i]["f"] - 1) < 1e-9 and all(abs(info[i]["new"][id(r[0])] - r[5]) < 1e-9 for r in info[i]["runs"]) for i in u):
                continue
            if a["anchored"]:
                align = "self"
            else:
                peers = [v for v in units if v is not u and self.panel_of[id(info[v[0]]["lf"])] == self.panel_of[id(a["lf"])]
                         and float(info[v[0]]["d"] @ a["d"]) > 0.999 and abs(info[v[0]]["base"] - a["base"]) < 0.01]
                same = lambda j: any(abs(spans[tuple(v)][j] - (U0, U1)[j]) < 0.3 for v in peers)
                same_c = any(abs((spans[tuple(v)][0] + spans[tuple(v)][1]) / 2 - (U0 + U1) / 2) < 0.3 for v in peers)
                align = "left" if same(0) else "right" if same(1) else "center" if same_c else "center"
            # new layout of the unit along its reading direction
            pos, cur = [], None
            for j, i in enumerate(u):
                c = info[i]; w_new = (c["s1"] - c["s0"]) * c["k"] * c["f"]
                if j == 0:
                    cur = U0
                else:
                    p = info[u[j - 1]]
                    cur += ((c["start"] - p["end"]) @ a["d"]) * p["f"]
                pos.append((cur, w_new)); cur += w_new
            W_new = cur - U0
            off = {"left": 0.0, "right": (U1 - U0) - W_new, "center": ((U1 - U0) - W_new) / 2, "self": 0.0}[align]
            for j, i in enumerate(u):
                c = info[i]
                if align != "self":
                    along = pos[j][0] + off - c["start"] @ a["d"]
                    perp = ((c["start"] - a["start"]) @ a["up"]) * (a["f"] - 1)
                    delta = a["d"] * along + a["up"] * perp
                    c["lf"].M = np.array([[1, 0, delta[0]], [0, 1, delta[1]], [0, 0, 1]]) @ c["lf"].M
                    shift[i] = along + pos[j][1] - (c["s1"] - c["s0"]) * c["k"]
                done = set()
                for e, xs, y, s, face, size in c["runs"]:
                    if id(e) in done: continue
                    done.add(id(e))
                    n = c["new"][id(e)]
                    (e if e.get("font-size") else c["lf"].el).set("font-size", f"{n:g}")   # where the size is set
                    if e is c["lf"].el or c["anchored"]:
                        continue
                    for at in ("dx", "dy"):                       # offsets of a run inside the text scale with it
                        v = e.get(at)
                        if v and _NUM.fullmatch(v.strip()):
                            e.set(at, f"{float(v) * c['f']:.3f}")
                    if e.get("x") and e is not c["runs"][0][0]:
                        e.set("x", " ".join(f"{c['s0'] + (float(v) - c['s0']) * c['f']:.3f}" for v in _NUM.findall(e.get("x"))))
        return shift

    # ------------------------------------------------------------ formulas
    def replace_formulas(self, pieces, shift=None):
        """pieces: [(box, tex, size, colour, anchor, follows)]. The glyph outlines inside box (filled paths
        without stroke) are removed and the formula is set as live text with its ink aligned on the outlines'
        left edge, centre or right edge (anchor) and vertical centre. follows: text of the label the formula
        continues, whose end shift (from resize) the formula takes over."""
        out = []
        for box, tex, size, colour, anchor, follows in pieces:
            x0, y0, x1, y1 = box
            hit = []
            for lf in self.leaves:
                if lf.tag() != "path" or lf.el.get("fill") in (None, "none") or lf.el.get("stroke") not in (None, "none"):
                    continue
                b = lf.bbox()
                if b and b[0] >= x0 - 0.8 and b[2] <= x1 + 0.8 and b[1] >= y0 - 0.8 and b[3] <= y1 + 0.8:
                    hit.append((lf, b))
            assert hit, f"no outlines for {tex} in {box}"
            bx0, by0 = min(b[0] for _, b in hit), min(b[1] for _, b in hit)
            bx1, by1 = max(b[2] for _, b in hit), max(b[3] for _, b in hit)
            el, ink = typeset(tex, size, colour)
            ax = {"left": bx0 - ink[0], "center": (bx0 + bx1) / 2 - (ink[0] + ink[2]) / 2, "right": bx1 - ink[2]}[anchor]
            if follows:
                lead = [lf for lf in self.leaves if lf.tag() == "text" and "".join(r[3] for r in lf.runs()).strip() == follows
                        and abs(lf.bbox()[1] - by0) < 8]
                assert len(lead) == 1, (tex, follows, len(lead))
                ax += (shift or {}).get(id(lead[0]), 0.0)
            base = (by0 + by1) / 2 + (ink[1] + ink[3]) / 2
            panel = self.panel_of[id(hit[0][0])]
            lf = Leaf(el, np.array([[1, 0, ax], [0, 1, base], [0, 0, 1]]), min(h.order for h, _ in hit))
            lf.inherited = {}
            gone = {id(h) for h, _ in hit}
            self.leaves = [l for l in self.leaves if id(l) not in gone] + [lf]
            self.panel_of[id(lf)] = panel
            out.append((tex, len(hit), round(bx0, 1), round(by0, 1), round(bx1, 1), round(by1, 1)))
        return out

    def svg(self):
        """The SVG text: defs first, then one group per panel."""
        root = ET.Element(self.root.tag, {k: v for k, v in self.root.attrib.items()})
        for d in self.defs:
            root.append(d)
        if self.nested_defs:
            extra = ET.SubElement(root, f"{NS}defs", {"id": "definitions"})
            for d in self.nested_defs:
                extra.append(d)
        order = sorted(self.panels, key=lambda n: (n.split("_")[-1][-1], len(n)))    # A, B, ..., a shared legend after its last panel
        groups = {n: ET.SubElement(root, f"{NS}g", {"id": f"panel_{n}"}) for n in order}
        for lf in sorted(self.leaves, key=lambda l: l.order):
            n = self.panel_of.get(id(lf))
            if n is None:
                raise ValueError(f"element without a panel: {ET.tostring(lf.el)[:120]}")
            if not np.allclose(lf.M, np.eye(3)):
                lf.el.set("transform", matrix_str(lf.M))
            groups[n].append(lf.el)
        return ET.tostring(root, encoding="unicode")
