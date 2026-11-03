"""Finish a figure SVG for the manuscript and for Illustrator.

Every figure is 7.2 in (518.4 pt, 183 mm) wide, the full page width of most journals, and gets explicit
width and height in pt, so every viewer shows all figures at the same scale. Text is live text in two
families. Words are Arial. Formulas are LM Roman 10, the OpenType version of LaTeX's Computer Modern, laid out
by matplotlib with Latin Modern Math; the few glyphs that only the math font has (subset and superset signs,
the transpose sign, greater-or-equal, italic psi) are outlines, so no math font is needed (Illustrator did not
accept Latin Modern Math).

What this script does to an SVG made by matplotlib or exported by Illustrator:

1. Font attributes. Sizes are rounded up to half points (matplotlib sets scripts at 70 %). matplotlib writes the CSS shorthand style="font: 700 6.5px 'Arial'", which Illustrator
   does not read; it becomes font-family, font-size, font-weight and font-style attributes in the form
   Illustrator itself writes, PostScript name first ("Arial-BoldMT, Arial"). Helvetica becomes Arial.
   matplotlib takes lowercase Greek from STIX, which LM Roman lacks; those glyphs move to the
   mathematical italic letters of Latin Modern Math (U+1D6FC on), at the same position.
2. Embedded fonts. Every face in use is subset to the characters it shows and embedded as WOFF2 in an
   @font-face rule named after its PostScript name. Browsers and Obsidian use these, so the figure looks
   the same on every machine. Illustrator ignores @font-face and uses the installed font of that name,
   so LM Roman 10 has to be installed there (Arial ships with macOS and Windows).
3. No clipping. Every clip-path is removed (Illustrator warns that clipping is lost in SVG Tiny).
   The figure scripts keep every mark inside its axes, and check() renders the figure with and
   without clipping to confirm that nothing changes.

    python figures/svg_finish.py <in.svg> <out.svg>          # finish an SVG exported elsewhere
"""
import base64, io, os, re, subprocess, sys, tempfile
import xml.etree.ElementTree as ET
from fontTools import subset
from fontTools.ttLib import TTFont

SVG, XLINK = "http://www.w3.org/2000/svg", "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG); ET.register_namespace("xlink", XLINK)
for _p, _u in (("dc", "http://purl.org/dc/elements/1.1/"), ("cc", "http://creativecommons.org/ns#"),
               ("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#")):
    ET.register_namespace(_p, _u)
LM = os.environ.get("LM_FONTS", "/usr/share/texlive/texmf-dist/fonts/opentype/public")
ARIAL = os.environ.get("ARIAL_FONTS", "/home/sadreyev_lab/pschneider/miniconda3/envs/bioinfo/fonts")
WIDTH_PT = 518.4

# PostScript name: (font file, family list as Illustrator writes it, weight, style)
FACES = {
    "ArialMT": (f"{ARIAL}/arial.ttf", "ArialMT, Arial", None, None),
    "Arial-BoldMT": (f"{ARIAL}/arialbd.ttf", "Arial-BoldMT, Arial", "700", None),
    "Arial-ItalicMT": (f"{ARIAL}/ariali.ttf", "Arial-ItalicMT, Arial", None, "italic"),
    "LatinModernMath-Regular": (f"{LM}/lm-math/latinmodern-math.otf", "LatinModernMath-Regular, 'Latin Modern Math'", None, None),
    "LMRoman10-Regular": (f"{LM}/lm/lmroman10-regular.otf", "LMRoman10-Regular, 'Latin Modern Roman'", None, None),
    "LMRoman10-Bold": (f"{LM}/lm/lmroman10-bold.otf", "LMRoman10-Bold, 'Latin Modern Roman'", "700", None),
    "LMRoman10-Italic": (f"{LM}/lm/lmroman10-italic.otf", "LMRoman10-Italic, 'Latin Modern Roman'", None, "italic"),
    "LMRoman10-BoldItalic": (f"{LM}/lm/lmroman10-bolditalic.otf", "LMRoman10-BoldItalic, 'Latin Modern Roman'", "700", "italic"),
}
MPL_FONTS = [FACES[k][0] for k in FACES if k.startswith(("LatinModern", "LMRoman"))]
TTF_CACHE = os.environ.get("KIMONU_FONT_CACHE", "/scratch/sadreyev_lab/pschneider/kimonu_font_cache")


def _as_truetype(otf):
    """A TrueType copy of a CFF-based OpenType font, cached. matplotlib embeds a font in PDF as TrueType (Type 42), and
    the CFF outlines of Latin Modern would make that embedding invalid. Same glyphs, names and advance widths."""
    from fontTools.pens.cu2quPen import Cu2QuPen
    from fontTools.pens.ttGlyphPen import TTGlyphPen
    from fontTools.ttLib import newTable
    dst = os.path.join(TTF_CACHE, os.path.basename(otf).replace(".otf", ".ttf"))
    if os.path.isfile(dst):
        return dst
    os.makedirs(TTF_CACHE, exist_ok=True)
    f = TTFont(otf)
    gs, order = f.getGlyphSet(), f.getGlyphOrder()
    glyf = newTable("glyf"); glyf.glyphOrder = order; glyf.glyphs = {}
    for name in order:
        pen = TTGlyphPen(gs)
        gs[name].draw(Cu2QuPen(pen, max_err=1.0, reverse_direction=True))
        glyf.glyphs[name] = pen.glyph()
    f["glyf"], f["loca"] = glyf, newTable("loca")
    del f["CFF "]
    if "VORG" in f: del f["VORG"]
    f["maxp"] = newTable("maxp"); m = f["maxp"]
    m.tableVersion, m.numGlyphs = 0x00010000, len(order)
    for k in ("maxZones", "maxTwilightPoints", "maxStorage", "maxFunctionDefs", "maxInstructionDefs", "maxStackElements",
              "maxSizeOfInstructions", "maxComponentElements"):
        setattr(m, k, 1 if k == "maxZones" else 0)
    for k in ("maxPoints", "maxContours", "maxCompositePoints", "maxCompositeContours", "maxComponentDepth"):
        setattr(m, k, 0)
    f["head"].glyphDataFormat = 0; f["head"].indexToLocFormat = 0
    f["post"].formatType = 2.0; f["post"].extraNames = []; f["post"].mapping = {}
    f.sfntVersion = "\x00\x01\x00\x00"
    f.save(dst)
    return dst


def _tex_spacing():
    """Space formulas as TeX does. matplotlib measures its em as the width of an italic 'm' (0.88 em in Latin
    Modern) and puts 0.2 of it around every relation and binary operator; TeX uses the font size as em, 5/18 em
    around relations and arrows and 4/18 em around binary operators."""
    import matplotlib._mathtext as M
    if getattr(M.Parser, "_tex_spacing", False):
        return
    def _make_space(self, percentage):
        st = self.get_state()
        return M.Kern(st.fontsize * st.dpi / 72 * percentage)
    symbol = M.Parser.symbol
    def tex_symbol(self, s, loc, toks):
        out = symbol(self, s, loc, toks)
        c = toks["sym"]
        if len(out) == 1 and isinstance(out[0], M.Hlist) and c in self._spaced_symbols:
            chars = [n for n in out[0].children if not isinstance(n, M.Kern)]
            if len(chars) == 1:
                k = (4 if c in self._binary_operators else 5) / 18
                out = [M.Hlist([self._make_space(k), chars[0], self._make_space(k)], do_kern=True)]
        return out
    M.Parser._make_space, M.Parser.symbol, M.Parser._tex_spacing = _make_space, tex_symbol, True


def mathtext_rc():
    """rcParams for matplotlib: Arial words, Latin Modern formulas, text kept as text in SVG and embedded in PDF."""
    from matplotlib import font_manager
    _tex_spacing()
    for p in MPL_FONTS:
        font_manager.fontManager.addfont(_as_truetype(p))
    return {"font.family": "Arial", "svg.fonttype": "none", "pdf.fonttype": 42, "mathtext.fontset": "custom",
            "mathtext.rm": "Latin Modern Math", "mathtext.it": "Latin Modern Roman:italic",
            "mathtext.bf": "Latin Modern Roman:bold", "mathtext.bfit": "Latin Modern Roman:bold:italic",
            "mathtext.fallback": "stix"}


def _face(family, bold, italic):
    """PostScript name of the face for a CSS or Illustrator family list."""
    f = family.lower()
    if "stix" in f or "latin modern math" in f or "latinmodernmath" in f:
        return "LatinModernMath-Regular"
    if "latin modern roman" in f or "lmroman" in f:
        return "LMRoman10-" + ({(0, 0): "Regular", (1, 0): "Bold", (0, 1): "Italic", (1, 1): "BoldItalic"}[(bold, italic)])
    if "arial" in f or "helvetica" in f or "sans" in f:
        bold = bold or "bold" in f
        italic = italic or "italic" in f or "oblique" in f
        return "Arial-BoldMT" if bold else ("Arial-ItalicMT" if italic else "ArialMT")
    raise ValueError(f"unexpected font family {family!r}")


def _greek_math_italic(ch, bold):
    """Lowercase Greek as a mathematical italic (bold italic) letter of Latin Modern Math."""
    o = ord(ch)
    if 0x3B1 <= o <= 0x3C9:
        return chr((0x1D736 if bold else 0x1D6FC) + o - 0x3B1)
    return ch


def _font_shorthand(v):
    """'italic 700 6.5px 'Arial'' -> (family, size, bold, italic)."""
    m = re.fullmatch(r"\s*(italic\s+|oblique\s+)?(?:(normal|bold|\d{3})\s+)?([\d.]+)px\s+(.+?)\s*", v)
    if not m:
        raise ValueError(f"cannot read font shorthand {v!r}")
    w = m.group(2) or "400"
    return m.group(4).strip("'\""), m.group(3), w in ("bold", "700", "800", "900"), bool(m.group(1))


def _fix_fonts(root):
    """Turn every font setting into Illustrator-style attributes; returns {face: used characters}."""
    used = {}

    def walk(el, inherited, inherited_size):
        face, size = inherited, inherited_size
        style = el.get("style")
        if style and "font:" in style:
            parts = [p.strip() for p in style.split(";") if p.strip()]
            rest = []
            for p in parts:
                k, _, v = p.partition(":")
                if k.strip() == "font":
                    fam, size, bold, italic = _font_shorthand(v)
                    face = _face(fam, bold, italic)
                    stix = "stix" in fam.lower()
                    el.set("font-size", size)
                    if stix:                                        # lowercase Greek from STIX: Latin Modern Math
                        for e in el.iter():
                            if e.text:
                                e.text = "".join(_greek_math_italic(c, bold) for c in e.text)
                else:
                    rest.append(p)
            if rest:
                el.set("style", "; ".join(rest))
            else:
                del el.attrib["style"]
        if el.get("font-family"):
            fam = el.get("font-family")
            bold = el.get("font-weight") in ("bold", "700", "800", "900") or "bold" in fam.lower()
            italic = el.get("font-style") == "italic"
            face = _face(fam, bold, italic)
        if face is not inherited or el.get("font-family"):
            _, families, weight, fstyle = FACES[face]
            el.set("font-family", families)
            for a, v in (("font-weight", weight), ("font-style", fstyle)):
                if v:
                    el.set(a, v)
                elif a in el.attrib:
                    del el.attrib[a]
        if el.get("font-size", "").endswith("px"):
            el.set("font-size", el.get("font-size")[:-2])
        tag = el.tag.rsplit("}", 1)[-1]
        if tag in ("text", "tspan") and el.text and el.text.strip():
            if face is None:
                raise ValueError(f"text without a font: {el.text!r}")
            used.setdefault(face, set()).update(el.text)
        for c in el:
            walk(c, face, size)
            if c.tail and c.tail.strip() and tag in ("text", "tspan"):
                used.setdefault(face, set()).update(c.tail)

    walk(root, None, None)
    return used


def _woff2(face, chars):
    """The face subset to chars, as base64 WOFF2."""
    opts = subset.Options(); opts.flavor = "woff2"; opts.name_IDs = ["*"]; opts.notdef_outline = True
    font = TTFont(FACES[face][0])
    sub = subset.Subsetter(opts); sub.populate(unicodes={ord(c) for c in chars} | {0x20}); sub.subset(font)
    buf = io.BytesIO(); font.flavor = "woff2"; font.save(buf)
    return base64.b64encode(buf.getvalue()).decode()


_GLYPHS = {}


def _outline(ch):
    """SVG path data of a Latin Modern Math glyph in font units (y up), and the units per em."""
    if ch not in _GLYPHS:
        from fontTools.pens.svgPathPen import SVGPathPen
        f = _font("LatinModernMath-Regular"); gs = f.getGlyphSet()
        pen = SVGPathPen(gs); gs[f.getBestCmap()[ord(ch)]].draw(pen)
        _GLYPHS[ch] = (pen.getCommands(), f["head"].unitsPerEm)
    return _GLYPHS[ch]


_FONTS = {}


def _font(face):
    if face not in _FONTS:
        _FONTS[face] = TTFont(FACES[face][0])
    return _FONTS[face]


def _math_as_text_fonts(root):
    """Formulas need no math font. matplotlib lays them out with Latin Modern Math; its letters, digits and brackets
    are the same Computer Modern glyphs as in LM Roman 10, so they are set in that text font, and the few glyphs it
    lacks (relations such as the subset sign, the transpose sign, italic psi) are drawn as outlines at the same place.
    Illustrator did not accept the math font."""
    roman = _font("LMRoman10-Regular").getBestCmap()
    _, fam_roman, _, _ = FACES["LMRoman10-Regular"]
    parent = {c: p for p in root.iter() for c in p}
    for text in list(root.iter(f"{{{SVG}}}text")):
        T = text.get("transform", "")
        st = dict(kv.split(":", 1) for kv in (text.get("style") or "").split(";") if ":" in kv)
        fill = text.get("fill") or st.get("fill", "").strip() or "#000000"
        after = []
        for ts in list(text):
            if not (ts.get("font-family") or "").startswith("LatinModernMath"):
                continue
            s, size = ts.text or "", float(ts.get("font-size") or text.get("font-size"))
            xs = [float(v) for v in re.findall(r"[-+\d.eE]+", ts.get("x", ""))]
            y = float(ts.get("y", "0"))
            assert len(xs) == len(s), f"formula run without one x per glyph: {s!r}"
            runs, i = [], 0
            while i < len(s):                                   # consecutive glyphs of one kind
                j = i
                kind = ord(s[i]) in roman
                while j < len(s) and (ord(s[j]) in roman) == kind:
                    j += 1
                runs.append((kind, i, j)); i = j
            k = list(text).index(ts); text.remove(ts)
            for kind, i, j in runs:
                if kind:
                    a = {kk: v for kk, v in ts.attrib.items() if kk not in ("x", "font-family")}
                    a.update(x=" ".join(f"{v:.3f}" for v in xs[i:j]), **{"font-family": fam_roman})
                    n = ET.Element(f"{{{SVG}}}tspan", a); n.text = s[i:j]
                    text.insert(k, n); k += 1
                else:
                    for c, x in zip(s[i:j], xs[i:j]):
                        d, upm = _outline(c)
                        after.append(ET.Element(f"{{{SVG}}}path", {
                            "d": d, "fill": ts.get("fill") or fill,
                            "transform": f"{T} translate({x:.3f} {y:.3f}) scale({size / upm:.6g} {-size / upm:.6g})".strip()}))
        if after:
            p = parent[text]; k = list(p).index(text)
            for n, el in enumerate(after):
                p.insert(k + 1 + n, el)
            if not "".join(text.itertext()).strip():
                p.remove(text)


def _plain_runs(root):
    """Runs of words as Illustrator writes them: one position per run, ordinary spaces. matplotlib sets the words
    inside a formula string with one x per character and no-break spaces, which Illustrator reported as missing
    alternate glyphs of Arial."""
    for text in list(root.iter(f"{{{SVG}}}text")):
        for el in list(text):
            if el.tag != f"{{{SVG}}}tspan" or not el.text or "Arial" not in (el.get("font-family") or text.get("font-family") or ""):
                continue
            s = el.text
            xs = [float(v) for v in re.findall(r"[-+\d.eE]+", el.get("x", ""))]
            if len(xs) <= 1:
                el.text = s.replace("\xa0", " ")
                continue
            assert len(xs) == len(s), f"run without one x per character: {s!r}"
            # matplotlib puts every character of one style in one run, across the formulas between them: split
            # where the next character does not follow on from this one
            bold = "Bold" in (el.get("font-family") or text.get("font-family") or "")
            f = _font("Arial-BoldMT" if bold else "ArialMT"); cmap, hm, upm = f.getBestCmap(), f["hmtx"].metrics, f["head"].unitsPerEm
            size = float(el.get("font-size") or text.get("font-size"))
            adv = lambda c: hm[cmap.get(ord(c), cmap[32])][0] * size / upm
            segs, a = [], 0
            for i in range(1, len(s) + 1):
                if i == len(s) or abs(xs[i] - xs[i - 1] - adv(s[i - 1])) > 0.25 * size:
                    segs.append((a, i)); a = i
            k = list(text).index(el); text.remove(el)
            for a, b in segs:
                seg = s[a:b]
                i = next((n for n, c in enumerate(seg) if c not in " \xa0"), None)
                if i is None:
                    continue                                    # only spaces: the formula positions carry the gap
                n = ET.Element(f"{{{SVG}}}tspan", {kk: v for kk, v in el.attrib.items() if kk != "x"})
                n.set("x", f"{xs[a + i]:.3f}"); n.text = seg[i:].replace("\xa0", " ").rstrip()
                text.insert(k, n); k += 1


def _missing(face, chars):
    cmap = TTFont(FACES[face][0]).getBestCmap()
    return sorted(c for c in chars if ord(c) not in cmap and c.strip())


def finish(svg_text):
    """The finished SVG text (see module docstring)."""
    root = ET.fromstring(svg_text)
    vb = [float(v) for v in root.get("viewBox").replace(",", " ").split()]
    assert abs(vb[2] - WIDTH_PT) < 0.05, f"figure is {vb[2]} pt wide, not {WIDTH_PT}"
    root.set("width", f"{vb[2]:g}pt"); root.set("height", f"{vb[3]:g}pt")
    # 3. no clipping
    for parent in root.iter():
        for c in list(parent):
            if c.tag == f"{{{SVG}}}clipPath":
                parent.remove(c)
    for el in root.iter():
        if "clip-path" in el.attrib:
            del el.attrib["clip-path"]
        if el.get("style") and "clip-path" in el.get("style"):
            el.set("style", re.sub(r"clip-path:\s*url\([^)]*\);?", "", el.get("style")).strip())
    # 1. fonts as attributes, 2. embedded
    _fix_fonts(root)
    for el in root.iter():                                     # sizes in half points (scripts are 70 % of their base)
        fs = el.get("font-size")
        if fs and abs(float(fs) * 2 - round(float(fs) * 2)) > 1e-6:
            el.set("font-size", f"{__import__('math').ceil(float(fs) * 2) / 2:g}")
    _math_as_text_fonts(root)
    _plain_runs(root)
    used = _fix_fonts(root)
    for face, chars in used.items():
        miss = _missing(face, chars)
        assert not miss, f"{face} has no glyph for {miss}"
    css = "".join(f'@font-face{{font-family:"{face}";src:url(data:font/woff2;base64,{_woff2(face, chars)}) format("woff2");'
                  f'font-weight:{FACES[face][2] or 400};font-style:{FACES[face][3] or "normal"}}}'
                  for face, chars in sorted(used.items()))
    for old in [e for e in root.iter(f"{{{SVG}}}style") if "@font-face" in (e.text or "")]:
        for p in root.iter():
            if old in list(p):
                p.remove(old)
    defs = root.find(f"{{{SVG}}}defs")
    if defs is None:
        defs = ET.Element(f"{{{SVG}}}defs"); root.insert(0, defs)
    st = ET.SubElement(defs, f"{{{SVG}}}style", {"type": "text/css"}); st.text = css
    defs.remove(st); defs.insert(0, st)
    out = ET.tostring(root, encoding="unicode")
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + out + "\n", {f: "".join(sorted(c)) for f, c in used.items()}


def _render(svg_text, png, zoom=3):
    with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False, encoding="utf-8") as f:
        f.write(svg_text)
    subprocess.run(["rsvg-convert", "-z", str(zoom), "-b", "white", f.name, "-o", png], check=True)
    os.unlink(f.name)


def check(before, after, tmpdir, tol=0):
    """Pixels that change when the clipping is removed: render both with clipping as the only difference."""
    import numpy as np
    from PIL import Image
    root = ET.fromstring(before)
    for el in root.iter():
        el.attrib.pop("clip-path", None)
    noclip = ET.tostring(root, encoding="unicode")
    a, b = os.path.join(tmpdir, "clip.png"), os.path.join(tmpdir, "noclip.png")
    _render(before, a); _render(noclip, b)
    A, B = (np.asarray(Image.open(p).convert("L"), dtype=int) for p in (a, b))
    diff = np.abs(A - B) > 40
    ys, xs = np.nonzero(diff)
    k = A.shape[1] / float(root.get("viewBox").replace(",", " ").split()[2])     # pixels per pt
    return int(diff.sum()), (None if not len(xs) else tuple(round(v / k, 1) for v in (xs.min(), ys.min(), xs.max(), ys.max())))


if __name__ == "__main__":
    src, dst = sys.argv[1:3]
    raw = open(src, encoding="utf-8").read()
    n, box = check(raw, raw, tempfile.mkdtemp())
    assert n == 0, f"removing the clipping changes {n} pixels in {box} (pt)"
    out, used = finish(raw)
    open(dst, "w", encoding="utf-8").write(out)
    for f, c in used.items():
        print(f"{f}: {len(c)} characters")
