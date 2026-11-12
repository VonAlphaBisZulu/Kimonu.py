"""Replace formula text runs of an SVG by LaTeX outlines (Computer Modern), so no math font is needed.

Arial has no glyph for some math signs (proportional to, subset of, superscript T), and Illustrator then
substitutes Myriad Pro or shows boxes. Each formula is set by pdflatex in its own box at the size of the run
it replaces, converted to SVG by pdftocairo, and its glyphs are inlined as plain paths. Words next to a
formula stay Arial text. A formula keeps the left edge of the runs it replaces, or their centre where the
label is centred.

A formula line is (baseline y, x range of the runs replaced, font size, anchor, parts); a part is
("tex", source[, fill]), ("text", words, fill) or ("gap", pt). Needs pdflatex (amsmath, bm) and pdftocairo.
"""
import os, re, subprocess, tempfile

INK = "#1b2430"
REG = {**dict(zip("abcdefghijklmnopqrstuvwxyz", [.556, .556, .5, .556, .556, .278, .556, .556, .222, .222, .5, .222, .833,
                                                  .556, .556, .556, .556, .333, .5, .278, .556, .5, .722, .5, .5, .5])),
       " ": .278, ",": .278, "(": .333, ")": .333, "·": .278, "=": .584, "/": .278, "*": .389, "_": .556, ":": .278,
       "+": .584, "²": .333, "ψ": .713, "∝": .584, "√": .549,
       **dict(zip("ABCKMNPSVWXY", [.667, .667, .722, .667, .833, .722, .667, .667, .667, .944, .667, .667])),
       **dict(zip("0123456789", [.556] * 10))}
arial = lambda s, fs: fs * sum(REG.get(c, .6) for c in s)          # Arial advance width, close enough to centre


def _typeset(sources, sizes):
    """pdflatex one page per formula, cropped to its box; returns (tmp dir, [(width, height, depth)] in bp)"""
    tmp = tempfile.mkdtemp()
    body = "\n".join(
        rf"\setbox0=\hbox{{\fontsize{{{fs}}}{{{fs * 1.2}}}\selectfont${src}$}}"
        rf"\pdfpagewidth=\wd0 \pdfpageheight=\dimexpr\ht0+\dp0\relax"
        rf"\immediate\write16{{BOX {i};\the\wd0;\the\ht0;\the\dp0;}}\shipout\box0"
        for i, (src, fs) in enumerate(zip(sources, sizes)))
    open(f"{tmp}/f.tex", "w").write(r"""\documentclass{article}
\usepackage{amsmath,amssymb,bm}
\pdfhorigin=0pt \pdfvorigin=0pt
\begin{document}
""" + body + "\n\\end{document}\n")
    log = subprocess.run(["pdflatex", "-interaction=nonstopmode", "f.tex"], cwd=tmp, capture_output=True, text=True,
                         env=dict(os.environ, max_print_line="100000")).stdout
    box = {int(i): [float(v) * 72 / 72.27 for v in (w, h, d)]                  # TeX pt -> PostScript pt
           for i, w, h, d in re.findall(r"BOX (\d+);([\d.]+)pt;([\d.]+)pt;([\d.]+)pt;", log)}
    assert len(box) == len(sources), log[-2000:]
    return tmp, [box[i] for i in range(len(sources))]


def _outlines(tmp, page, height, fill):
    """the glyphs and rules of one page as plain paths, with the baseline at y = 0 and the left edge at x = 0"""
    subprocess.run(["pdftocairo", "-svg", "-f", str(page + 1), "-l", str(page + 1), "f.pdf", f"p{page}.svg"], cwd=tmp, check=True)
    g = open(f"{tmp}/p{page}.svg").read()
    sym = {m.group(1): m.group(2) for m in re.finditer(r'<(?:symbol|g)\b[^>]*id="(glyph[^"]+)"[^>]*>\s*<path[^>]*\sd="([^"]*)"', g)}
    o = []
    for m in re.finditer(r'<use xlink:href="#([^"]+)" x="([-\d.e]+)" y="([-\d.e]+)"', g):
        if sym[m.group(1)].strip():
            o.append(f'<path transform="translate({float(m.group(2)):.3f} {float(m.group(3)) - height:.3f})" d="{sym[m.group(1)].strip()}" fill="{fill}"/>')
    for m in re.finditer(r'<path\b([^>]*)/>', g[g.find("</defs>"):]):          # rules: fraction bars, roots, overlines
        d, tr = re.search(r'\sd="([^"]*)"', m.group(1)), re.search(r'transform="([^"]*)"', m.group(1))
        t = f' transform="{tr.group(1)}"' if tr else ""
        o.append(f'<g transform="translate(0 {-height:.3f})"><path{t} d="{d.group(1)}" fill="{fill}"/></g>')
    return o


def _runs(texts, y, x0, x1, fs):
    hit = []
    for m in texts:
        tr = re.search(r'translate\(([-\d.]+)[ ,]+([-\d.]+)\)', m.group(1)); size = re.search(r'font-size="([^"]+)"', m.group(1))
        if tr and abs(float(tr.group(2)) - y) < .2 and x0 <= float(tr.group(1)) < x1 and size and float(size.group(1)) == fs:
            hit.append(m)
    return hit


def replace(s, F, group_id="formulas_latex"):
    """s with the text runs of every formula line in F replaced by outlines; returns (new svg, runs replaced)"""
    parts = [(p, f[2]) for f in F for p in f[4] if p[0] == "tex"]
    tmp, box = _typeset([p[1] for p, _ in parts], [fs for _, fs in parts])
    texts = list(re.finditer(r'<text\b([^>]*)>(.*?)</text>', s, re.S))
    cut, new, page = [], [], 0
    for y, (x0, x1), fs, anchor, line in F:
        hit = _runs(texts, y, x0, x1, fs)
        assert hit, (y, x0, x1)
        xs = [float(re.search(r'translate\(([-\d.]+)', m.group(1)).group(1)) for m in hit]
        last = hit[xs.index(max(xs))]; last_txt = re.sub(r"<[^>]+>", "", last.group(2))
        left, right = min(xs), max(xs) + arial(last_txt, fs)
        cut += hit
        items, x = [], 0.0                                    # lay the parts out from 0, then shift into place
        for part in line:
            if part[0] == "gap":
                x += part[1]; continue
            if part[0] == "text":
                items.append(f'<text transform="translate({x:.3f} 0)" fill="{part[2]}" font-family="ArialMT, Arial" font-size="{fs}">'
                             f'<tspan x="0" y="0" xml:space="preserve">{part[1]}</tspan></text>')
                x += arial(part[1], fs); continue
            fill = part[2] if len(part) > 2 else INK
            items.append(f'<g transform="translate({x:.3f} 0)">' + "".join(_outlines(tmp, page, box[page][1], fill)) + "</g>")
            x += box[page][0]; page += 1
        ox = left if anchor == "left" else (left + right) / 2 - x / 2
        new.append(f'<g transform="translate({ox:.3f} {y})">' + "".join(items) + "</g>")
    for m in sorted(cut, key=lambda m: -m.start()):
        s = s[:m.start()] + s[m.end():]
    s = s.replace("</svg>", f'<g id="{group_id}">\n' + "\n".join(new) + "\n</g>\n</svg>")
    return s, len(cut)


def place(lines, group_id=None):
    """formula lines drawn at given anchors, for generators: each line is (x, y, size, anchor, parts) with anchor
    "left", "center" or "right" at baseline y; returns one <g> per line (all in one group if group_id is set)"""
    parts = [(p, f[2]) for f in lines for p in f[4] if p[0] == "tex"]
    tmp, box = _typeset([p[1] for p, _ in parts], [fs for _, fs in parts]) if parts else (None, [])
    out, page = [], 0
    for x0, y, fs, anchor, line in lines:
        items, x = [], 0.0
        for part in line:
            if part[0] == "gap":
                x += part[1]; continue
            if part[0] == "text":
                items.append(f'<text transform="translate({x:.3f} 0)" fill="{part[2]}" font-family="ArialMT, Arial" font-size="{fs}">'
                             f'<tspan x="0" y="0" xml:space="preserve">{part[1]}</tspan></text>')
                x += arial(part[1], fs); continue
            fill = part[2] if len(part) > 2 else INK
            items.append(f'<g transform="translate({x:.3f} 0)">' + "".join(_outlines(tmp, page, box[page][1], fill)) + "</g>")
            x += box[page][0]; page += 1
        ox = x0 if anchor == "left" else x0 - x / 2 if anchor == "center" else x0 - x
        out.append(f'<g transform="translate({ox:.3f} {y})">' + "".join(items) + "</g>")
    return [f'<g id="{group_id}">', *out, "</g>"] if group_id else out
