"""Figures 3 and 4 in the words and network names of the text, edited in place on the vault files.

Figure 4 carries Phil's Illustrator edits and Figure 3 those of its chain (scripts/README.md), so both are edited
where they are instead of drawn again. Only these texts change.

Figure 3  A, B  the iIS312 stages named as in every other item (figures/names.py), "iIS312 Try" for "iIS312 T" and so
                on, each label's end kept at the axis.
          E, F  "nullspace" for "kernel" and "networks" for "models", in the titles, the y label and the note.
Figure 4  A, B  "networks" for "models" in the legends.

Checks: every old text occurs as often as expected, no new text box runs into another text or out of the page, and
rendering both files (rsvg-convert) changes no pixel outside the boxes of the changed texts. A file that already reads
the new words is reported and copied unchanged.

    python figures/fig3_fig4_update_terms.py <vault figures dir> <out dir>     # from scripts/; writes Figure3.svg, Figure4.svg
"""
import os, re, sys, tempfile
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish, svg_panels
from names import network

SVG = svg_finish.SVG
ET.register_namespace("", SVG); ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
NUM = r"-?(?:\d+\.?\d*(?:e-?\d+)?|\.\d+)"
STAGE = {"iIS312 T": network("iIS312_Trypomastigote"), "iIS312 A": network("iIS312_Amastigote"),
         "iIS312 E": network("iIS312_Epimastigote")}
# figure: [(panel, old, new, keep, count)]; keep "end" holds the end of the text in place, otherwise its anchor stays
EDITS = {"Figure3.svg": [(p, o, n, "end", 1) for p in ("panel_A", "panel_B") for o, n in STAGE.items()] + [
             ("panel_E", "models", "networks", "anchor", 1),
             ("panel_E", "exact kernel", "exact nullspace", "anchor", 1),
             ("panel_F", "SVD kernel", "SVD nullspace", "anchor", 1),
             (None, "30 of 33 models, and the other 3 are", "30 of 33 networks, and the other 3 are", "anchor", 1),
             (None, "repairs this. The SVD kernel", "repairs this. The SVD nullspace", "anchor", 1)],
         "Figure4.svg": [("panel_A", f"{n} models, slope {s}", f"{n} networks, slope {s}", "anchor", 1)
                         for n, s in (("11", "2.57"), ("11", "3.61"), ("30", "3.14"), ("31", "2.87"), ("33", "1.15"), ("33", "1.88"))] + [
             ("panel_B", "at most 439 s on 33 models, one core", "at most 439 s on 33 networks, one core", "anchor", 1)]}


def matrix(t):
    """the 3x3 matrix of a transform attribute (translate, matrix or none)"""
    M = np.eye(3)
    for kind, args in re.findall(r"(translate|matrix)\(([^)]*)\)", t or ""):
        v = [float(x) for x in re.findall(NUM, args)]
        if kind == "translate":
            m = np.array([[1, 0, v[0]], [0, 1, v[1] if len(v) > 1 else 0], [0, 0, 1]])
        else:
            m = np.array([[v[0], v[2], v[4]], [v[1], v[3], v[5]], [0, 0, 1]])
        M = M @ m
    return M


def string(t):
    return "".join(t.itertext())


def anchor(t):
    st = dict(kv.split(":", 1) for kv in (t.get("style") or "").split(";") if ":" in kv)
    return (t.get("text-anchor") or st.get("text-anchor", "start")).strip()


def width(s, t):
    face = "Arial-BoldMT" if "Bold" in (t.get("font-family") or "") else "ArialMT"
    return svg_panels.text_width(s, face, float(t.get("font-size")))


def local_xy(t):
    ts = t.find(f"{{{SVG}}}tspan")
    el = ts if ts is not None and ts.get("x") is not None else t
    return float(el.get("x") or 0), float(el.get("y") or 0)


def box(t, s=None):
    """axis-aligned box of the glyph run of text t (with string s instead of its own), in page pt"""
    s = string(t) if s is None else s
    w, size = width(s, t), float(t.get("font-size"))
    x, y = local_xy(t)
    x0 = x - w * {"start": 0, "middle": 0.5, "end": 1}[anchor(t)]
    corners = np.array([[x0, y - 0.75 * size, 1], [x0 + w, y - 0.75 * size, 1], [x0, y + 0.22 * size, 1], [x0 + w, y + 0.22 * size, 1]])
    P = (matrix(t.get("transform")) @ corners.T).T[:, :2]
    return np.r_[P.min(0), P.max(0)]


def fmt(M, old):
    if old.startswith("translate"):
        return f"translate({M[0, 2]:.3f} {M[1, 2]:.3f})"
    a, b, c, d = (old_v for old_v in re.findall(NUM, old)[:4])          # keep the rotation as written
    return f"matrix({a} {b} {c} {d} {M[0, 2]:.3f} {M[1, 2]:.3f})"


def edit(src, dest, edits):
    text = open(src, encoding="utf-8").read()
    root = ET.fromstring(text)
    panels = {g.get("id"): g for g in root if g.tag == f"{{{SVG}}}g"}
    texts = [t for t in root.iter(f"{{{SVG}}}text")]
    if all(sum(string(t) == o for t in texts) == 0 and sum(string(t) == n for t in texts) >= k for _, o, n, _, k in edits):
        open(dest, "w", encoding="utf-8").write(text)
        return [f"{os.path.basename(src)}: reads the new words already, copied"]
    report, changed, boxes = [], [], []
    for pan, old, new, keep, count in edits:
        scope = panels[pan] if pan else root
        hits = [t for t in scope.iter(f"{{{SVG}}}text") if string(t) == old]
        assert len(hits) == count, (src, pan, old, len(hits))
        for t in hits:
            b_old = box(t)
            if keep == "end" and anchor(t) != "end":
                M = matrix(t.get("transform"))
                dw = (width(new, t) - width(old, t)) * (1 if anchor(t) == "start" else 0.5)
                M[:2, 2] -= dw * M[:2, 0]                               # along the direction of the text
                t.set("transform", fmt(M, t.get("transform")))
            leaf = next((e for e in t.iter() if (e.text or "") == old), None)
            assert leaf is not None, old
            leaf.text = new
            b_new = box(t)
            assert b_new[0] >= 0 and b_new[2] <= 518.4, (new, b_new)
            for u in texts:                                             # no new overlap with another text
                if u is t or not string(u).strip() or u.get("font-size") is None:    # (tick labels set in tspans)
                    continue
                o = box(u)
                hit = lambda b: b[0] < o[2] - 0.3 and o[0] < b[2] - 0.3 and b[1] < o[3] - 0.3 and o[1] < b[3] - 0.3
                assert not hit(b_new) or hit(b_old), f"{new!r} runs into {string(u)!r}"
            changed.append(t); boxes += [b_old, b_new]
            report.append(f"{pan or 'note'}: {old!r} -> {new!r}, box {b_old[2] - b_old[0]:.1f} x {b_old[3] - b_old[1]:.1f} -> "
                          f"{b_new[2] - b_new[0]:.1f} x {b_new[3] - b_new[1]:.1f} pt")
    svg, used = svg_finish.finish(ET.tostring(root, encoding="unicode"))
    open(dest, "w", encoding="utf-8").write(svg)
    # every changed pixel lies in the box of a changed text (old or new), padded by 1.5 pt
    tmp = tempfile.mkdtemp()
    a, b = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
    svg_finish._render(text, a); svg_finish._render(svg, b)
    A = np.asarray(Image.open(a).convert("L"), dtype=int); Bm = np.asarray(Image.open(b).convert("L"), dtype=int)
    assert A.shape == Bm.shape
    k = A.shape[1] / 518.4
    ch = np.abs(A - Bm) > 40
    allowed = np.zeros_like(ch)
    for x0, y0, x1, y1 in boxes:
        allowed[max(int((y0 - 1.5) * k), 0):int((y1 + 1.5) * k) + 1, max(int((x0 - 1.5) * k), 0):int((x1 + 1.5) * k) + 1] = True
    stray = ch & ~allowed
    assert not stray.any(), f"{src}: {int(stray.sum())} pixels changed outside the edited texts, e.g. at {np.argwhere(stray)[0] / k} pt"
    for name, g in panels.items():
        own = set(g.iter(f"{{{SVG}}}text"))
        n = sum(t in own for t in changed)
        if n:
            report.append(f"{name}: {n} texts changed")
    report.append(f"{os.path.basename(src)}: {int(ch.sum())} pixels changed, all inside the edited texts; fonts {sorted(used)}")
    return report


if __name__ == "__main__":
    src_dir, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for name, edits in EDITS.items():
        print("\n".join(edit(os.path.join(src_dir, name), os.path.join(out, name), edits)))
