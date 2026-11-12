"""Compare two SVG figures drawn by the scripts of figures/, for the rebuild check of zenodo_bundle.sh.

    python figures/svg_same.py <a.svg> <b.svg> <label>        # prints one line, exit 0 when both render the same

Two runs of the same script differ in what matplotlib and fontTools make new each time: the dc:date, the random ids of
markers, clip paths and glyphs, and the head.modified timestamp of every embedded font. The text is compared with these
normalised ("same text"), and both files are rendered by rsvg-convert at zoom 3 and compared pixel by pixel ("same
pixels"). An in-place script that redraws panels of a published file (fig3CD_coverage.py) writes new ids into the
redrawn panels, so for it only the pixels can agree.
"""
import base64, hashlib, io, os, re, sys, tempfile
import numpy as np
from PIL import Image
from fontTools.ttLib import TTFont
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish


def font_digest(m):
    f = TTFont(io.BytesIO(base64.b64decode(m.group(2))), recalcTimestamp=False)
    f["head"].modified = 0; f["head"].checkSumAdjustment = 0
    h = hashlib.sha1()
    for tag in sorted(t for t in f.keys() if t != "GlyphOrder"):
        h.update(tag.encode() + (f["head"].compile(f) if tag == "head" else f.getTableData(tag)))
    return f"{m.group(1)}{h.hexdigest()}"


def normalised(t):
    t = re.sub(r"<dc:date>[^<]*</dc:date>", "<dc:date/>", t)
    t = re.sub(r"(url\(data:font/woff2;base64,)([A-Za-z0-9+/=]+)", font_digest, t)
    ids = {}
    for i in re.findall(r'id="([^"]+)"', t):
        if re.fullmatch(r"[mp][0-9a-f]{10}|DejaVu[^\"]*|[A-Za-z]+-[0-9a-f]+|g[0-9a-f]{8,}", i) and i not in ids:
            ids[i] = f"id{len(ids)}"
    if ids:
        t = re.sub("|".join(re.escape(k) for k in sorted(ids, key=len, reverse=True)), lambda m: ids[m.group(0)], t)
    return t


a, b, label = sys.argv[1:4]
A, B = (open(p, encoding="utf-8").read() for p in (a, b))
tmp = tempfile.mkdtemp()
pa, pb = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(A, pa); svg_finish._render(B, pb)
IA, IB = (np.asarray(Image.open(p).convert("RGB"), dtype=int) for p in (pa, pb))
pixels = IA.shape == IB.shape and not (IA != IB).any()
text = normalised(A) == normalised(B)
print(f"{'same pixels' if pixels else 'DIFFERENT pixels'}{', same text' if text else ''}  {label}")
sys.exit(0 if pixels else 1)
