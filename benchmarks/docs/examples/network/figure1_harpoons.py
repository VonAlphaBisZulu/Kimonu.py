"""Figure 1: harpoons with the flush-edge heads of Figure 2 (figure2.py, pair()), in A, B and the D thumbnails.

Each shaft now runs to its tip, and its head sits half a stroke out from the shaft's centre line, towards the barb,
so the head's flat edge closes out with the outer edge of the shaft. Heads keep their size and the longitudinal
position of their tip; a curved shaft is continued along its end tangent, and its head turned onto that tangent.
Shafts keep their stroke widths (1.2 pt in A and B, .94 / .564 in the thumbnails). A head that is flush already
is left alone, so a second run changes nothing."""
import math, re, sys
src, out = sys.argv[1:3]
s = open(src).read()
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
ELEM = r"<(polygon|path|line)\b[^>]*/>"

def attr(el, a):
    m = re.search(r'\s%s="([^"]*)"' % a, el); return m.group(1) if m else None
def f3(x):
    return f"{x:.3f}".rstrip("0").rstrip(".")
def dot(a, b): return a[0] * b[0] + a[1] * b[1]
def sub(a, b): return (a[0] - b[0], a[1] - b[1])
def add(a, b, k=1): return (a[0] + k * b[0], a[1] + k * b[1])
def unit(a): L = math.hypot(*a); return (a[0] / L, a[1] / L)

def segments(d):
    """absolute cubic segments [p0, c1, c2, p3] of the M/L/H/V/C/S path data Illustrator writes"""
    toks = re.findall(r"[MmLlHhVvCcSsZz]|" + NUM, d)
    if any(re.match(r"[QqTtAa]", t) for t in toks) or re.search(r"[QqTtAa]", d): raise ValueError(d)
    i, cmd, cur, out_ = 0, None, (0, 0), []
    def take(k):
        nonlocal i
        v = [float(t) for t in toks[i:i + k]]; i += k; return v
    while i < len(toks):
        if re.match(r"[A-Za-z]", toks[i]): cmd = toks[i]; i += 1
        o = cur if cmd.islower() else (0, 0); C = cmd.upper()
        if C == "Z": continue
        if C == "M": cur = add(o, take(2)); cmd = "l" if cmd == "m" else "L"; continue
        if C in "LHV":
            if C == "L": p = add(o, take(2))
            elif C == "H": x, = take(1); p = (o[0] + x if cmd == "h" else x, cur[1])
            else: y, = take(1); p = (cur[0], o[1] + y if cmd == "v" else y)
            out_.append([cur, cur, p, p])
        elif C == "C":
            a = take(6); out_.append([cur, add(o, a[:2]), add(o, a[2:4]), add(o, a[4:])])
        else:                                                # S: first control point mirrors the last one
            a = take(4); out_.append([cur, add(cur, sub(cur, out_[-1][2])), add(o, a[:2]), add(o, a[2:])])
        cur = out_[-1][3]
    return out_
def bez(seg, u):
    return tuple((1 - u) ** 3 * seg[0][k] + 3 * (1 - u) ** 2 * u * seg[1][k] + 3 * (1 - u) * u ** 2 * seg[2][k] + u ** 3 * seg[3][k]
                 for k in (0, 1))

def shaft_ends(el):
    """[(which, end point, direction of travel out of that end, segment ending there)]"""
    if el.startswith("<line"):
        p0 = (float(attr(el, "x1")), float(attr(el, "y1"))); p1 = (float(attr(el, "x2")), float(attr(el, "y2")))
        t = unit(sub(p1, p0)); return [("start", p0, (-t[0], -t[1]), [p1, p1, p0, p0]), ("end", p1, t, [p0, p0, p1, p1])]
    sg = segments(attr(el, "d"))
    if not sg: return []
    first, last = sg[0], sg[-1]
    t0 = next(unit(sub(first[0], first[k])) for k in (1, 2, 3) if math.hypot(*sub(first[k], first[0])) > 1e-9)
    t1 = next(unit(sub(last[3], last[k])) for k in (2, 1, 0) if math.hypot(*sub(last[3], last[k])) > 1e-9)
    return [("start", first[0], t0, first[::-1]), ("end", last[3], t1, last)]

def harpoon(el):
    """tip, base, barb, length and barb of a half head (right angle at the base), or None"""
    v = [float(x) for x in re.findall(NUM, attr(el, "points") or "")]
    p = list(zip(v[::2], v[1::2]))
    if len(p) == 4 and p[0] == p[3]: p = p[:3]
    if len(p) != 3: return None
    for k in range(3):
        B, X, Y = p[k], p[(k + 1) % 3], p[(k + 2) % 3]
        for T, P in ((X, Y), (Y, X)):
            a, b = sub(T, B), sub(P, B)
            if abs(dot(a, b)) < .03 * math.hypot(*a) * math.hypot(*b) and math.hypot(*a) > math.hypot(*b):
                return T, B, P, math.hypot(*a), math.hypot(*b)

els = [(m.start(), m.end(), m.group(0)) for m in re.finditer(ELEM, s)]
shafts = []
for i, j, el in els:
    if el.startswith("<polygon") or attr(el, "fill") != "none" or not attr(el, "stroke"): continue
    try: ends = shaft_ends(el)
    except (ValueError, IndexError, StopIteration): continue
    shafts += [dict(i=i, el=el, w=float(attr(el, "stroke-width") or 1), stroke=attr(el, "stroke"), which=w_, e=e, t=t, seg=sg)
               for w_, e, t, sg in ends]

edits, report, done = [], [], 0
for i, j, el in els:
    h = el.startswith("<polygon") and harpoon(el)
    if not h: continue
    T, B, P, L, b = h
    fill = attr(el, "fill") or "#000"
    v, nb = unit(sub(T, B)), unit(sub(P, B))
    # its shaft: same colour, travelling into the head, centre line within a stroke of the flat edge on the barb side
    cand = [(math.hypot(*sub(sh["e"], add(B, sub(T, B), .5))), sh) for sh in shafts if sh["stroke"] == fill and dot(sh["t"], v) > .95]
    flush = [sh for _, sh in cand if abs(dot(sub(T, sh["e"]), sh["t"])) < 1e-2 and abs(dot(sub(sh["e"], T), nb) + sh["w"] / 2) < 1e-2]
    if flush: done += 1; continue
    cand = [(dist, sh) for dist, sh in cand if dist < L and 0 < dot(sub(sh["e"], B), nb) < sh["w"]]
    if not cand: sys.exit(f"no shaft for the head {el}")
    _, sh = min(cand, key=lambda c: c[0])
    e, t, w = sh["e"], sh["t"], sh["w"]
    n = (t[1], -t[0]) if dot((t[1], -t[0]), nb) > 0 else (-t[1], t[0])  # barb side, square to the shaft
    run = dot(sub(T, e), t)                                              # shaft end to tip, along the shaft
    e2 = add(e, t, run)
    T2 = add(e2, n, w / 2); B2 = add(T2, t, -L); P2 = add(B2, n, b)
    pts = " ".join(f"{f3(x)} {f3(y)}" for x, y in (T2, B2, P2, T2))
    edits.append((i, j, re.sub(r'points="[^"]*"', f'points="{pts}"', el)))
    # how far the shaft's centre line strays from straight under the head (curved shafts that ran into it)
    stray = max([abs(dot(sub(q, e), n)) for q in (bez(sh["seg"], k / 400) for k in range(401))
                 if 0 <= dot(sub(e, q), t) <= L - run] or [0])
    report.append((el, sh, run, math.degrees(math.asin(max(-1, min(1, v[0] * t[1] - v[1] * t[0])))), stray,
                   math.hypot(*sub(T2, T))))
    sh["new_end"] = e2

# shafts: lines get the new end point, paths a straight piece along their end tangent
for sh in shafts:
    if "new_end" not in sh: continue
    el = sh.get("el2", sh["el"]); x, y = sh["new_end"]
    if el.startswith("<line"):
        k = "1" if sh["which"] == "start" else "2"
        el = re.sub(rf'x{k}="[^"]*" y{k}="[^"]*"', f'x{k}="{f3(x)}" y{k}="{f3(y)}"', el)
    else:
        d = attr(el, "d")
        if sh["which"] == "start":
            m = re.match(r"\s*M\s*(" + NUM + r")[\s,]*(" + NUM + r")", d)
            d = f"M{f3(x)},{f3(y)}L{m.group(1)},{m.group(2)}" + d[m.end():]
        else:
            d = d + f"L{f3(x)},{f3(y)}"
        el = re.sub(r'\sd="[^"]*"', f' d="{d}"', el, count=1)
    for o in shafts:
        if o["i"] == sh["i"]: o["el2"] = el
for o in {o["i"]: o for o in shafts if "el2" in o}.values():
    edits.append((o["i"], o["i"] + len(o["el"]), o["el2"]))
for i, j, new in sorted(edits, reverse=True):
    s = s[:i] + new + s[j:]
open(out, "w").write(s)

print(f"{len(report)} heads moved, {done} flush already")
for el, sh, run, ang, stray, moved in report:
    print(f"  {attr(el, 'points')[:24]:24s} {attr(el, 'fill') or '#000':8s} w={sh['w']:<5} shaft {sh['which']:5s} +{run:.2f} pt, "
          f"turned {ang:+.1f} deg, tip moved {moved:.2f} pt, shaft strays {stray:.3f} pt under the head")
