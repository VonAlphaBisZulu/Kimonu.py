"""Re-fit every gray reverse path of an Illustrator network SVG as a true offset of its black partner.

Each gray path (stroke #bbb) is matched to the nearest black path, its ends are projected onto that
path, and the offset of the black curve between those projections, at one common distance (the
median of the current gaps unless given), is fitted by cubics with the black curve's tangents,
split until each is within 0.05 of the exact offset. Gray arrowheads move with the end they sit on. Everything else in the file is untouched.

    python even_offsets.py <in.svg> <out.svg> [distance]
"""
import re, sys
import numpy as np

src, out = sys.argv[1], sys.argv[2]
s = open(src).read()
GRAY = r'stroke="#bbb"'

def parse(d):
    """'M x,y c ...' / 'M x,y C ...' / 'M x,y h ...' -> absolute cubic control points (4 x 2)."""
    nums = [float(v) for v in re.findall(r'-?\d*\.?\d+(?:e-?\d+)?', d)]
    p0 = np.array(nums[:2]); cmd = re.findall(r'[a-zA-Z]', d)[1]
    if cmd == 'c': return np.array([p0, p0 + nums[2:4], p0 + nums[4:6], p0 + nums[6:8]])
    if cmd == 'C': return np.array([p0, nums[2:4], nums[4:6], nums[6:8]])
    if cmd == 'h': p3 = p0 + [nums[2], 0]
    elif cmd == 'v': p3 = p0 + [0, nums[2]]
    else: raise ValueError(d)
    return np.array([p0, p0 + (p3 - p0) / 3, p0 + 2 * (p3 - p0) / 3, p3])

def line(x1, y1, x2, y2):
    a, b = np.array([x1, y1]), np.array([x2, y2]); return np.array([a, a + (b - a) / 3, a + 2 * (b - a) / 3, b])
B = lambda P, t: ((1 - t) ** 3)[:, None] * P[0] + (3 * (1 - t) ** 2 * t)[:, None] * P[1] + (3 * (1 - t) * t ** 2)[:, None] * P[2] + (t ** 3)[:, None] * P[3]
dB = lambda P, t: (3 * (1 - t) ** 2)[:, None] * (P[1] - P[0]) + (6 * (1 - t) * t)[:, None] * (P[2] - P[1]) + (3 * t ** 2)[:, None] * (P[3] - P[2])
T = np.linspace(0, 1, 801)
def normal(P, t):
    d = dB(P, t); d /= np.linalg.norm(d, axis=1, keepdims=True); return np.stack([-d[:, 1], d[:, 0]], 1)

black = []
for m in re.finditer(r'<path d="([^"]+)"[^>]*stroke="#000"[^>]*/>', s): black.append(parse(m.group(1)))
for m in re.finditer(r'<line x1="([\d.]+)" y1="([\d.]+)" x2="([\d.]+)" y2="([\d.]+)"[^>]*stroke="#000"', s):
    black.append(line(*map(float, m.groups())))
grays = [(m, parse(m.group(1))) for m in re.finditer(r'<path d="([^"]+)"[^>]*' + GRAY + r'[^>]*/>', s)]

def nearest(P, Q):
    """black curve closest to gray curve Q, the projected parameters of Q's ends, the side and the gaps"""
    q = B(Q, T[::20])
    best = None
    for Pb in black:
        pts = B(Pb, T); dist = np.linalg.norm(q[:, None] - pts[None], axis=2)
        score = dist.min(1).mean()
        if best is None or score < best[0]: best = (score, Pb, dist.argmin(1))
    _, Pb, idx = best
    t0, t1 = T[idx[0]], T[idx[-1]]
    nb = normal(Pb, T[idx]); side = np.sign(np.mean(np.sum((q - B(Pb, T[idx])) * nb, 1)))
    return Pb, t0, t1, side, np.abs(np.sum((q - B(Pb, T[idx])) * nb, 1))

fits = [nearest(None, Q) for _, Q in grays]
gap = float(sys.argv[3]) if len(sys.argv) > 3 else float(np.median(np.concatenate([f[4] for f in fits])))

def fit(Pb, t0, t1, side):
    t = np.linspace(t0, t1, 200); pts = B(Pb, t) + side * gap * normal(Pb, t)
    u = np.linspace(0, 1, len(t))                    # chord parameters would do too; the offset is smooth
    d0 = dB(Pb, np.array([t0]))[0]; d1 = dB(Pb, np.array([t1]))[0]
    d0 *= np.sign(t1 - t0) / np.linalg.norm(d0); d1 *= np.sign(t1 - t0) / np.linalg.norm(d1)
    p0, p3 = pts[0], pts[-1]
    b0, b1, b2, b3 = (1 - u) ** 3, 3 * (1 - u) ** 2 * u, 3 * (1 - u) * u ** 2, u ** 3
    A = np.stack([(b1[:, None] * d0).ravel(), (-b2[:, None] * d1).ravel()], 1)
    rhs = (pts - b0[:, None] * p0 - (b1 + b2)[:, None] * 0 - b1[:, None] * p0 - b2[:, None] * p3 - b3[:, None] * p3).ravel()
    a, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    a = np.maximum(a, 1e-3)
    C = np.array([p0, p0 + a[0] * d0, p3 - a[1] * d1, p3])
    err = np.abs(np.linalg.norm(B(C, np.linspace(0, 1, 200))[:, None] - pts[None], axis=2).min(1)).max()
    return C, err

def pieces(Pb, t0, t1, side, tol=0.05, depth=0):
    """split at the parameter midpoint until every cubic is within tol of the exact offset"""
    C, err = fit(Pb, t0, t1, side)
    if err <= tol or depth == 4: return [C], err
    tm = (t0 + t1) / 2
    (a, ea), (b, eb) = pieces(Pb, t0, tm, side, tol, depth + 1), pieces(Pb, tm, t1, side, tol, depth + 1)
    return a + b, max(ea, eb)

new = s; report = []
for (m, Q), (Pb, t0, t1, side, gaps) in zip(grays, fits):
    Cs, err = pieces(Pb, t0, t1, side); C = np.array([Cs[0][0], Cs[0][1], Cs[-1][2], Cs[-1][3]])
    d = f"M{Cs[0][0,0]:.3f},{Cs[0][0,1]:.3f}" + "".join(f"C{c[1,0]:.3f},{c[1,1]:.3f} {c[2,0]:.3f},{c[2,1]:.3f} {c[3,0]:.3f},{c[3,1]:.3f}" for c in Cs)
    new = new.replace(m.group(0), m.group(0).replace(m.group(1), d), 1)
    # a gray head at either end moves with that end
    for old_end, new_end in ((Q[0], C[0]), (Q[3], C[3])):
        for h in re.finditer(r'<polygon points="([^"]+)" fill="#bbb"/>', new):
            p = np.array([float(v) for v in h.group(1).split()]).reshape(-1, 2)
            if np.linalg.norm(p - old_end, axis=1).min() < 9:
                p = p + (new_end - old_end)
                new = new.replace(h.group(0), '<polygon points="' + " ".join(f"{v:.3f}" for v in p.ravel()) + '" fill="#bbb"/>', 1)
                break
    report.append(f"gap was {gaps.min():.2f}-{gaps.max():.2f}, now {gap:.2f} +- {err:.3f} ({len(Cs)} cubic{'s' if len(Cs) > 1 else ''})")
open(out, "w").write(new)
print(f"common distance {gap:.2f}"); print("\n".join(report))
