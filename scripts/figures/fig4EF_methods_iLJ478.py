"""Figure 4E, F: the map of relations and methods (E) and what each tool computes, on one model (F).
fig4_figS1_S2_benchmark.py draws them into Figure 4 through draw(); run alone, this script writes them as Figure4EF_draft.

E  Rows: flux balance (N r = 0, linear) and mass action (polynomial in x); columns: reaction or metabolite space
   and complex space. Each cell names the relation and the matrix it is read from; A and Y^T carry reaction and
   metabolite space into complex space. Langary et al.'s pipeline is the dashed path: concordance, then the
   reaction graph (upstream sets) to kinetic modules, then the Shinar-Feinberg rule (coupled complexes that
   differ in one or two metabolites) to ACR/ACRR. The readout is zero rows (fixed) and proportional (flux
   balance) or equal (mass action, log space) rows (coupled). Below, the mass-action methods nested by what they
   report; the reference-equivalent mode is light green with a dashed edge, as in A and B.
F  The computations per tool and stage on iLJ478 (3440 complexes), each cell in the same slots: concordance
   (test, count, pairs skipped, time) and kinetic modules with ACR/ACRR (test, input, count, ACR rule, time).
   Sources: reference_final/reference_iLJ478.json (pairs tested, concordance time); the reference R stages
   timed on scratch (ref_rstage_timing_2026-11-03, median of three); sweep/cocoa_eff/cell_iLJ478_seed1234.log
   (COCOA.jl efficient=true pair funnel and concordance time); COCOA's kinetic stage rerun alone from the
   saved concordance labels (cocoa_kin_timing_2026-11-03: 2 rounds, warm median of three, same ACR, ACRR and
   giant module as the full run); Kimonu from figdata.json (timing_exclusive.sbatch: iterated route 2 rounds, median
   stage times; reference-equivalent mode in all; read from figdata.json). Kimonu's kinetic cell names what the iterated route
   computes: equal rows of ker W and of Y^T ker((Y Delta)^T), which contains B, not B itself.
   The dashed green line marks Kimonu's reference-equivalent mode: Kimonu's concordance, then the reference's
   kinetic stage, labelled on the right. It runs in a lane above the wall time of each cell.

    python figures/fig4EF_methods_iLJ478.py <out_dir>     # from benchmarks/; fig4_figS1_S2_benchmark.py imports draw()
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import to_rgb

INK, MUTE, RULE = "#1b2430", "#6b7480", "#c9ced4"
GREEN, BLUE, RED, GREY = "#4DAF4A", "#377EB8", "#E41A1C", "#8a929b"
tint = lambda c, t=0.85: tuple(t + (1 - t) * v for v in to_rgb(c))
HCD = 2.95                                # height of the C, D block in inches
REF_KIN = "3.5 min"                       # reference R stages, kinetic 209 s (incl. 3.4 s load) + ACR 2.5 s, median of three

def draw(ax, W=7.2, H=HCD):
    fd = json.load(open("figdata.json")); KD, KR = fd["kimonu_det"]["iLJ478"], fd["kimonu_ref"]["iLJ478"]   # Kimonu's times on iLJ478
    """Draw C and D on ax, whose data coordinates are inches (0..W, 0..H); returns the texts that carry math."""
    def box(x, y, w, h, face="white", edge=RULE, lw=0.6, r=0.04, z=1, ls="-"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=face, ec=edge, lw=lw, ls=ls, zorder=z))

    def text(x, y, s, size=6, color=INK, ha="center", va="center", weight="normal", z=5, **kw):
        ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight, zorder=z, linespacing=1.15, **kw)

    def arrow(p, q, color=INK, ls="-", lw=0.8, rad=0.0, z=4):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=6, color=color, lw=lw, ls=ls,
                                     connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0, zorder=z))

    # ---------------- E: the map
    x0 = 0.52; cw, ch, gx, gy = 1.08, 0.46, 0.46, 0.30
    ytop = H - 0.42 - ch                      # top row of cells
    cols = [x0, x0 + cw + gx]; rows = [ytop, ytop - ch - gy]
    text(cols[0] + cw / 2, H - 0.3, "reactions, metabolites", size=6, color=MUTE)
    text(cols[1] + cw / 2, H - 0.3, "complexes", size=6, color=MUTE)
    text(0.27, rows[0] + ch / 2, "flux\nbalance\n$\\mathbf{N}\\mathbf{r}=\\mathbf{0}$", size=6, color=MUTE)
    text(0.27, rows[1] + ch / 2, "mass\naction\npolynomial", size=6, color=MUTE)
    cells = {(0, 0): ("flux coupling", "$\\mathbf{K}=\\ker\\mathbf{N}$"), (0, 1): ("balanced, concordance", "$\\mathbf{M}=\\mathbf{A}\\mathbf{K}$"),
             (1, 0): ("ACR, ACRR", "$\\mathbf{B}$"), (1, 1): ("kinetic modules, ACRR", "$\\mathbf{Y}^\\top\\mathbf{B}$")}
    for (i, j), (a, b) in cells.items():
        box(cols[j], rows[i], cw, ch, face="white", edge="#9aa3ad")
        text(cols[j] + cw / 2, rows[i] + ch * 0.66, a, size=7)
        text(cols[j] + cw / 2, rows[i] + ch * 0.3, b, size=7)
    mid = lambda i: rows[i] + ch / 2
    arrow((cols[0] + cw, mid(0)), (cols[1], mid(0))); text(cols[0] + cw + gx / 2, mid(0) + 0.08, "$\\mathbf{A}$", size=7)
    arrow((cols[0] + cw, mid(1) + 0.1), (cols[1], mid(1) + 0.1)); text(cols[0] + cw + gx / 2, mid(1) + 0.19, "$\\mathbf{Y}^\\top$", size=7)
    # Langary et al.'s path: concordance -> reaction graph (upstream sets) -> kinetic modules -> Shinar-Feinberg rule -> ACR/ACRR
    arrow((cols[1] + cw * 0.5, rows[0]), (cols[1] + cw * 0.5, rows[1] + ch), color=BLUE, ls="--")
    text(cols[1] + cw * 0.5 + 0.05, rows[1] + ch + gy / 2, "reaction graph", size=6, color=BLUE, ha="left")
    arrow((cols[1], mid(1) - 0.04), (cols[0] + cw, mid(1) - 0.04), color=BLUE, ls="--")
    text(cols[0] + cw + gx / 2, mid(1) - 0.09, "Shinar–\nFeinberg", size=6, color=BLUE, va="top")
    text(x0 + cw + gx / 2, rows[1] - 0.15, "zero rows fixed, proportional (top) or equal (bottom) rows coupled", size=6, color=MUTE)

    # nested methods, mass action
    yl, hl = 0.1, 0.66
    names = [("reference-\nequivalent", "upstream sets,\nShinar–Feinberg", GREEN), ("iterated\nroute", "proven,\nevery $k$", GREEN),
             ("true\nrelations", ("two-point", "criterion", "(polynomial) or", "MILP sign", "certificate"), GREY), ("tangent\nroute", "complete,\ngeneric $k$", GREEN)]
    bw, sep = 0.66, 0.15
    xs = [0.05 + k * (bw + sep) for k in range(4)]
    for k, (nm, sub, col) in enumerate(names):                 # reference-equivalent: light green, dashed, as in A and B
        box(xs[k], yl, bw, hl, face=tint(col, 0.85 if k == 0 else 0.88), edge=col, lw=0.8 if k == 0 else 0.7, ls="--" if k == 0 else "-")
        text(xs[k] + bw / 2, yl + hl - 0.15, nm, size=7, weight="bold" if k == 1 else "normal")
        if isinstance(sub, tuple):                             # five lines 6 pt apart, the last 1.9 pt above the foot of the box
            for i, line in enumerate(reversed(sub)):
                text(xs[k] + bw / 2, yl + (1.9 + 6 * i) / 72, line, size=6, color=INK, va="baseline")
        else:
            text(xs[k] + bw / 2, yl + 0.18, sub, size=6, color=INK)
        if k < 3:
            text(xs[k] + bw + sep / 2, yl + hl / 2, "$\\subseteq$", size=9)
    text(0.05, yl + hl + 0.1, "mass action, from what is proven (left) to what is complete (right)", size=6, color=MUTE, ha="left")

    # ---------------- F: what each tool computes, iLJ478. Every cell fills the same slots, in the same order:
    # concordance: test, count, pairs skipped, -, wall time; kinetic: test, input, count, ACR rule, wall time.
    # A lane above the wall time carries the dashed line of Kimonu's reference-equivalent mode.
    xd = 3.2; wl, wc1, gut, wc2 = 0.55, 1.3, 0.12, 1.53
    yh = H - 0.42                              # header baseline
    text(xd + wl + wc1 / 2, yh + 0.12, "concordance", size=7, weight="bold")
    text(xd + wl + wc1 + gut + wc2 / 2, yh + 0.12, "kinetic modules, ACR, ACRR", size=7, weight="bold")
    n = fd["n_complexes"]["iLJ478"]                         # with the zero complex, as table S2
    text(xd + wl + (wc1 + gut + wc2) / 2, H - 0.12, f"iLJ478, {n} complexes, {n * (n - 1) / 2e6:.1f} M complex pairs", size=6, color=MUTE)
    rh = 0.74; gap = 0.06; step, lane = 0.112, 0.185      # slot spacing; lane height above the cell bottom
    SF = "ACR, ACRR by Shinar–Feinberg rule"
    ROWS = [("Reference", "(Langary et\nal., 2025)", BLUE,
             ["LPs per pair, float", "LPs on 3.54 M pairs", "2.37 M balanced pairs skipped", "", "5.1 d"],
             ["upstream sets, merged by overlap", "reuses concordance modules, no LP", "1 graph search per module, 1 round", SF, REF_KIN]),
            ("COCOA.jl", "", RED,
             ["LPs per pair, float, after sampling", "LPs on 90 k pairs", "5.82 M skipped by three filters", "", "3.1 h"],
             ["upstream sets, merged by overlap", "reuses concordance modules, no LP", "1 graph search per module, 2 rounds", SF, "4.5 s"]),
            ("Kimonu", "", GREEN,
             ["proportional rows of $\\mathbf{A}\\mathbf{K}$, exact", "1 nullspace, $\\ker\\mathbf{N}$", "no pair tested, one readout", "", f"{KD['sec_conc']:.1f} s"],
             ["equal rows of $\\ker\\mathbf{W}$, $\\mathbf{Y}^\\top\\ker((\\mathbf{Y}\\Delta)^\\top)$", "all of $\\ker\\mathbf{N}$ at one flux $\\mathbf{r}^{*}$, 1 LP",
              "2 nullspaces per round, 2 rounds", "ACR, ACRR from zero and equal rows", f"{KD['sec_kin'] + KD['sec_acr']:.1f} s"])]
    cell_y = []
    for k, (lab, cite, col, c1, c2) in enumerate(ROWS):
        y = yh - 0.06 - (k + 1) * rh - k * gap; cell_y.append(y)
        text(xd + 0.04, y + rh / 2 + (0.1 if cite else 0), lab, size=7, color=col, ha="left", weight="bold")
        if cite:
            text(xd + 0.04, y + rh / 2 - 0.06, cite, size=6, color=col, ha="left")
        for x, w, lines in ((xd + wl, wc1, c1), (xd + wl + wc1 + gut, wc2, c2)):
            box(x, y, w - 0.03, rh, face=tint(col, 0.9), edge=tint(col, 0.55), lw=0.5)
            top = y + rh - 0.06 - step / 2           # four slots from the top, the wall time at the bottom, so lines align
            for n, ln in enumerate(lines[:4]):
                text(x + 0.05, top - n * step, ln, size=6, ha="left")
            text(x + 0.05, y + 0.08, lines[4], size=6, ha="left")

    # Kimonu's reference-equivalent mode: Kimonu's concordance, then the reference's kinetic stage (all of it,
    # median of three, figdata.json). The line runs in the lane above the wall times.
    yK, yR = cell_y[2] + lane, cell_y[0] + lane
    xg = xd + wl + wc1 - 0.03 + (gut + 0.03) / 2              # middle of the gutter between the two columns
    xe = xd + wl + wc1 + gut + wc2 - 0.03 + 0.09              # just right of the kinetic column
    ax.plot([xd + wl - 0.04, xg, xg, xg + 0.05], [yK, yK, yR, yR], color=GREEN, lw=1.0, ls=(0, (3, 1.6)),
            solid_joinstyle="miter", dash_joinstyle="miter", zorder=6)
    arrow((xg + 0.05, yR), (xe, yR), color=GREEN, ls=(0, (3, 1.6)), lw=1.0, z=6)
    text(xe + 0.03, yR + 0.05, "Kimonu,\nreference-\nequivalent\nmode", size=6, color=GREEN, ha="left", va="bottom")
    text(xe + 0.03, yR + 0.01, f"{KR['seconds']:.1f} s", size=6, ha="left", va="top")

    for L, x in (("E", 0.0), ("F", xd - 0.05)):
        ax.text(x, H - 0.02, L, fontsize=18, fontweight="bold", family="Arial", va="top", ha="left", color=INK)
    return [t for t in ax.texts if "$" in t.get_text()]

if __name__ == "__main__":
    import svg_finish
    plt.rcParams.update({"font.size": 6.5, **svg_finish.mathtext_rc()})
    fig = plt.figure(figsize=(7.2, HCD))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 7.2); ax.set_ylim(0, HCD); ax.axis("off")
    draw(ax)
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(out, f"Figure4EF_draft.{ext}"), dpi=300)
