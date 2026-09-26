"""Figure 2 of the paper: M1-M3 stacked with panel letters A-C, as PNG and PDF.

    python figures/make_figure2.py <out_dir>    # run from benchmarks/ after make_main.py
"""
import os, sys, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.image as mpimg

out = sys.argv[1]
panels = [mpimg.imread(f"figures/main/{n}.png") for n in ("M1_kinetic_modules", "M2_acrr", "M3_runtime")]
heights = [p.shape[0] / p.shape[1] for p in panels]
fig, axs = plt.subplots(3, 1, figsize=(10, 10 * sum(heights)), gridspec_kw=dict(height_ratios=heights, hspace=0.02))
for ax, img, letter in zip(axs, panels, "ABC"):
    ax.imshow(img); ax.axis("off")
    ax.text(0.003, 0.995, letter, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="left")
fig.subplots_adjust(0, 0, 1, 1)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(out, f"Figure2_draft.{ext}"), dpi=200)
