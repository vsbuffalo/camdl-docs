"""Prequential (one-step-ahead predictive) model comparison for the three fits.

Reads ll_increment from `camdl pfilter --trace` output for each model and plots:
  1) Per-day one-step-ahead log predictive score
  2) Cumulative prequential log-lik over time
  3) Total prequential log-score with per-model ESS floor (tail-depletion flag)
"""
import sys
from pathlib import Path
import numpy as np, polars as pl, matplotlib.pyplot as plt

HERE = Path(__file__).parent
BOOK = HERE.parent.parent
sys.path.insert(0, str(BOOK))
from styles import apply_rc
apply_rc()

models = [
    ("Poisson",                  "trace_poisson.tsv",   "#bdc3c7", 2),
    ("NegBin",                   "trace_negbin.tsv",    "#3498db", 3),
    ("Poisson + β proc-noise",   "trace_procnoise.tsv", "#e67e22", 3),
]
data = {}
for name, f, color, k in models:
    df = pl.read_csv(HERE/f, separator="\t", comment_prefix="#")
    data[name] = dict(df=df, color=color, k=k)

# Observed data for obs overlay
obs = pl.read_csv(BOOK/"data/in_bed.tsv", separator="\t")

from matplotlib.ticker import MaxNLocator
# Widen slightly + reserve right margin for the shared legend.
# Use 6 logical columns so row 0 can be 3 panels (2 cols each) and
# row 1 can be 2 panels (3 cols each) — gives equal-width bottom row.
fig = plt.figure(figsize=(14.5, 6.4))
gs = fig.add_gridspec(2, 6, height_ratios=[1, 1],
                       hspace=0.45, wspace=1.1, right=0.86)

# Panel A: observed counts for orientation
axA = fig.add_subplot(gs[0, 0:2])
axA.plot(obs["time"], obs["in_bed"], "o-", color="black", linewidth=1.2,
         clip_on=False, zorder=5)
axA.fill_between(obs["time"], 0, obs["in_bed"], alpha=0.08, color="black")
axA.set_xlabel("day"); axA.set_ylabel("boys in bed")
axA.set_title("Observed outbreak", fontsize=10, pad=4)
axA.set_ylim(bottom=0)
axA.margins(x=0)
axA.spines[["top","right"]].set_visible(False)

# Panel B: per-day one-step-ahead log predictive score
axB = fig.add_subplot(gs[0, 2:4])
for name, d in data.items():
    axB.plot(d["df"]["time"], d["df"]["ll_increment"], "-o",
             color=d["color"], linewidth=1.4, markersize=5,
             markerfacecolor="white", markeredgewidth=1.2,
             markeredgecolor=d["color"], label=name)
axB.axhline(0, color="#555", linewidth=0.5)
axB.set_xlabel("day t")
axB.set_ylabel(r"$\log\,\hat{p}(y_t \mid y_{1:t-1}, \hat\theta)$")
axB.set_title("Per-day prequential log score\n"
              "(one-step-ahead predictive, higher = better)",
              fontsize=10, pad=4)
axB.spines[["top","right"]].set_visible(False)

# Panel C: ESS per step — tail-depletion flag for log-score fragility
axC = fig.add_subplot(gs[0, 4:6])
for name, d in data.items():
    axC.plot(d["df"]["time"], d["df"]["ESS"], "-o",
             color=d["color"], linewidth=1.4, markersize=5,
             markerfacecolor="white", markeredgewidth=1.2,
             markeredgecolor=d["color"], label=name)
axC.axhline(500, color="#c0392b", linestyle=":", linewidth=1.0,
            label="ESS = 500 (10%)")
axC.set_xlabel("day t")
axC.set_ylabel("effective sample size")
axC.set_title("PF effective sample size per step\n"
              "(low ESS = log-score unstable)",
              fontsize=10, pad=4)
axC.spines[["top","right"]].set_visible(False)

# Panel D: cumulative prequential ll over time
axD = fig.add_subplot(gs[1, 0:3])
for name, d in data.items():
    ll_cum = d["df"]["ll_increment"].cum_sum().to_numpy()
    axD.plot(d["df"]["time"], ll_cum, "-o",
             color=d["color"], linewidth=1.6, markersize=5,
             markerfacecolor="white", markeredgewidth=1.2,
             markeredgecolor=d["color"], label=name)
axD.set_xlabel("day t")
axD.set_ylabel("cumulative prequential log-lik")
axD.set_title("Cumulative one-step-ahead log score over time "
              "(Dawid prequential decomposition)",
              fontsize=10, pad=4)
axD.spines[["top","right"]].set_visible(False)

# Panel E: total prequential score bar chart with ΔScore labels
axE = fig.add_subplot(gs[1, 3:6])
totals = np.array([d["df"]["ll_increment"].sum() for d in data.values()])
colors = [d["color"] for d in data.values()]
names  = list(data.keys())
ks     = [d["k"] for d in data.values()]
axE.bar(range(3), totals, color=colors)
axE.set_xticks(range(3))
axE.set_xticklabels([n.replace(" + ", "\n+ ") for n in names], fontsize=8.5)
axE.set_ylabel("prequential log-lik")
axE.set_title("Total prequential score\n(higher = better)",
              fontsize=10, pad=4)
pad = (totals.max() - totals.min()) * 0.15 + 0.5
axE.set_ylim(totals.min() - pad*1.2, totals.max() + pad*1.5)
for i, v in enumerate(totals):
    axE.text(i, v + pad*0.15, f"{v:.2f}", ha="center", fontsize=9)
best = int(np.argmax(totals))
axE.scatter([best], [totals[best] + pad*0.9], marker="*", s=200,
            color="red", edgecolor="white", linewidth=0.8, zorder=5)
axE.spines[["top","right"]].set_visible(False)

# Integer day ticks on every time-axis panel
for ax in (axA, axB, axC, axD):
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# Shared legend to the right of all panels
model_handles = []
model_labels  = []
for name, d in data.items():
    h, = plt.plot([], [], "-o", color=d["color"], linewidth=1.6,
                  markersize=6, markerfacecolor="white",
                  markeredgewidth=1.3, markeredgecolor=d["color"])
    model_handles.append(h); model_labels.append(name)
# Include the ESS threshold line in the shared legend
thr, = plt.plot([], [], color="#c0392b", linestyle=":", linewidth=1.0)
model_handles.append(thr); model_labels.append("ESS = 500 (10%)")

fig.legend(model_handles, model_labels, frameon=False, fontsize=9.5,
           loc="center left", bbox_to_anchor=(0.87, 0.5))

fig.suptitle("Prequential model comparison — one-step-ahead predictive scoring (Dawid 1984)",
             fontsize=11.5, y=1.00)
fig.savefig(HERE.parent/"figures/preq_comparison.png", bbox_inches="tight")
plt.close(fig)
print("wrote figures/preq_comparison.png")

# Text summary
print()
print(f"{'model':<28} {'k':>3} {'preq log-lik':>14} {'ΔPreq':>8} {'min ESS':>8}")
print("-"*72)
best_v = totals.max()
for name, d, tot in zip(names, data.values(), totals):
    min_ess = float(d["df"]["ESS"].min())
    print(f"{name:<28} {d['k']:>3} {tot:>14.3f} {best_v-tot:>8.2f} {min_ess:>8.0f}")

# Days where each model does best/worst
print("\nPer-day one-step-ahead log score (higher = better):")
tbl = pl.DataFrame({"day": data["Poisson"]["df"]["time"]})
for name, d in data.items():
    tbl = tbl.with_columns(d["df"]["ll_increment"].alias(name))
print(tbl)
