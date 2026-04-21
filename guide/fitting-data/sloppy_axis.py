"""Predictions across the sloppy k-axis.

We established that k is independently weakly identified in the NegBin fit
(the (β, γ) optimum doesn't shift with k). That means many (β̂, γ̂, k̂)
points along the k-ridge have essentially equivalent log-likelihood.
If k is truly sloppy, the *predictions* at those points should also be
nearly identical — a sloppy direction by definition contributes nothing
to the observables.

This script:
1. Picks 5 points along the k-ridge at (β, γ) fixed near the grid max.
2. Simulates 200 chain-binomial forward replicates at each.
3. Plots 5 overlaid PPC envelopes to see whether predictions overlap.
"""
import subprocess, shlex, tempfile
from pathlib import Path
import polars as pl, numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"

GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

# Anchor (β, γ) at the grid-max basin — not refine's reported MLE, since we
# want the sloppy direction in its cleanest location.
BETA_FIXED  = 2.047
GAMMA_FIXED = 0.632
K_VALUES = [20, 50, 100, 200, 400]
N_REPS = 200
SEED = 1
COLORS = ["#6a51a3", "#3182bd", "#31a354", "#fd8d3c", "#de2d26"]

obs = pl.read_csv(ROOT/"data/in_bed.tsv", separator="\t")

def ribbon_stats(df, col="in_bed"):
    return df.group_by("time").agg(
        pl.col(col).quantile(0.05).alias("lo"),
        pl.col(col).quantile(0.5).alias("med"),
        pl.col(col).quantile(0.95).alias("hi"),
    ).sort("time")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8),
                         gridspec_kw={"width_ratios": [1.3, 1]})
ax = axes[0]

# For each k, simulate forward, plot ribbon
summary_rows = []
for k_val, color in zip(K_VALUES, COLORS):
    with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as f:
        f.write(f"beta = {BETA_FIXED}\ngamma = {GAMMA_FIXED}\nk = {k_val}\n"
                f"I0 = 5\nN0 = 763\n")
        params = f.name
    out_tsv = OUT/f"sloppy_k_{k_val}.tsv"
    subprocess.run(shlex.split(
        f"camdl simulate boarding_school_sir_negbin.ir.json --params {params} "
        f"--replicates {N_REPS} --seed {SEED} --scenario baseline "
        f"--backend chain_binomial --dt 1.0 --obs {out_tsv}"),
        check=True, capture_output=True, cwd=str(ROOT))
    df = pl.read_csv(out_tsv, separator="\t")
    q = ribbon_stats(df)
    ax.fill_between(q["time"], q["lo"], q["hi"], color=color, alpha=0.15, linewidth=0, edgecolor="none")
    ax.plot(q["time"], q["med"], color=color, linewidth=1.6,
            label=f"k = {k_val}")
    for t in q["time"].to_list():
        row = q.filter(pl.col("time") == t).to_dicts()[0]
        summary_rows.append({"k": k_val, "time": t,
                             "lo": row["lo"], "med": row["med"], "hi": row["hi"]})

ax.scatter(obs["time"], obs["in_bed"], color="black", s=28, zorder=6,
           label="observed")
ax.set_xlabel("Day"); ax.set_ylabel("Boys in bed")
ax.set_title("PPC envelopes at 5 points along the k-ridge\n"
             f"β = {BETA_FIXED} and γ = {GAMMA_FIXED} fixed at the grid-max basin",
             fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.spines[["top","right"]].set_visible(False)

# Right panel: PPC envelope WIDTH at peak vs k — quantifies how much (if
# any) the predictive spread changes along the sloppy direction
ax = axes[1]
summ = pl.DataFrame(summary_rows)
summ = summ.with_columns((pl.col("hi") - pl.col("lo")).alias("width"))
peak_day = int(obs.filter(pl.col("in_bed") == obs["in_bed"].max())["time"][0])
peak_widths = summ.filter(pl.col("time") == peak_day)
ax.plot(peak_widths["k"], peak_widths["width"], marker="o",
        color="#2c3e50", linewidth=2, markersize=8,
        markeredgecolor="white", markeredgewidth=1)
ax.set_xlabel("k (along the sloppy axis)")
ax.set_ylabel(f"90% PPC envelope width at day {peak_day}")
ax.set_xscale("log")
ax.set_title("How much do predictions change along the k-ridge?", fontsize=10, loc="left")
ax.spines[["top","right"]].set_visible(False)
ax.annotate("small k → larger obs variance\nlarge k → Poisson limit",
            xy=(0.97, 0.92), xycoords="axes fraction",
            fontsize=8, color=GREY, verticalalignment="top",
            horizontalalignment="right")

fig.suptitle("Sloppy k axis — predictions vary modestly across the ridge",
             fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(FIG/"sloppy_axis.png", bbox_inches="tight")
plt.close(fig)
print(f"wrote {FIG/'sloppy_axis.png'}")
print("\npeak-day envelope widths by k:")
print(peak_widths.select(["k", "lo", "med", "hi", "width"]))
