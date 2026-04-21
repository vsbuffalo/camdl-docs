"""Two-panel diagnostic for the simple SIR IF2 fit.

Panel A: unconditional forward posterior predictive.
Panel B: smoothing paths over latent I from the particle filter.

Both use --backend chain_binomial --dt 1.0 to match the fit.
"""
import subprocess, shlex
from pathlib import Path
import polars as pl, matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"; FIG.mkdir(exist_ok=True)
OUT  = ROOT/"output";  OUT.mkdir(exist_ok=True)

# Locate the refine MLE under the unified output tree.
mle_candidates = sorted((ROOT/"results/fits").glob("fit-*/real/fit_42/refine/mle_params.toml"))
if not mle_candidates:
    raise SystemExit("No fit found. Run: camdl fit run fit.toml --seed 42")
MLE = mle_candidates[-1]
DATA= ROOT/"data/in_bed.tsv"

GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

UC_TSV = OUT/"ppc_unconditional.tsv"
SM_TSV = OUT/"smoothing_paths.tsv"

# Panel A — unconditional forward, backend matches fit
subprocess.run(shlex.split(
    f"camdl simulate boarding_school_sir.ir.json --params {MLE} "
    f"--replicates 200 --seed 1 --scenario baseline "
    f"--backend chain_binomial --dt 1.0 --obs {UC_TSV}"),
    check=True, capture_output=True, cwd=str(ROOT))

# Panel B — smoothing paths from PF
subprocess.run(shlex.split(
    f"camdl pfilter boarding_school_sir.ir.json --params {MLE} "
    f"--data {DATA} --particles 5000 --seed 42 "
    f"--save-paths 200 {SM_TSV} --dt 1.0"),
    check=True, capture_output=True, cwd=str(ROOT))

obs = pl.read_csv(DATA, separator="\t")
uc  = pl.read_csv(UC_TSV, separator="\t")
sm  = pl.read_csv(SM_TSV, separator="\t")

def ribbon(df, time_col, value_col, ax, color, label):
    q = df.group_by(time_col).agg(
        pl.col(value_col).quantile(0.05).alias("lo"),
        pl.col(value_col).quantile(0.5).alias("med"),
        pl.col(value_col).quantile(0.95).alias("hi"),
    ).sort(time_col)
    ax.fill_between(q[time_col], q["lo"], q["hi"], color=color, alpha=0.22, linewidth=0, edgecolor="none")
    ax.plot(q[time_col], q["med"], color=color, linewidth=1.8, label=label)

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), sharey=True)
ribbon(uc, "time", "in_bed", axes[0], "#3366cc",
       "unconditional forward\n(5–95% + median)")
axes[0].scatter(obs["time"], obs["in_bed"], color="black", s=30, zorder=5, label="observed")
axes[0].set_title("A. Unconditional posterior predictive\n(camdl simulate --backend chain_binomial, 200 replicates)",
                  fontsize=10, loc="left")
axes[0].set_xlabel("Day"); axes[0].set_ylabel("Boys in bed")
axes[0].legend(frameon=False, fontsize=9, loc="upper right")
axes[0].spines[["top","right"]].set_visible(False)

ribbon(sm, "time", "I", axes[1], "#e74c3c",
       "smoothing latent I\n(5–95% + median)")
axes[1].scatter(obs["time"], obs["in_bed"], color="black", s=30, zorder=5, label="observed")
axes[1].set_title("B. Smoothing paths over latent I\n(camdl pfilter --save-paths, 5000 particles → 200 paths)",
                  fontsize=10, loc="left")
axes[1].set_xlabel("Day")
axes[1].legend(frameon=False, fontsize=9, loc="upper right")
axes[1].spines[["top","right"]].set_visible(False)

fig.suptitle("Simple SIR + Poisson obs, IF2 MLE — backend-matched chain-binomial throughout",
             fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(FIG/"two_panel.png", bbox_inches="tight")
plt.close(fig)
print(f"wrote {FIG/'two_panel.png'}")
