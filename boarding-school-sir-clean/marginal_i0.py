"""1D profile log-lik vs I(0), for a few R(0) values.

β, γ anchored at the fixed-IC Poisson-SIR MLE (the only converged fit we
have for this data). Each curve fixes R(0) and sweeps I(0).
"""
import subprocess, tempfile, time
from pathlib import Path
import numpy as np, polars as pl
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"
MODEL = "boarding_school_sir_init.ir.json"
DATA  = ROOT/"data/in_bed.tsv"
GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

BETA_FIXED  = 1.9058652207
GAMMA_FIXED = 0.6559714166
R0_VALUES   = [0, 50, 100, 200, 300]
I0_GRID     = np.linspace(1, 30, 30)
N_PARTICLES = 1000
SEED = 42

def eval_ll(params):
    with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as f:
        for k,v in params.items():
            f.write(f"{k} = {v}\n")
        path = f.name
    try:
        res = subprocess.run(
            ["camdl", "pfilter", MODEL, "--params", path,
             "--data", str(DATA), "--particles", str(N_PARTICLES),
             "--seed", str(SEED), "--dt", "1.0"],
            capture_output=True, text=True, check=False, cwd=str(ROOT))
        return float([ln for ln in res.stdout.strip().splitlines() if ln.strip()][-1])
    except (ValueError, IndexError):
        return float("nan")

cache = OUT/"marginal_i0.npz"
if cache.exists():
    d = np.load(cache)
    lls_by_r0 = {r0: d[f"r0_{r0}"] for r0 in R0_VALUES}
else:
    t0 = time.time()
    lls_by_r0 = {}
    for r0 in R0_VALUES:
        lls = np.array([
            eval_ll({"beta": BETA_FIXED, "gamma": GAMMA_FIXED,
                     "I0": float(i0), "R_init": float(r0), "N0": 763})
            for i0 in I0_GRID
        ])
        lls_by_r0[r0] = lls
        print(f"  R(0)={r0}: ll range [{lls.min():.2f}, {lls.max():.2f}] "
              f"argmax at I0={I0_GRID[np.nanargmax(lls)]:.2f}, "
              f"max={lls.max():.2f}  ({time.time()-t0:.0f}s)")
    np.savez(cache, **{f"r0_{r0}": lls_by_r0[r0] for r0 in R0_VALUES})

fig, ax = plt.subplots(figsize=(9, 5))
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(R0_VALUES)))
for r0, c in zip(R0_VALUES, colors):
    lls = lls_by_r0[r0]
    i_max = np.nanargmax(lls)
    ax.plot(I0_GRID, lls, color=c, linewidth=1.6, marker="o",
            markersize=4, label=f"R(0) = {r0}")
    ax.scatter([I0_GRID[i_max]], [lls[i_max]], color=c, s=120, marker="*",
               edgecolor="white", linewidth=1.0, zorder=5)

# Vertical line at fixed-IC baseline's I(0) = 5
ax.axvline(5, color="#3498db", linestyle="--", alpha=0.7, linewidth=1.0,
           label="baseline I(0) = 5")
# Explicit legend entry for the star marker
import matplotlib.lines as mlines
star_handle = mlines.Line2D([], [], marker="*", markersize=12, linestyle="",
                            markerfacecolor="none", markeredgecolor=GREY,
                            markeredgewidth=1.2,
                            label="argmax I(0) per curve")
handles, labels = ax.get_legend_handles_labels()
handles.append(star_handle); labels.append("argmax I(0) per curve")
# Insert an invisible spacer at position 3 — the 4th row of the left
# column — so that with ncol=2 col-major fill, left column is
# [R0=0, R0=50, R0=100, blank], right column is [R0=200, R0=300,
# baseline, argmax]. Top rows of both columns (R0=0 and R0=200) are
# then visible and vertically aligned.
blank = mlines.Line2D([0], [0], color="white", linewidth=1.6)
handles = handles[:3] + [blank] + handles[3:]
labels  = labels[:3]  + ["          "] + labels[3:]
ax.set_xlabel("I(0)   initial infectious")
ax.set_ylabel("log-likelihood (1000-particle PF)")
ax.set_title(f"ll vs I(0) at β = {BETA_FIXED:.3f}, γ = {GAMMA_FIXED:.3f} "
             "(fixed-IC MLE), varying R(0)",
             fontsize=10, loc="left")
ax.legend(handles, labels, frameon=False, fontsize=9, loc="upper right",
          ncol=2, bbox_to_anchor=(1.0, 1.04))
ax.spines[["top","right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG/"marginal_i0.png", bbox_inches="tight")
plt.close(fig)
print(f"wrote {FIG/'marginal_i0.png'}")

# Print the argmax table
print("\nargmax I(0) per R(0):")
print(f"{'R(0)':>5}  {'best I(0)':>10}  {'ll at max':>10}  {'ll at I(0)=5':>14}")
for r0 in R0_VALUES:
    lls = lls_by_r0[r0]
    i_max = np.nanargmax(lls)
    # find I(0) = 5 index
    j5 = int(np.argmin(np.abs(I0_GRID - 5.0)))
    print(f"{r0:>5d}  {I0_GRID[i_max]:>10.2f}  {lls[i_max]:>10.2f}  {lls[j5]:>14.2f}")
