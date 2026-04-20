"""IF2 chain trajectories overlaid on the 2D log-lik surface.

Three panels: β×γ, β×k, γ×k. The β×k and γ×k panels use a *wider* grid
(k up to 4000) than the main profile-likelihood figure so the full scout
range is visible — scout chains reach k > 1500 during exploration.
The β×γ panel reuses the existing grid (no k axis).

Each panel shows:
  - 8 scout chains, every 100th iteration, light gray — dense exploration.
  - The best refine chain's full trace, darker gray — the "descent" onto
    the sloppy ridge.
  - Iter-0 circle, iter-final star, grid-max diamond.
"""
import subprocess, tempfile
from pathlib import Path
import numpy as np, polars as pl
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"
GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

MLE_PATH = next((ROOT/"results/fits").glob("fit_negbin-*/real/fit_42/refine/mle_params.toml"))
fit_dir  = MLE_PATH.parents[1]  # .../fit_42/

# ── Regenerate wider-k grids for the trajectory plot only (k up to 4000) ──
MODEL = "boarding_school_sir_negbin.ir.json"
DATA  = ROOT/"data/in_bed.tsv"
N_GRID_TRAJ = 20
N_PARTICLES_TRAJ = 1000
SEED_TRAJ = 42

def _eval_ll(params: dict) -> float:
    with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as f:
        for k, v in params.items():
            f.write(f"{k} = {v}\n")
        path = f.name
    try:
        res = subprocess.run(
            ["camdl", "pfilter", MODEL, "--params", path,
             "--data", str(DATA), "--particles", str(N_PARTICLES_TRAJ),
             "--seed", str(SEED_TRAJ), "--dt", "1.0"],
            capture_output=True, text=True, check=False, cwd=str(ROOT))
        return float([ln for ln in res.stdout.strip().splitlines() if ln.strip()][-1])
    except (ValueError, IndexError):
        return float("nan")

def _run_grid_traj(x_name, x_range, y_name, y_range, fixed, tag):
    p = OUT/f"profile_traj_{tag}.npz"
    if p.exists():
        print(f"  [{tag}] using cached {p.name}")
        return
    import time
    xs = np.linspace(*x_range, N_GRID_TRAJ)
    ys = np.linspace(*y_range, N_GRID_TRAJ)
    lls = np.full((len(xs), len(ys)), np.nan)
    t0 = time.time()
    for i, xv in enumerate(xs):
        for j, yv in enumerate(ys):
            lls[i, j] = _eval_ll({**fixed, x_name: float(xv), y_name: float(yv)})
        if (i+1) % 5 == 0:
            print(f"  [{tag}] {(i+1)*len(ys)}/{len(xs)*len(ys)}  {time.time()-t0:.0f}s")
    np.savez(p, xs=xs, ys=ys, lls=lls, x_name=x_name, y_name=y_name)

mle = {}
for ln in MLE_PATH.read_text().splitlines():
    s = ln.split("#",1)[0].strip()
    if s.startswith("[") or "=" not in s: continue
    k, v = s.split("=",1)
    try: mle[k.strip()] = float(v.strip())
    except: pass

# Wider-k grids for traj plot: k ∈ [5, 4000]
print("ensuring wider-k trajectory grids exist...")
_run_grid_traj("beta", (1.4, 3.0), "k", (5.0, 2000.0),
               {"I0": 5, "N0": 763, "gamma": mle["gamma"]}, "beta_k")
_run_grid_traj("gamma", (0.4, 1.0), "k", (5.0, 2000.0),
               {"I0": 5, "N0": 763, "beta": mle["beta"]}, "gamma_k")

def load_chains(stage_dir):
    dfs = []
    for cdir in sorted(stage_dir.glob("chain_*")):
        if not cdir.is_dir():
            continue
        d = pl.read_csv(cdir/"parameter_traces.tsv", separator="\t",
                        comment_prefix="#").sort("iteration")
        cid = int(cdir.name.split("_")[1])
        d = d.with_columns(pl.lit(cid).alias("chain"))
        dfs.append(d)
    return pl.concat(dfs) if dfs else pl.DataFrame()

refine = load_chains(fit_dir/"refine")
scout  = load_chains(fit_dir/"scout")

# Show the single best scout chain by final if2_perturbed_loglik
N_SCOUT_SHOW = 1
scout_final = scout.group_by("chain").agg(
    pl.col("if2_perturbed_loglik").last().alias("final_ll"),
    pl.col("iteration").max().alias("max_iter"),
)
top_scout_chains = (scout_final
                    .sort("final_ll", descending=True)
                    .head(N_SCOUT_SHOW)["chain"].to_list())
scout_top = scout.filter(pl.col("chain").is_in(top_scout_chains))
print(f"showing {len(top_scout_chains)} best scout chains: {sorted(top_scout_chains)}")
final_iter = refine["iteration"].max()
endpoints = refine.filter(pl.col("iteration") == final_iter)
# best = highest if2_perturbed_loglik (least negative)
best_chain = int(endpoints.sort("if2_perturbed_loglik", descending=True)["chain"][0])
print(f"best refine chain: {best_chain}")

best = refine.filter(pl.col("chain") == best_chain).sort("iteration")
print(f"plotting all {best.height} iterations of chain {best_chain}")
print(f"  iter range: {best['iteration'].min()} → {best['iteration'].max()}")
print(f"  β:  {float(best['beta'].min()):.3f} → {float(best['beta'].max()):.3f}")
print(f"  γ:  {float(best['gamma'].min()):.3f} → {float(best['gamma'].max()):.3f}")
print(f"  k:  {float(best['k'].min()):.1f} → {float(best['k'].max()):.1f}")
print(f"  ll: {float(best['if2_perturbed_loglik'].min()):.2f} → "
      f"{float(best['if2_perturbed_loglik'].max()):.2f}")


def plot_panel(ax, tag, x_lab, y_lab, x_col, y_col, x_mle, y_mle,
               use_traj_grid=False):
    fname = f"profile_traj_{tag}.npz" if use_traj_grid else f"profile_{tag}.npz"
    d = np.load(OUT/fname, allow_pickle=True)
    xs, ys, lls = d["xs"], d["ys"], d["lls"]
    ll_max = np.nanmax(lls)
    delta = np.clip(ll_max - lls, 0, 30)
    im = ax.pcolormesh(xs, ys, delta.T, cmap="viridis_r", shading="auto",
                       vmin=0, vmax=30)
    cs = ax.contour(xs, ys, (ll_max - lls).T, levels=[1.92, 4.61, 9.21],
                    colors=["cornflowerblue","#e67e22","#c0392b"], linewidths=1.0)
    ax.clabel(cs, fmt={1.92:"95%", 4.61:"99%", 9.21:"99.9%"}, fontsize=8)

    # scout trajectory — full trace (no thinning), mid-gray so it's
    # lighter than the refine trace but still readable.
    for cid, grp in scout_top.group_by("chain", maintain_order=True):
        g = grp.sort("iteration")
        xv = g[x_col].to_numpy(); yv = g[y_col].to_numpy()
        ax.plot(xv, yv, color="#808080", alpha=0.75, linewidth=0.7,
                zorder=3)

    # refine trajectory — solid dark gray for contrast against the white scout lines
    xvals = best[x_col].to_numpy()
    yvals = best[y_col].to_numpy()
    ax.plot(xvals, yvals, color="#2c2c2c", alpha=0.95, linewidth=1.4,
            zorder=5)

    # iter-0 marker (unfilled circle) for the refine chain
    ax.scatter([xvals[0]], [yvals[0]], marker="o", s=60,
               facecolor="white", edgecolor="#2c3e50", linewidth=1.2,
               zorder=6)

    # grid max
    gi, gj = np.unravel_index(np.nanargmax(lls), lls.shape)
    ax.scatter([xs[gi]], [ys[gj]], marker="D", s=110, color="#2ecc71",
               edgecolor="white", linewidth=1.4, zorder=7)
    # refine MLE
    ax.scatter([x_mle], [y_mle], marker="*", s=260, color="red",
               edgecolor="white", linewidth=1.0, zorder=8)

    # clip to grid extent so off-panel points don't skew axes
    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    ax.set_xlabel(x_lab); ax.set_ylabel(y_lab)
    ax.spines[["top","right"]].set_visible(False)
    return im

# Pair-plot layout — axes share both orientations so margins line up:
#
#               β x-axis        k x-axis
#   γ y-axis   [   β × γ   ]   [   k × γ   ]   ← γ shared along the row
#   k y-axis   [   β × k   ]   [  legend   ]   ← β shared down left col
#
# To align γ between β×γ (γ on y) and γ×k (γ on y), the γ×k panel is
# *transposed* into a "k×γ" panel: k on x-axis, γ on y-axis.
fig, axes = plt.subplots(
    2, 2, figsize=(12.5, 10.5),
    gridspec_kw={"hspace": 0.08, "wspace": 0.08},
    sharex=False, sharey=False,
)
axes[0,0].sharex(axes[1,0])   # β shared down column 0
axes[0,0].sharey(axes[0,1])   # γ shared across row 0

best_final = best.tail(1)
b_mle = float(best_final["beta"][0])
g_mle = float(best_final["gamma"][0])
k_mle = float(best_final["k"][0])

# Upper-left: β × γ (β x, γ y)
im = plot_panel(axes[0,0], "beta_gamma", "β", "γ", "beta", "gamma",
                b_mle, g_mle, use_traj_grid=False)
axes[0,0].set_title(f"β × γ  (k = {mle['k']:.1f} fixed)",
                    fontsize=10, loc="left")
plt.setp(axes[0,0].get_xticklabels(), visible=False)
axes[0,0].set_xlabel("")

# Upper-right: k × γ — transpose of γ × k so γ is on y-axis (shared w/ β×γ)
#   We call plot_panel_transposed which loads profile_traj_gamma_k.npz
#   but swaps the x and y coordinates.
def plot_panel_transposed(ax, tag, x_lab, y_lab, x_col, y_col,
                          x_mle, y_mle):
    """Same as plot_panel(..., use_traj_grid=True) but swaps x and y
    (so γ×k → k×γ)."""
    d = np.load(OUT/f"profile_traj_{tag}.npz", allow_pickle=True)
    xs, ys, lls = d["xs"], d["ys"], d["lls"]
    ll_max = np.nanmax(lls)
    delta = np.clip(ll_max - lls, 0, 30)
    # Transpose: original has xs on first axis, ys on second; now plot
    # ys on x-axis and xs on y-axis (so lls is used without transpose)
    im_ = ax.pcolormesh(ys, xs, delta, cmap="viridis_r", shading="auto",
                        vmin=0, vmax=30)
    cs = ax.contour(ys, xs, (ll_max - lls),
                    levels=[1.92, 4.61, 9.21],
                    colors=["cornflowerblue","#e67e22","#c0392b"], linewidths=1.0)
    ax.clabel(cs, fmt={1.92:"95%", 4.61:"99%", 9.21:"99.9%"}, fontsize=8)
    for cid, grp in scout_top.group_by("chain", maintain_order=True):
        g = grp.sort("iteration")
        # plot k on x, original-x-col (gamma) on y
        ax.plot(g["k"].to_numpy(), g[y_col].to_numpy(),
                color="#808080", alpha=0.75, linewidth=0.7, zorder=3)
    ax.plot(best["k"].to_numpy(), best[y_col].to_numpy(),
            color="#2c2c2c", alpha=0.95, linewidth=1.4, zorder=5)
    ax.scatter([best["k"].to_numpy()[0]], [best[y_col].to_numpy()[0]],
               marker="o", s=60, facecolor="white", edgecolor="#2c3e50",
               linewidth=1.2, zorder=6)
    # grid max (transposed coordinates)
    gi, gj = np.unravel_index(np.nanargmax(lls), lls.shape)
    ax.scatter([ys[gj]], [xs[gi]], marker="D", s=110, color="#2ecc71",
               edgecolor="white", linewidth=1.4, zorder=7)
    ax.scatter([x_mle], [y_mle], marker="*", s=260, color="red",
               edgecolor="white", linewidth=1.0, zorder=8)
    ax.set_xlim(ys.min(), ys.max())
    ax.set_ylim(xs.min(), xs.max())
    ax.set_xlabel(x_lab); ax.set_ylabel(y_lab)
    ax.spines[["top","right"]].set_visible(False)
    return im_

plot_panel_transposed(axes[0,1], "gamma_k",  "k", "γ", "k", "gamma",
                      k_mle, g_mle)
axes[0,1].set_title(f"k × γ  (β = {mle['beta']:.3f} fixed)",
                    fontsize=10, loc="left")
# hide γ yticks since shared with (0,0)
plt.setp(axes[0,1].get_yticklabels(), visible=False)
axes[0,1].set_ylabel("")
# hide k xticks since they're duplicated on (1,0)'s y-axis... actually k is
# x-axis here but y-axis on (1,0) — different orientation, so keep
# these ticks visible.

# Lower-left: β × k (β x, k y)
plot_panel(axes[1,0], "beta_k",  "β", "k", "beta", "k", b_mle, k_mle,
           use_traj_grid=True)
axes[1,0].set_title(f"β × k  (γ = {mle['gamma']:.3f} fixed)",
                    fontsize=10, loc="left")

# Lower-right: legend + colorbar
axes[1,1].axis("off")

cbar_ax = fig.add_axes([0.615, 0.10, 0.014, 0.30])
fig.colorbar(im, cax=cbar_ax, label="Δ log-lik from panel max (clipped at 30)")

import matplotlib.lines as mlines
handles = [
    mlines.Line2D([], [], color="#808080", alpha=0.85, linewidth=1.1,
                  label="scout chain (500 iter)"),
    mlines.Line2D([], [], color="#2c2c2c", linewidth=1.6,
                  label=f"refine chain {best_chain} (400 iter)"),
    mlines.Line2D([], [], marker="o", markersize=7, linestyle="",
                  markerfacecolor="white", markeredgecolor="#2c3e50",
                  markeredgewidth=1.1, label="refine iter 0"),
    mlines.Line2D([], [], marker="*", markersize=13, linestyle="",
                  markerfacecolor="red", markeredgecolor="white",
                  markeredgewidth=1.2, label="refine MLE"),
    mlines.Line2D([], [], marker="D", markersize=8, linestyle="",
                  markerfacecolor="#2ecc71", markeredgecolor="white",
                  markeredgewidth=1.2, label="grid max"),
    mlines.Line2D([], [], color="cornflowerblue", linewidth=1.1,
                  label="95% Δ-ll contour"),
    mlines.Line2D([], [], color="#e67e22", linewidth=1.1,
                  label="99% Δ-ll contour"),
    mlines.Line2D([], [], color="#c0392b", linewidth=1.1,
                  label="99.9% Δ-ll contour"),
]
axes[1,1].legend(handles=handles, loc="upper right", frameon=False,
                 fontsize=10, bbox_to_anchor=(0.98, 0.98))

fig.suptitle("IF2 chain trajectories on the NegBin log-lik surface — pair-plot layout",
             y=0.995, fontsize=11)
fig.savefig(FIG/"trajectory_overlay.png", bbox_inches="tight")
plt.close(fig)
print(f"\nwrote {FIG/'trajectory_overlay.png'}")
