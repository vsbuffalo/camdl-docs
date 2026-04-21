"""Poisson SIR chain-trajectory overlay — analog of trajectory_overlay.py
(NegBin model) but for the 4-parameter Poisson free-IC fit
(fit_free_init_v2), using the two 2D profile likelihood surfaces we
already computed as the backgrounds.

Two panels: (β, γ) profile with chain trajectories in β×γ space, and
(I(0), R(0)) profile with chain trajectories in I(0)×R(0) space.
Scout chains thinned to every 20 iterations, one chain highlighted
(best by final if2_perturbed_loglik)."""
from pathlib import Path
import polars as pl, numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"
GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

EDGE = GREY
EDGE_LW_STAR = 0.6
EDGE_LW_SQ   = 0.4

# ── 1. Load scout chains from the failed 4-param Poisson fit ──
fit_dir = next((ROOT/"results/fits").glob("fit_free_init_v2-*/real/fit_42"))
scout_dir = fit_dir/"scout"

dfs = []
for c in sorted(scout_dir.glob("chain_*")):
    if not c.is_dir(): continue
    # loglik column has "NA" strings where PF eval was skipped; parse as
    # null so it stays numeric instead of being inferred as a string.
    d = pl.read_csv(c/"parameter_traces.tsv", separator="\t",
                    comment_prefix="#", null_values="NA").sort("iteration")
    cid = int(c.name.split("_")[1])
    dfs.append(d.with_columns(pl.lit(cid).alias("chain")))
scout = pl.concat(dfs)
scout_thin = scout.filter(pl.col("iteration") % 20 == 0)

# Best chain = highest clean PF log-lik at the *final* PF-eval iteration
# (loglik is NA except every 10 iters). Matches camdl's "best chain"
# reporting in the fit log.
last_ll_by_chain = (scout
    .filter(pl.col("loglik").is_not_null())
    .sort("iteration")
    .group_by("chain", maintain_order=True)
    .agg(pl.col("loglik").last().alias("final_ll")))
best_chain = int(last_ll_by_chain.sort("final_ll", descending=True)["chain"][0])
best_ll_clean = float(last_ll_by_chain.sort("final_ll", descending=True)["final_ll"][0])
best = scout.filter(pl.col("chain") == best_chain).sort("iteration")
best_final = best.tail(1)
print(f"best scout chain: {best_chain}  (best clean ll = {best_ll_clean:.2f})")
print(f"  endpoint: β={float(best_final['beta'][0]):.3f}, "
      f"γ={float(best_final['gamma'][0]):.3f}, "
      f"I(0)={float(best_final['I0'][0]):.2f}, "
      f"R(0)={float(best_final['R_init'][0]):.1f}")

# ── 2. Load 2D profile surfaces ──
def load_profile(tsv_name, x_col, y_col):
    df = pl.read_csv(OUT/tsv_name, separator="\t", comment_prefix="#")
    xs = np.array(sorted(set(df[x_col].to_list())))
    ys = np.array(sorted(set(df[y_col].to_list())))
    lls = np.full((len(xs), len(ys)), np.nan)
    for row in df.iter_rows(named=True):
        i = int(np.argmin(np.abs(xs - row[x_col])))
        j = int(np.argmin(np.abs(ys - row[y_col])))
        lls[i, j] = row["max_loglik"]
    return xs, ys, lls

bx, by, bll = load_profile("profile_true_bg_2d.tsv", "beta", "gamma")
ix, iy, ill = load_profile("profile_true_2d.tsv", "I0", "R_init")

def plot_panel(ax, xs, ys, lls, x_col, y_col, x_lab, y_lab, title):
    ll_max = np.nanmax(lls)
    delta = np.clip(ll_max - lls, 0, 30)
    im = ax.pcolormesh(xs, ys, delta.T, cmap="viridis_r", shading="auto",
                       vmin=0, vmax=30)
    # 2D Wilks contours (k=2 focal params)
    cs = ax.contour(xs, ys, (ll_max - lls).T, levels=[3.00, 4.61, 6.91],
                    colors=["cornflowerblue","#e67e22","#c0392b"], linewidths=1.0)
    ax.clabel(cs, fmt={3.00:"95%", 4.61:"99%", 6.91:"99.9%"}, fontsize=8)

    # All 128 scout chains — thinned light-grey trajectories
    for cid, grp in scout_thin.group_by("chain", maintain_order=True):
        if int(cid[0]) == best_chain: continue
        g = grp.sort("iteration")
        ax.plot(g[x_col].to_numpy(), g[y_col].to_numpy(),
                color="#808080", alpha=0.35, linewidth=0.5, zorder=3)
    ax.plot(best[x_col].to_numpy(), best[y_col].to_numpy(),
            color="#2c2c2c", alpha=0.95, linewidth=1.2, zorder=5,
            label=f"best scout (chain {best_chain})")
    # Best scout chain's endpoint — helps the eye see where that chain
    # actually lands vs where the profile max says the MLE is. The two
    # may differ because the scout chain and the profile max can sit on
    # different points of the same identifiability ridge.
    bx = float(best[x_col].to_numpy()[-1])
    by = float(best[y_col].to_numpy()[-1])
    ax.scatter([bx], [by], marker="o", s=60, facecolor="white",
               edgecolor="#2c2c2c", linewidth=1.0, zorder=6,
               label="chain 94 endpoint")

    # Profile max
    i_star, j_star = np.unravel_index(np.nanargmax(lls), lls.shape)
    ax.scatter([xs[i_star]], [ys[j_star]], marker="*", s=260, color="red",
               edgecolor="white", linewidth=1.0, zorder=7,
               label="profile max")

    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    ax.set_xlabel(x_lab); ax.set_ylabel(y_lab)
    ax.set_title(title, fontsize=10, pad=4)
    ax.spines[["top","right"]].set_visible(False)
    return im

fig, axes = plt.subplots(1, 2, figsize=(14, 5.6),
                         gridspec_kw={"wspace": 0.28})

im = plot_panel(axes[0], bx, by, bll, "beta", "gamma", "β", "γ",
                "(β, γ) profile — scout chains on top\n"
                "[I(0), R(0) IF2-optimised at each cell]")
axes[0].legend(frameon=False, fontsize=9, loc="upper right")

plot_panel(axes[1], ix, iy, ill, "I0", "R_init", "I(0)", "R(0)",
           "(I(0), R(0)) profile — scout chains on top\n"
           "[β, γ IF2-optimised at each cell]")

# Shared colorbar on the right
cbar_ax = fig.add_axes([0.915, 0.18, 0.012, 0.66])
fig.colorbar(im, cax=cbar_ax, label="Δ log-lik from profile max (clipped at 30)")

fig.suptitle("Poisson free-IC fit — scout chain trajectories overlaid on 2D profile surfaces",
             fontsize=11, y=0.99)
fig.subplots_adjust(top=0.88, bottom=0.13, left=0.06, right=0.90)
fig.savefig(FIG/"trajectory_poisson.png", bbox_inches="tight")
plt.close(fig)
print(f"wrote {FIG/'trajectory_poisson.png'}")
