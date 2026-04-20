"""2D β × σ_se profile Δll heatmap with scout chain endpoints overlaid.
Shows where the 64 scout chains actually landed vs where the true ridge is."""
import sys
from pathlib import Path
import numpy as np, polars as pl
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
BOOK = HERE.parent
sys.path.insert(0, str(BOOK))
from profile_plot import plot_dll_heatmap
from plot_traces import load_chains, best_chain_by_final_ll
from styles import apply_rc, MARKER_STAR

apply_rc()

# --- Profile surface
df = pl.read_csv(HERE/"output/profile_beta_sigma.tsv",
                 separator="\t", comment_prefix="#")
betas  = np.array(sorted(set(df["beta"].to_list())))
sigmas = np.array(sorted(set(df["sigma_se"].to_list())))
lls = np.full((len(betas), len(sigmas)), np.nan)
for row in df.iter_rows(named=True):
    i = int(np.argmin(np.abs(betas  - row["beta"])))
    j = int(np.argmin(np.abs(sigmas - row["sigma_se"])))
    lls[i, j] = row["max_loglik"]

ll_max = np.nanmax(lls)
i_star, j_star = np.unravel_index(np.nanargmax(lls), lls.shape)

# --- Scout chains
SCOUT = next(HERE.glob("results/fits/fit_procnoise-*/real/fit_42/scout"))
scout = load_chains(SCOUT).sort(["chain", "iteration"])
best_c = best_chain_by_final_ll(scout)

# Endpoint per chain: last row of each chain
endpoints = (scout.group_by("chain", maintain_order=True)
             .agg(pl.col("beta").last().alias("beta_end"),
                  pl.col("sigma_se").last().alias("sigma_end"),
                  pl.col("loglik").drop_nulls().last().alias("ll_end")))
b_end = endpoints["beta_end"].to_numpy()
s_end = endpoints["sigma_end"].to_numpy()
ll_end = endpoints["ll_end"].to_numpy()
cids   = endpoints["chain"].to_numpy()

# Best chain's trajectory (for the line overlay)
best_df = scout.filter(pl.col("chain") == best_c).sort("iteration")
b_traj = best_df["beta"].to_numpy()
s_traj = best_df["sigma_se"].to_numpy()

# --- Plot: single zoomed panel on the ridge
fig, ax_single = plt.subplots(1, 1, figsize=(8.5, 5.0))
axes = np.array([None, ax_single])  # keep zoom panel at index 1

def draw_panel(ax, bslice=None, sslice=None, show_traj=True):
    if bslice is None:
        bs, bi0, bi1 = betas, 0, len(betas)
    else:
        bi0, bi1 = bslice
        bs = betas[bi0:bi1]
    if sslice is None:
        ss, si0, si1 = sigmas, 0, len(sigmas)
    else:
        si0, si1 = sslice
        ss = sigmas[si0:si1]
    Z = lls[bi0:bi1, si0:si1]

    im = plot_dll_heatmap(
        ax, bs, ss, Z, dll_clip=15,
        xlabel="β", ylabel="σ_se",
        markers=[(betas[i_star], sigmas[j_star], MARKER_STAR, "profile max")],
    )

    # Scout endpoints — plain gray, standard non-highlight chain color
    sc = ax.scatter(b_end, s_end, color="#999999",
                    s=40, edgecolor="white", linewidth=0.5, zorder=4,
                    label=f"scout chain endpoints (n={len(cids)})")
    # Best chain trajectory — dark gray to match other chain diagnostics
    if show_traj:
        ax.plot(b_traj, s_traj, "-", color="#2c3e50",
                linewidth=1.2, alpha=0.9, zorder=3.5,
                label=f"best chain #{best_c} trajectory")
        ax.scatter([b_traj[0]], [s_traj[0]], marker="o", s=40,
                   facecolor="white", edgecolor="#2c3e50",
                   linewidth=1.2, zorder=5, label="chain start")

    ax.legend(frameon=False, fontsize=8,
              loc="center left", bbox_to_anchor=(1.18, 0.5))
    return im, sc

# Zoom: find the ridge region
zi, zj = np.where((ll_max - lls) < 4.0)
bi0, bi1 = max(0, zi.min()-1), min(len(betas), zi.max()+2)
si0, si1 = max(0, zj.min()-1), min(len(sigmas), zj.max()+2)
im2, sc2 = draw_panel(axes[1], bslice=(bi0, bi1), sslice=(si0, si1))
axes[1].set_title(f"β × σ_se profile likelihood (ridge zoom, Δll < 4) "
                  f"with 64 scout chain endpoints",
                  fontsize=10.5, pad=6)
fig.colorbar(im2, ax=axes[1], label="Δ log-lik from profile max (clipped at 15)",
             pad=0.02)
fig.savefig(HERE/"figures/profile_beta_sigma_with_chains.png",
            bbox_inches="tight")
plt.close(fig)

# Print diagnostics
print(f"profile max: β = {betas[i_star]:.3f}, σ_se = {sigmas[j_star]:.3f}, "
      f"ll = {ll_max:.2f}")
print(f"best scout chain: #{best_c}")
print(f"  endpoint: β = {b_end[cids == best_c][0]:.3f}, "
      f"σ_se = {s_end[cids == best_c][0]:.3f}, "
      f"ll = {ll_end[cids == best_c][0]:.2f}")
dll_end = ll_max - ll_end
print(f"\nScout endpoint Δll distribution (smaller is closer to profile max):")
for q in (0.1, 0.25, 0.5, 0.75, 0.9):
    print(f"  q{int(q*100):2d}: {np.quantile(dll_end, q):5.2f}")
print(f"\nEndpoints within Δll < 1 of profile max: "
      f"{int((dll_end < 1).sum())}/{len(dll_end)}")
print(f"Endpoints within Δll < 2 of profile max: "
      f"{int((dll_end < 2).sum())}/{len(dll_end)}")
print(f"Endpoints within Δll < 5 of profile max: "
      f"{int((dll_end < 5).sum())}/{len(dll_end)}")
