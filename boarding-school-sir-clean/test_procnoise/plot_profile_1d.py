"""Plot the 1D σ_se profile likelihood with Wilks 1D CI."""
import sys
from pathlib import Path
import numpy as np, polars as pl
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
BOOK = HERE.parent
sys.path.insert(0, str(BOOK))
from styles import (apply_rc, GREY, CI_LEVELS_1D, CI_CONTOUR_COLORS,
                    CI_LABELS_BY_LEVEL)

apply_rc()
df = pl.read_csv(HERE/"output/profile_sigma_1d.tsv",
                 separator="\t", comment_prefix="#").sort("sigma_se")

sigmas = df["sigma_se"].to_numpy()
lls    = df["max_loglik"].to_numpy()
betas  = df["beta"].to_numpy()
gammas = df["gamma"].to_numpy()

ll_max = lls.max()
i_star = int(np.argmax(lls))

NB_LL  = -58.0857824
POI_LL = -61.1391851

fig, axes = plt.subplots(1, 3, figsize=(13, 3.8),
                         gridspec_kw={"wspace": 0.32})

# Panel 1: 1D profile with Wilks CI
ax = axes[0]
ax.plot(sigmas, lls, "-o", color="#2c3e50", linewidth=1.5,
        markerfacecolor="white", markeredgecolor="#2c3e50",
        markeredgewidth=1.2, markersize=6)
ax.scatter([sigmas[i_star]], [lls[i_star]], marker="*", s=260,
           color="red", edgecolor="white", linewidth=1.0, zorder=6,
           label=f"profile max (σ_se={sigmas[i_star]:.2f})")
for lvl, col in zip(CI_LEVELS_1D, CI_CONTOUR_COLORS):
    ax.axhline(ll_max - lvl, color=col, linestyle=":", linewidth=1.0,
               label=f"{CI_LABELS_BY_LEVEL.get(lvl, str(lvl))} Wilks CI")
ax.axhline(NB_LL, color="#3498db", linestyle="--", linewidth=1.0,
           alpha=0.9, label=f"NegBin fixed-IC ll = {NB_LL:.2f}")
ax.axhline(POI_LL, color="#95a5a6", linestyle="--", linewidth=1.0,
           alpha=0.9, label=f"Poisson fixed-IC ll = {POI_LL:.2f}")
ax.set_xscale("log")
ax.set_xlabel("σ_se (process-noise SD on β)")
ax.set_ylabel("profile log-likelihood")
ax.set_title("1D profile over σ_se  (β, γ optimised per cell)",
             fontsize=10, pad=4)
ax.legend(frameon=False, fontsize=7.5, loc="lower left")
ax.spines[["top","right"]].set_visible(False)

# Panel 2: β̂ along the profile
ax = axes[1]
ax.plot(sigmas, betas, "-o", color="#2c3e50", linewidth=1.5,
        markerfacecolor="white", markeredgecolor="#2c3e50",
        markeredgewidth=1.2, markersize=6)
ax.scatter([sigmas[i_star]], [betas[i_star]], marker="*", s=260,
           color="red", edgecolor="white", linewidth=1.0, zorder=6)
ax.set_xscale("log")
ax.set_xlabel("σ_se")
ax.set_ylabel("β̂  (MLE at each σ_se)")
ax.set_title("β̂ traces the ridge", fontsize=10, pad=4)
ax.spines[["top","right"]].set_visible(False)

# Panel 3: γ̂ along the profile
ax = axes[2]
ax.plot(sigmas, gammas, "-o", color="#2c3e50", linewidth=1.5,
        markerfacecolor="white", markeredgecolor="#2c3e50",
        markeredgewidth=1.2, markersize=6)
ax.scatter([sigmas[i_star]], [gammas[i_star]], marker="*", s=260,
           color="red", edgecolor="white", linewidth=1.0, zorder=6)
ax.set_xscale("log")
ax.set_xlabel("σ_se")
ax.set_ylabel("γ̂")
ax.set_title("γ̂ is stable (well-identified)", fontsize=10, pad=4)
ax.spines[["top","right"]].set_visible(False)

fig.suptitle("Process-noise fit: σ_se 1D profile likelihood "
             "(camdl profile, IF2 per cell)",
             fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(HERE/"figures/profile_sigma_1d.png", bbox_inches="tight")
plt.close(fig)

# Print summary stats
print(f"profile max ll = {ll_max:.2f} at σ_se = {sigmas[i_star]:.3f}  "
      f"(β = {betas[i_star]:.3f}, γ = {gammas[i_star]:.3f})")
mask_95 = (ll_max - lls) < CI_LEVELS_1D[0]
if mask_95.any():
    print(f"95% Wilks CI on σ_se (Δll < {CI_LEVELS_1D[0]}): "
          f"[{sigmas[mask_95].min():.3f}, {sigmas[mask_95].max():.3f}]")
print(f"Δll to NegBin = {ll_max - NB_LL:.2f}")
print(f"Δll to Poisson = {ll_max - POI_LL:.2f}")
