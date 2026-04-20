"""Plot the 2D β × σ_se profile likelihood for the process-noise fit."""
import sys
from pathlib import Path
import numpy as np, polars as pl

HERE = Path(__file__).parent
BOOK = HERE.parent
sys.path.insert(0, str(BOOK))
from profile_plot import plot_profile_2d
from styles import MARKER_SQUARE

df = pl.read_csv(HERE/"output/profile_beta_sigma.tsv",
                 separator="\t", comment_prefix="#")

# Scout best-chain reference (from the earlier fit)
# β=2.746, σ_se=0.326 — put a square marker
ref = [(2.746, 0.326, MARKER_SQUARE, "scout best-chain")]

plot_profile_2d(
    df,
    focal_x="beta", focal_y="sigma_se",
    nuisance_params=["gamma"],
    labels={"beta": "β", "sigma_se": "σ_se", "gamma": "γ̂"},
    reference_points=ref,
    outpath=HERE/"figures/profile_beta_sigma.png",
    title="2D profile likelihood over (β, σ_se) — γ IF2-optimised per cell (process-noise fit, Poisson obs)",
    dll_clip=12,
)

betas  = np.array(sorted(set(df["beta"].to_list())))
sigmas = np.array(sorted(set(df["sigma_se"].to_list())))
lls = np.full((len(betas), len(sigmas)), np.nan)
gammas = np.full_like(lls, np.nan)
for row in df.iter_rows(named=True):
    i = int(np.argmin(np.abs(betas - row["beta"])))
    j = int(np.argmin(np.abs(sigmas - row["sigma_se"])))
    lls[i, j] = row["max_loglik"]
    gammas[i, j] = row["gamma"]

ll_max = np.nanmax(lls)
i_star, j_star = np.unravel_index(np.nanargmax(lls), lls.shape)
print(f"profile max ll = {ll_max:.2f} at β = {betas[i_star]:.3f}, "
      f"σ_se = {sigmas[j_star]:.3f} (γ̂ = {gammas[i_star, j_star]:.3f})")

# NegBin reference for context
NB_LL = -58.0857824
POI_LL = -61.1391851
print(f"\nFor comparison: NegBin ll = {NB_LL:.2f}  (Δ from profile max: {ll_max - NB_LL:.2f})")
print(f"                Poisson ll = {POI_LL:.2f}  (Δ from profile max: {ll_max - POI_LL:.2f})")

# Wilks fractions
for level, label in ((3.00, "95%"), (4.61, "99%"), (6.91, "99.9%")):
    frac = float(((ll_max - lls) < level).sum()) / lls.size
    print(f"  grid fraction within Δ {level:.2f} ({label} 2D CI): {100*frac:.0f}%")
