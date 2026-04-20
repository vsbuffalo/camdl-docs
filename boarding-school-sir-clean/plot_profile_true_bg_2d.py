"""2D profile on (β, γ) with I(0), R(0) nuisance-optimized at each cell."""
from pathlib import Path
import polars as pl, numpy as np

from styles import MARKER_SQUARE
from profile_plot import plot_profile_2d

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"

BETA_FIXED_IC  = 1.9058652207
GAMMA_FIXED_IC = 0.6559714166

df = pl.read_csv(OUT/"profile_true_bg_2d.tsv", separator="\t", comment_prefix="#")

betas  = np.array(sorted(set(df["beta"].to_list())))
gammas = np.array(sorted(set(df["gamma"].to_list())))
lls = np.full((len(betas), len(gammas)), np.nan)
for row in df.iter_rows(named=True):
    i = int(np.argmin(np.abs(betas - row["beta"])))
    j = int(np.argmin(np.abs(gammas - row["gamma"])))
    lls[i, j] = row["max_loglik"]
ll_max = np.nanmax(lls)
i_star, j_star = np.unravel_index(np.nanargmax(lls), lls.shape)
print(f"profile max ll = {ll_max:.2f} at β = {betas[i_star]:.3f}, "
      f"γ = {gammas[j_star]:.3f}")

plot_profile_2d(
    df,
    focal_x="beta", focal_y="gamma",
    nuisance_params=["I0", "R_init"],
    labels={"beta": "β", "gamma": "γ", "I0": "Î(0)", "R_init": "R̂(0)"},
    reference_points=[(BETA_FIXED_IC, GAMMA_FIXED_IC, MARKER_SQUARE,
                       "fixed-IC baseline")],
    outpath=FIG/"profile_true_bg_2d.png",
    title="True 2D profile likelihood over (β, γ) — camdl profile (IF2 per cell)",
)
print(f"wrote {FIG/'profile_true_bg_2d.png'}")

for level, label in ((3.00, "95%"), (4.61, "99%"), (6.91, "99.9%")):
    frac = float(((ll_max - lls) < level).sum()) / lls.size
    print(f"  grid fraction within Δ {level:.2f} ({label} 2D CI): {100*frac:.0f}%")
