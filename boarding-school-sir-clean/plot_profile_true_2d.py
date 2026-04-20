"""Plot the 2D profile log-lik over (I(0), R(0)) from camdl profile.
True profile — β and γ IF2-optimized at each grid cell."""
from pathlib import Path
import polars as pl, numpy as np

from styles import MARKER_SQUARE
from profile_plot import plot_profile_2d

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"

df = pl.read_csv(OUT/"profile_true_2d.tsv", separator="\t", comment_prefix="#")

# Report profile max + reference points before plotting
I0s = np.array(sorted(set(df["I0"].to_list())))
R0s = np.array(sorted(set(df["R_init"].to_list())))
lls = np.full((len(I0s), len(R0s)), np.nan)
for row in df.iter_rows(named=True):
    i = int(np.argmin(np.abs(I0s - row["I0"])))
    j = int(np.argmin(np.abs(R0s - row["R_init"])))
    lls[i, j] = row["max_loglik"]
ll_max = np.nanmax(lls)
i_star, j_star = np.unravel_index(np.nanargmax(lls), lls.shape)
print(f"profile max ll = {ll_max:.2f} at I(0) = {I0s[i_star]:.0f}, "
      f"R(0) = {R0s[j_star]:.0f}")

plot_profile_2d(
    df,
    focal_x="I0", focal_y="R_init",
    nuisance_params=["beta", "gamma"],
    labels={"I0": "I(0)", "R_init": "R(0)", "beta": "β̂", "gamma": "γ̂"},
    reference_points=[(5, 0, MARKER_SQUARE, "fixed-IC baseline")],
    outpath=FIG/"profile_true_2d.png",
    title="True 2D profile likelihood over (I(0), R(0)) — camdl profile (IF2 per cell)",
)
print(f"wrote {FIG/'profile_true_2d.png'}")

for level, label in ((3.00, "95%"), (4.61, "99%"), (6.91, "99.9%")):
    frac = float(((ll_max - lls) < level).sum()) / lls.size
    print(f"  grid fraction within Δ {level:.2f} ({label} 2D CI): {100*frac:.0f}%")
