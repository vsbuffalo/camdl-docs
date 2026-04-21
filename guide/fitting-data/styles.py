"""Shared plot styling — extract-ready for a future `camdl-diag` package.

Centralizes colormaps, confidence-level thresholds, marker styles, and
matplotlib rcParams used across the chapter's diagnostic plots.
"""
import matplotlib as mpl

GREY = "#555555"

# Colormap conventions — different families for Δll vs nuisance-MLE so
# the two are never visually confused.
CMAP_DLL     = "viridis_r"   # Δ log-lik from profile max (dark = bad)
CMAP_NUIS    = "cividis"     # nuisance MLE at each grid cell

# Wilks χ² thresholds for likelihood-ratio confidence regions.
# 1D (one focal param):
CI_LEVELS_1D = [1.92, 4.61, 9.21]   # 95 / 99 / 99.9 %
# 2D (two focal params):
CI_LEVELS_2D = [3.00, 4.61, 6.91]   # 95 / 99 / 99.9 %
CI_CONTOUR_COLORS = ["cornflowerblue", "#e67e22", "#c0392b"]
CI_LABELS_BY_LEVEL = {
    1.92: "95%", 4.61: "99%", 9.21: "99.9%",
    3.00: "95%", 6.91: "99.9%",
}

# Marker styles for profile-max / baseline / grid-max points
MARKER_STAR   = dict(marker="*", s=260, color="red",
                     edgecolor="white", linewidth=1.0, zorder=5)
MARKER_SQUARE = dict(marker="s", s=90, color="#3498db",
                     edgecolor="white", linewidth=0.8, zorder=5)
MARKER_DIAMND = dict(marker="D", s=110, color="#2ecc71",
                     edgecolor="white", linewidth=1.0, zorder=5)


def apply_rc():
    mpl.rcParams.update({
        "figure.dpi": 150, "font.size": 10,
        "axes.edgecolor": GREY, "axes.labelcolor": GREY,
        "xtick.color": GREY, "ytick.color": GREY, "text.color": GREY,
    })
