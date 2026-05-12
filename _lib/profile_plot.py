"""Profile-likelihood figures for camdl-book chapters.

Two helpers, one per profile shape:

  - ``plot_profile_panels`` — stacked 1D profile (one swept parameter
    on x, profile log-lik + nuisance MLEs stacked vertically).

  - ``plot_profile_2d`` — multi-panel 2D profile (two swept
    parameters on a grid; Δ log-lik heatmap with Wilks contours plus
    one nuisance-MLE heatmap per nuisance parameter).

Both consume the long-form TSV output of ``camdl profile``. Per the
post-`b8878e2` schema:

  - focal swept param(s) (e.g. ``s0``, ``R0``, ``beta``)
  - ``loglik`` — central log-likelihood at each grid cell (always
    present; equals the per-seed value when n_seeds=1, the mean
    across seeds when n_seeds>1)
  - one bare nuisance column per IF2-optimised parameter (``R0``,
    ``alpha``, ``gamma``, etc.) — the conditional MLE at that grid
    cell

When ``--seeds 1:N`` runs (n_seeds>1), additional spread-diagnostic
columns are appended: ``loglik_sd``, ``loglik_min``, ``loglik_max``
for the loglik, plus ``<param>_sd`` paired with each nuisance.
``loglik_sd`` is the per-cell trustworthiness diagnostic — high
spread → that cell's MLE is not pinned across seeds.

Style conventions (centralised here so 1D and 2D plots match):
  - Δll uses ``viridis_r`` (dark = bad, light = good).
  - Nuisance MLE heatmaps use ``cividis`` so they're never visually
    confused with Δll.
  - Wilks 2D contours at Δll = 3.00, 4.61, 6.91 (95 / 99 / 99.9 %).
  - Profile-max marked with a red star; user-supplied reference
    points (e.g. fixed-IC MLE for cross-comparison) get their own
    style via the ``reference_points`` kwarg.
"""
from __future__ import annotations
import numpy as np
import polars as pl
import matplotlib.pyplot as plt

PROFILE_COLOR = "#444444"
TRUTH_COLOR   = "black"
MLE_COLOR     = "#e67e22"

# 2D-profile colormap conventions
CMAP_DLL  = "viridis_r"
CMAP_NUIS = "cividis"

# Wilks χ² thresholds for likelihood-ratio confidence regions.
CI_LEVELS_2D       = [3.00, 4.61, 6.91]   # 95 / 99 / 99.9 %
CI_CONTOUR_COLORS  = ["cornflowerblue", "#e67e22", "#c0392b"]
CI_LABELS_BY_LEVEL = {3.00: "95%", 4.61: "99%", 6.91: "99.9%"}

# Marker styles
MARKER_STAR   = dict(marker="*", s=260, color="red",
                     edgecolor="white", linewidth=1.0, zorder=5)
MARKER_SQUARE = dict(marker="s", s=90,  color="#3498db",
                     edgecolor="white", linewidth=0.8, zorder=5)
MARKER_DIAMND = dict(marker="D", s=110, color="#2ecc71",
                     edgecolor="white", linewidth=1.0, zorder=5)


def _pivot_grid(df, x_col, y_col, v_col):
    """Long-form profile DataFrame → (xs, ys, matrix[i, j])."""
    xs = np.array(sorted(set(df[x_col].to_list())))
    ys = np.array(sorted(set(df[y_col].to_list())))
    Z  = np.full((len(xs), len(ys)), np.nan)
    for row in df.iter_rows(named=True):
        i = int(np.argmin(np.abs(xs - row[x_col])))
        j = int(np.argmin(np.abs(ys - row[y_col])))
        Z[i, j] = row[v_col]
    return xs, ys, Z


def _place_markers(ax, markers):
    """markers: list of (x, y, style_dict, label_or_None)."""
    for x, y, style, label in markers:
        kw = dict(style)
        if label:
            kw["label"] = label
        ax.scatter([x], [y], **kw)


def plot_dll_heatmap(ax, xs, ys, lls, *, dll_clip=12,
                     contour_levels=None, contour_colors=None,
                     markers=None, xlabel="", ylabel="", title=""):
    """Single panel: Δ log-lik heatmap (clipped) + Wilks contours."""
    contour_levels = contour_levels or CI_LEVELS_2D
    contour_colors = contour_colors or CI_CONTOUR_COLORS
    ll_max = np.nanmax(lls)
    delta  = np.clip(ll_max - lls, 0, dll_clip)
    im = ax.pcolormesh(xs, ys, delta.T, cmap=CMAP_DLL, shading="auto",
                       vmin=0, vmax=dll_clip)
    cs = ax.contour(xs, ys, (ll_max - lls).T, levels=contour_levels,
                    colors=contour_colors, linewidths=1.0)
    ax.clabel(cs, fmt={lvl: CI_LABELS_BY_LEVEL.get(lvl, f"{lvl:.1f}")
                       for lvl in contour_levels}, fontsize=8)
    if markers:
        _place_markers(ax, markers)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=10, pad=4)
    ax.spines[["top", "right"]].set_visible(False)
    return im


def plot_nuisance_heatmap(ax, xs, ys, Z, *, cmap=None, markers=None,
                          xlabel="", ylabel="", title=""):
    """Single panel: nuisance-MLE heatmap at each (xs, ys) grid cell."""
    im = ax.pcolormesh(xs, ys, Z.T, cmap=cmap or CMAP_NUIS, shading="auto")
    if markers:
        _place_markers(ax, markers)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=10, pad=4)
    ax.spines[["top", "right"]].set_visible(False)
    return im


def plot_profile_2d(df, *, focal_x, focal_y, nuisance_params,
                    labels=None, reference_points=None,
                    title="", figsize=None, dll_clip=12,
                    ll_col="loglik"):
    """Multi-panel 2D profile figure.

    Layout: Δll heatmap + one nuisance-MLE heatmap per ``nuisance_params``
    entry. Auto-laid out in up to 3 columns.

    Args:
        df: long-form profile DataFrame (polars). Must contain
            ``focal_x``, ``focal_y``, ``ll_col``, and one column per
            entry of ``nuisance_params``.
        focal_x, focal_y: str — names of the two swept parameters.
        nuisance_params: list[str] — nuisance MLE columns to render.
        labels: optional dict mapping column name → display label.
        reference_points: optional list of ``(x, y, style, label)`` tuples
            drawn on every panel (e.g. fixed-IC MLE marker).
        title: optional figure-level suptitle.
        figsize: override default ``(5.5 * n_cols, 4.8 * n_rows)``.
        dll_clip: colormap clipped to this many nats from profile max.
        ll_col: log-likelihood column name (default ``loglik``).

    Returns:
        ``(fig, axes_array)`` — caller can call ``fig.tight_layout()``
        and either save or render. Profile-max coordinates are
        annotated on the figure via the legend on the Δll panel.
    """
    labels = labels or {}
    def lbl(col): return labels.get(col, col)

    xs, ys, lls = _pivot_grid(df, focal_x, focal_y, ll_col)
    nuis_grids  = {p: _pivot_grid(df, focal_x, focal_y, p)[2]
                   for p in nuisance_params}

    ll_max = np.nanmax(lls)
    i_, j_ = np.unravel_index(np.nanargmax(lls), lls.shape)
    fmax_x = float(xs[i_]); fmax_y = float(ys[j_])
    # Short label — values belong in the figure caption, not legend
    star_label = "profile max"
    star = (fmax_x, fmax_y, MARKER_STAR, star_label)
    markers = [star] + list(reference_points or [])
    # bare-marker variant for nuisance panels: no labels (legend only on Δll)
    bare_markers = [(x, y, st, None) for (x, y, st, _) in markers]

    n_panels = 1 + len(nuisance_params)
    n_cols   = min(3, n_panels)
    n_rows   = (n_panels + n_cols - 1) // n_cols
    if figsize is None:
        # 5.0 wide × 4.0 tall per panel — for 3 panels in a row that's
        # a 15 × 4.0 figure (aspect 3.75:1), close to v1's ~2.85:1
        # published figures while a touch more landscape because we
        # also render the in-panel legend.
        figsize = (5.0 * n_cols, 4.0 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize,
                             gridspec_kw={"wspace": 0.55, "hspace": 0.32})
    axes = np.atleast_1d(axes).flatten()

    # Panel 1: Δll. Right-margin layout — colorbar in lower 40% of the
    # right strip, legend in upper 50%. Colorbar label as a small
    # title above the bar (using set_title on the cax) so the label
    # doesn't rotate-extend up into the legend area, which was the
    # failure mode of the default `colorbar(label=…)` placement.
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    im0 = plot_dll_heatmap(
        axes[0], xs, ys, lls, dll_clip=dll_clip, markers=markers,
        xlabel=lbl(focal_x), ylabel=lbl(focal_y),
        title=f"Δ log-likelihood   (max ll = {ll_max:.1f})")
    cax0 = inset_axes(axes[0], width="4%", height="40%",
                      loc="lower left",
                      bbox_to_anchor=(1.04, 0.02, 1, 1),
                      bbox_transform=axes[0].transAxes,
                      borderpad=0)
    cb0 = fig.colorbar(im0, cax=cax0)
    cb0.ax.tick_params(labelsize=7)
    cax0.set_title(f"Δ ll\n(≤ {dll_clip})", fontsize=7, pad=3, loc="left")
    # Legend in the upper 50% of the right strip.
    axes[0].legend(
        loc="upper left", fontsize=7,
        bbox_to_anchor=(1.02, 1.0), bbox_transform=axes[0].transAxes,
        frameon=True, facecolor="white", edgecolor="#888888",
        framealpha=1.0,
    )

    # Nuisance panels — same shrunk-colorbar layout for visual consistency.
    for ax, p in zip(axes[1:], nuisance_params):
        Z = nuis_grids[p]
        im = plot_nuisance_heatmap(
            ax, xs, ys, Z, markers=bare_markers,
            xlabel=lbl(focal_x), ylabel=lbl(focal_y),
            title=f"{lbl(p)} — IF2-optimised per cell")
        cax = inset_axes(ax, width="4%", height="40%",
                         loc="lower left",
                         bbox_to_anchor=(1.04, 0.02, 1, 1),
                         bbox_transform=ax.transAxes,
                         borderpad=0)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=7)
        cax.set_title(lbl(p), fontsize=7, pad=3, loc="left")

    for ax in axes[n_panels:]:
        ax.axis("off")
    if title:
        fig.suptitle(title, fontsize=11, y=0.99)
    return fig, axes


def plot_profile_panels(
    p: pl.DataFrame,
    x_col: str,
    panel_params: list[str],
    panel_labels: list[str],
    truth: dict[str, float],
    truth_ll: float,
    x_log: bool = True,
    x_range: tuple[float, float] | None = None,
    title: str = "",
    color: str = PROFILE_COLOR,
    ll_col: str = "loglik",
    show_truth_ll: bool = False,
):
    """Stacked-panel profile figure.

    Args:
      p: profile output dataframe (columns: x_col, ll_col, *panel_params).
      x_col: name of the swept parameter (e.g. "s0", "R0").
      panel_params: list whose first entry is ``ll_col`` (the loglik panel)
        and remaining entries are trading-off conditional-MLE panels (e.g.
        ``"R0"``, ``"alpha"``, ``"gamma"``).
      panel_labels: matching y-axis labels.
      truth: dict of param → truth value (must include x_col + every
        non-loglik param in panel_params).
      truth_ll: truth log-likelihood (for the top panel's horizontal line).
      x_log: log scale on x (use True for full-range views, False for
        narrow zooms).
      x_range: optional (xmin, xmax) — if given, filter p to this range
        before plotting.
      title: figure-level title.
      color: line/marker color.
      ll_col: log-likelihood column name (default ``loglik``; the
        first entry of ``panel_params`` is matched against this to
        identify which panel renders the loglik surface).
      show_truth_ll: if True, draw a horizontal dashed line at
        ``truth_ll`` on the loglik panel. Default False — drawing this
        line forces matplotlib to expand the y-axis to include
        ``truth_ll``, which can flatten the profile curvature
        visually when the truth ll sits far above (or below) the
        profile-MLE neighbourhood.
    """
    if x_range is not None:
        p = p.filter((pl.col(x_col) >= x_range[0]) & (pl.col(x_col) <= x_range[1]))

    # Profile MLE = grid point with max ll.
    ll = p[ll_col].to_numpy()
    mle_idx = int(np.argmax(ll))
    x_mle  = float(p[x_col][mle_idx])
    ll_mle = float(ll[mle_idx])

    fig, axes = plt.subplots(len(panel_params), 1, figsize=(7, 8.5),
                              sharex=True)
    if len(panel_params) == 1:
        axes = np.array([axes])

    for ax, col, lab in zip(axes, panel_params, panel_labels):
        ax.plot(p[x_col], p[col], color=color, linewidth=1.5,
                marker="o", markersize=6)

        if col == ll_col:
            if show_truth_ll:
                ax.axhline(truth_ll, color=TRUTH_COLOR, linestyle="--",
                           linewidth=1, label=f"truth ℓ = {truth_ll:.0f}")
            ax.scatter([x_mle], [ll_mle], marker="*", color=MLE_COLOR,
                       s=200, zorder=10, edgecolor="black", linewidth=0.6,
                       label=f"profile MLE ({x_col}={x_mle:.3g}, ℓ={ll_mle:.0f})")
            ax.legend(frameon=False, fontsize=8, loc="center left",
                      bbox_to_anchor=(1.02, 0.5))
        elif col in truth:
            ax.axhline(truth[col], color=TRUTH_COLOR, linestyle="--",
                       linewidth=1, label=f"truth {col}={truth[col]:.3g}")
            ax.legend(frameon=False, fontsize=8, loc="center left",
                      bbox_to_anchor=(1.02, 0.5))

        if x_col in truth:
            ax.axvline(truth[x_col], color=TRUTH_COLOR, linestyle=":",
                       linewidth=1)
        ax.axvline(x_mle, color=MLE_COLOR, linestyle="-.", linewidth=1)
        ax.set_ylabel(lab, fontsize=10)

    axes[-1].set_xlabel(f"{x_col} (profile)" + (", zoomed" if x_range else ""))
    if x_log:
        axes[-1].set_xscale("log")
    if title:
        axes[0].set_title(title)
    fig.tight_layout()
    return fig, axes
