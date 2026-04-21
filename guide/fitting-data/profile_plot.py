"""Plotting helpers for `camdl profile` output — extract-ready for the
future `camdl-diag` package.

Functions:
  - `plot_profile_2d(df, focal_x, focal_y, nuisance_params, ...)` —
    the main entry: takes a 2D-profile TSV DataFrame and produces the
    canonical figure (Δll heatmap + one nuisance-MLE heatmap per
    nuisance param).
  - `plot_nuisance_heatmap(ax, xs, ys, Z, ...)` — single panel,
    nuisance-MLE surface with consistent cividis colormap.
  - `plot_dll_heatmap(ax, xs, ys, lls, ...)` — single panel,
    Δll surface with 2D Wilks contours.
"""
from pathlib import Path
import numpy as np, polars as pl
import matplotlib.pyplot as plt

from styles import (GREY, CMAP_DLL, CMAP_NUIS,
                    CI_LEVELS_2D, CI_CONTOUR_COLORS, CI_LABELS_BY_LEVEL,
                    MARKER_STAR, MARKER_SQUARE, apply_rc)


def _pivot_grid(df, x_col, y_col, v_col):
    """Turn long-form profile DataFrame into (xs, ys, matrix[i,j])."""
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
        if label: kw["label"] = label
        ax.scatter([x], [y], **kw)


def plot_dll_heatmap(ax, xs, ys, lls, *,
                    dll_clip=30, contour_levels=CI_LEVELS_2D,
                    contour_colors=CI_CONTOUR_COLORS,
                    markers=None, xlabel="", ylabel="", title=""):
    """Panel: Δ log-likelihood heatmap (clipped) + Wilks contours.

    Args:
        ax: matplotlib Axes to draw on.
        xs, ys: 1D coordinate arrays (x on x-axis, y on y-axis).
        lls: 2D array of log-lik values, shape (len(xs), len(ys)).
        dll_clip: colormap clipped to [0, dll_clip] nats from max.
        contour_levels: Δll levels for overlaid contours.
        markers: list of (x, y, style_dict, label_or_None).
        xlabel, ylabel, title: panel labels.
    Returns: the QuadMesh for reuse in a shared colorbar.
    """
    ll_max = np.nanmax(lls)
    delta  = np.clip(ll_max - lls, 0, dll_clip)
    im = ax.pcolormesh(xs, ys, delta.T, cmap=CMAP_DLL, shading="auto",
                       vmin=0, vmax=dll_clip)
    cs = ax.contour(xs, ys, (ll_max - lls).T, levels=contour_levels,
                    colors=contour_colors, linewidths=1.0)
    labels = {lvl: CI_LABELS_BY_LEVEL.get(lvl, f"{lvl:.1f}")
              for lvl in contour_levels}
    ax.clabel(cs, fmt=labels, fontsize=8)
    if markers: _place_markers(ax, markers)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if title: ax.set_title(title, fontsize=10, pad=4)
    ax.spines[["top","right"]].set_visible(False)
    return im


def plot_nuisance_heatmap(ax, xs, ys, Z, *,
                          cmap=CMAP_NUIS, markers=None,
                          xlabel="", ylabel="", title="", cbar_label=""):
    """Panel: nuisance-MLE heatmap at each (xs, ys) grid cell.

    Uses a colormap distinct from Δll plots so the two panel types are
    never visually confused (viridis_r for Δll, cividis for nuisance).

    Returns: the QuadMesh.
    """
    im = ax.pcolormesh(xs, ys, Z.T, cmap=cmap, shading="auto")
    if markers: _place_markers(ax, markers)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if title: ax.set_title(title, fontsize=10, pad=4)
    ax.spines[["top","right"]].set_visible(False)
    return im


def plot_profile_2d(df, *, focal_x, focal_y, nuisance_params,
                     labels=None, outpath, title="",
                     reference_points=None, fig_figsize=None,
                     dll_clip=30):
    """Canonical 2D-profile figure: Δll panel + one nuisance-MLE panel
    per entry of `nuisance_params`. Auto-laid-out in a 2×N_cols grid.

    Args:
        df: polars DataFrame from `camdl profile` output; must have
            columns `focal_x`, `focal_y`, `max_loglik`, and each nuisance
            param.
        focal_x, focal_y: column names of focal (swept) params.
        nuisance_params: list of column names to plot as MLE heatmaps.
        labels: optional dict {col: "display label"} — if absent,
            col names are used.
        outpath: output path.
        reference_points: list of (x, y, style, label) tuples; the
            "profile max" marker is added automatically on top.
        fig_figsize: override figure size.
        dll_clip: Δll colormap clip.
    """
    apply_rc()
    labels = labels or {}
    def lbl(col): return labels.get(col, col)

    # Build Δll grid + per-nuisance grids
    xs, ys, lls = _pivot_grid(df, focal_x, focal_y, "max_loglik")
    nuis_grids = {p: _pivot_grid(df, focal_x, focal_y, p)[2]
                  for p in nuisance_params}

    ll_max = np.nanmax(lls)
    i_star, j_star = np.unravel_index(np.nanargmax(lls), lls.shape)
    focal_mle_x = float(xs[i_star])
    focal_mle_y = float(ys[j_star])

    # Build marker list: profile max star + user refs
    star_marker = (focal_mle_x, focal_mle_y, MARKER_STAR, "profile max")
    markers = [star_marker] + list(reference_points or [])

    # Layout: 1 row if only Δll, else 2 rows if 1-2 nuisance, else 2 cols
    n_panels = 1 + len(nuisance_params)
    n_cols = min(3, n_panels)
    n_rows = (n_panels + n_cols - 1) // n_cols
    if fig_figsize is None:
        fig_figsize = (5.5 * n_cols, 4.8 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=fig_figsize,
                              gridspec_kw={"wspace": 0.28, "hspace": 0.26})
    axes = np.atleast_1d(axes).flatten()

    # Panel 1: Δll
    im_dll = plot_dll_heatmap(
        axes[0], xs, ys, lls,
        dll_clip=dll_clip, markers=markers,
        xlabel=lbl(focal_x), ylabel=lbl(focal_y),
        title="Δ log-likelihood")
    # Only show legend on the Δll panel (first)
    if any(m[3] for m in markers):
        axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    fig.colorbar(im_dll, ax=axes[0],
                  label=f"Δ log-lik from profile max (clipped at {dll_clip})")

    # Nuisance panels
    for ax, p in zip(axes[1:], nuisance_params):
        Z = nuis_grids[p]
        im = plot_nuisance_heatmap(
            ax, xs, ys, Z, markers=markers,
            xlabel=lbl(focal_x), ylabel=lbl(focal_y),
            title=f"{lbl(p)} (IF2-optimised at each cell)",
            cbar_label=lbl(p))
        fig.colorbar(im, ax=ax, label=lbl(p))

    # hide any spare axes
    for ax in axes[n_panels:]:
        ax.axis("off")

    if title:
        fig.suptitle(title, fontsize=11, y=0.99)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath
