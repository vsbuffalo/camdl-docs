"""Pair-plot of IF2 chain trajectories — diagonal histograms + off-diagonal
scatters, last fraction of iterations, colored by chain.

Used in he2010-synthetic.qmd (R1 and R2). Designed to reveal:
  - Identifiability ridges (elongated correlations)
  - Multimodality (separated clusters by chain)
  - Bound pinning (mass at parameter limits)
  - Cross-chain agreement / disagreement on each parameter
"""
from __future__ import annotations
import io
from pathlib import Path
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def chain_color_palette(chain_ids):
    """Return a stable {chain_id: rgba} dict.

    Use the same call across every plot in a chapter that displays
    per-chain data, so a chain's color is the same on the ll trace,
    the param trajectory facets, the loglik-eval bar, and the pair plot.

    Args:
      chain_ids: iterable of chain IDs (any sortable type, typically int).

    Returns:
      dict mapping each id to an rgba tuple.
    """
    cids = sorted(set(chain_ids))
    if len(cids) <= 10:
        cmap = plt.get_cmap("tab10")
        return {c: cmap(i % 10) for i, c in enumerate(cids)}
    if len(cids) <= 20:
        cmap = plt.get_cmap("tab20")
        return {c: cmap(i % 20) for i, c in enumerate(cids)}
    cmap = plt.get_cmap("viridis")
    return {c: cmap(i / max(1, len(cids)-1)) for i, c in enumerate(cids)}


def load_chain_traces(scout_dir: Path) -> pl.DataFrame:
    """Long-format dataframe of all chains' parameter traces."""
    dfs = []
    chain_dirs = [p for p in scout_dir.glob("chain_*")
                  if p.is_dir() and p.name.split("_")[1].isdigit()]
    for cd in sorted(chain_dirs, key=lambda p: int(p.name.split("_")[1])):
        cid = int(cd.name.split("_")[1])
        lines = [l for l in (cd/"parameter_traces.tsv").read_text().splitlines()
                 if not l.startswith("#")]
        df = pl.read_csv(io.StringIO("\n".join(lines)), separator="\t",
                         infer_schema_length=10000)
        dfs.append(df.with_columns(pl.lit(cid).alias("chain")))
    return pl.concat(dfs)


def pair_plot_chains(
    traces: pl.DataFrame,
    params: list[str],
    truth: dict[str, float] | None = None,
    bounds: dict[str, tuple[float, float]] | None = None,
    last_frac: float = 0.5,
    title: str = "",
    fig_size: tuple[float, float] | None = None,
    point_size: float = 6,
    chain_alpha: float = 0.6,
):
    """Pair plot: diagonal = per-chain marginal histograms, off-diagonal =
    chain-colored scatter of the last `last_frac` of iterations.

    Args:
      traces: long-format DataFrame with columns [iteration, chain, *params].
      params: list of parameter names to plot.
      truth: optional {param: value} dict — drawn as black × on scatters,
             vertical line on diagonals.
      bounds: optional {param: (lo, hi)} — drawn as red dashed lines on
              diagonals to make bound-pinning visible.
      last_frac: fraction of iterations to plot (default 0.5 = last half).
      title: figure-level title.
      fig_size: (w, h); auto-sized from #params if None.
      point_size: scatter marker size.
      chain_alpha: scatter alpha.

    Returns:
      (fig, axes) tuple. axes is the (n × n) array of axes.
    """
    n = len(params)
    max_iter = int(traces["iteration"].max())
    iter_cut = int(max_iter * (1 - last_frac))
    tail = traces.filter(pl.col("iteration") >= iter_cut)
    chain_ids = sorted(tail["chain"].unique().to_list())
    colors = chain_color_palette(chain_ids)

    if fig_size is None:
        fig_size = (1.9 * n + 1, 1.9 * n + 1)
    fig, axes = plt.subplots(n, n, figsize=fig_size)
    if n == 1:
        axes = np.array([[axes]])

    for i, p_y in enumerate(params):
        for j, p_x in enumerate(params):
            ax = axes[i, j]
            if i == j:
                # Diagonal: one histogram per chain, overlaid
                vals_all = tail[p_y].to_numpy()
                vals_all = vals_all[np.isfinite(vals_all)]
                if vals_all.size == 0:
                    ax.set_visible(False); continue
                lo, hi = float(vals_all.min()), float(vals_all.max())
                # Pad if degenerate
                if hi - lo < 1e-9:
                    pad = max(abs(lo) * 0.01, 1e-6)
                    lo, hi = lo - pad, hi + pad
                bins = np.linspace(lo, hi, 25)
                for c in chain_ids:
                    sub = tail.filter(pl.col("chain") == c)[p_y].to_numpy()
                    sub = sub[np.isfinite(sub)]
                    if sub.size:
                        ax.hist(sub, bins=bins, color=colors[c],
                                alpha=0.5, edgecolor="none")
                if truth and p_y in truth:
                    ax.axvline(truth[p_y], color="black", linestyle="--",
                               linewidth=1.5, zorder=10)
                if bounds and p_y in bounds:
                    for b in bounds[p_y]:
                        ax.axvline(b, color="#cc3333", linestyle=":",
                                   linewidth=1.0, alpha=0.7)
                ax.set_yticks([])
            elif i > j:
                # Lower triangle: scatter colored by chain
                for c in chain_ids:
                    sub = tail.filter(pl.col("chain") == c)
                    if sub.height:
                        ax.scatter(sub[p_x], sub[p_y],
                                   color=colors[c], s=point_size,
                                   alpha=chain_alpha,
                                   edgecolor="none")
                if truth and p_x in truth and p_y in truth:
                    ax.scatter([truth[p_x]], [truth[p_y]], marker="x",
                               color="black", s=80, linewidth=2.5,
                               zorder=10)
                if bounds:
                    if p_x in bounds:
                        for b in bounds[p_x]:
                            ax.axvline(b, color="#cc3333", linestyle=":",
                                       linewidth=0.7, alpha=0.4)
                    if p_y in bounds:
                        for b in bounds[p_y]:
                            ax.axhline(b, color="#cc3333", linestyle=":",
                                       linewidth=0.7, alpha=0.4)
            else:
                ax.set_visible(False)

            # Labels only on outer edges
            if i == n - 1:
                ax.set_xlabel(p_x, fontsize=9)
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel(p_y, fontsize=9)
            else:
                ax.set_yticklabels([])
            ax.tick_params(labelsize=7)

    # Legend: one entry per chain, on the upper-right empty space
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=colors[c], markersize=6,
                          label=f"chain {c}") for c in chain_ids]
    if n >= 2:
        legend_ax = axes[0, n-1]
        legend_ax.set_visible(True)
        legend_ax.axis("off")
        legend_ax.legend(handles=handles, frameon=False, fontsize=7,
                         loc="center", ncol=2 if len(chain_ids) > 8 else 1)

    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig, axes
