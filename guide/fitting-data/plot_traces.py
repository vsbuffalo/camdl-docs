"""Reusable trace-plot function for IF2 chain diagnostics.

Usage:
    from plot_traces import load_chains, plot_traces

    chains = load_chains(fit_dir/"scout")
    plot_traces(chains,
                params=["beta","gamma","I0","R_init","if2_perturbed_loglik"],
                labels=["β","γ","I(0)","R(0)","log-lik"],
                log_params={"R_init"},
                outpath="figures/traces.png",
                title="Poisson free-IC scout traces")
"""
from pathlib import Path
import polars as pl, numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

GREY = "#555555"


def _apply_rc():
    mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
        "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,
        "text.color":GREY})


def load_chains(stage_dir):
    """Load all chain parameter traces under `stage_dir/chain_*/parameter_traces.tsv`
    into one polars DataFrame with an added 'chain' column. Parses "NA" as null."""
    dfs = []
    for c in sorted(Path(stage_dir).glob("chain_*")):
        if not c.is_dir(): continue
        d = pl.read_csv(c/"parameter_traces.tsv", separator="\t",
                        comment_prefix="#", null_values="NA").sort("iteration")
        cid = int(c.name.split("_")[1])
        dfs.append(d.with_columns(pl.lit(cid).alias("chain")))
    return pl.concat(dfs) if dfs else pl.DataFrame()


def best_chain_by_final_ll(chains, ll_col="loglik"):
    """Return the chain ID with the highest *final* non-null ll value
    (matches camdl's 'best chain' convention). Falls back to perturbed
    ll if clean ll column is all null."""
    last = (chains.filter(pl.col(ll_col).is_not_null())
            .sort("iteration")
            .group_by("chain", maintain_order=True)
            .agg(pl.col(ll_col).last().alias("final_ll")))
    if last.height == 0 and ll_col == "loglik":
        return best_chain_by_final_ll(chains, ll_col="if2_perturbed_loglik")
    return int(last.sort("final_ll", descending=True)["chain"][0])


def plot_traces(chains, params, labels, *, outpath,
                title="", log_params=None,
                best_chain=None, highlight_color="#e74c3c",
                color_mode="tab20", bounds=None,
                n_cols=2, figsize=None):
    """Plot chain traces for each parameter in `params`.

    Args:
        chains: polars DataFrame with columns `iteration`, `chain`, and
            each of `params`.
        params: list of column names to plot.
        labels: list of display labels, same length as `params`.
        outpath: figure output path.
        title: suptitle.
        log_params: iterable of param names to plot with symlog y-axis.
        best_chain: chain ID to highlight; None to disable.
        highlight_color: color for the best chain.
        color_mode: "tab20" (arbitrary colours, one per chain) or "ll"
            (viridis keyed on each chain's final clean log-lik — bright
            = good, dark = bad; requires `loglik` column in `chains`).
        bounds: dict {param_name: (lo, hi)} to draw as dashed horizontal
            reference lines (shows when chains pin at parameter bounds).
        n_cols: subplot grid columns.
        figsize: override figure size.
    """
    _apply_rc()
    log_params = set(log_params or ())
    bounds = bounds or {}
    chain_ids = sorted(int(c) for c in chains["chain"].unique().to_list())

    if color_mode == "ll":
        last = (chains.filter(pl.col("loglik").is_not_null())
                .sort("iteration")
                .group_by("chain", maintain_order=True)
                .agg(pl.col("loglik").last().alias("final_ll")))
        by_chain = {int(r["chain"]): float(r["final_ll"])
                    for r in last.iter_rows(named=True)}
        lls = np.array([by_chain.get(c, np.nan) for c in chain_ids])
        # Normalise so best (max ll) → 1.0, clip at 30 nats below best
        lo = np.nanmax(lls) - 30
        norm = np.clip((lls - lo) / (np.nanmax(lls) - lo + 1e-9), 0, 1)
        cmap = {c: plt.cm.viridis(norm[i]) for i, c in enumerate(chain_ids)}
    else:
        palette = plt.cm.tab20(np.linspace(0, 1, 20))
        cmap = {c: palette[i % 20] for i, c in enumerate(chain_ids)}

    n = len(params)
    n_rows = (n + n_cols - 1) // n_cols
    if figsize is None:
        figsize = (6 * n_cols, 2.6 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=True)
    axes = axes.flatten()

    for ax, p, lab in zip(axes, params, labels):
        for cid, grp in chains.group_by("chain", maintain_order=True):
            c = int(cid[0])
            g = grp.sort("iteration")
            is_best = (best_chain is not None and c == best_chain)
            col = highlight_color if is_best else cmap[c]
            alpha = 0.95 if is_best else (0.65 if color_mode == "ll" else 0.35)
            lw    = 1.3  if is_best else 0.5
            z     = 5    if is_best else 3
            ax.plot(g["iteration"], g[p], color=col, alpha=alpha,
                    linewidth=lw, zorder=z)
        # bound reference lines (dashed grey)
        if p in bounds:
            for y in bounds[p]:
                ax.axhline(y, color=GREY, linestyle=":", linewidth=0.8,
                           alpha=0.7, zorder=2)
        ax.set_ylabel(lab)
        ax.spines[["top","right"]].set_visible(False)
        if p in log_params:
            ax.set_yscale("symlog")

    for ax in axes[n:]:
        ax.axis("off")

    # if using ll colouring, add a colorbar in the empty panel if one exists
    if color_mode == "ll" and n < n_rows * n_cols:
        sm = plt.cm.ScalarMappable(
            cmap="viridis",
            norm=mpl.colors.Normalize(vmin=lo, vmax=np.nanmax(lls)))
        cbar = fig.colorbar(sm, ax=axes[n:], orientation="vertical",
                             fraction=0.05, pad=0.05, shrink=0.6,
                             label="final clean log-lik")

    bottom_row = axes[n_cols * (n_rows - 1):n_cols * n_rows]
    for ax in bottom_row[:n - n_cols * (n_rows - 1)]:
        ax.set_xlabel("iteration")

    if title:
        fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath


def plot_rank_overlay(chains, params, labels, *, outpath,
                      title="", n_bins=20, n_cols=2, figsize=None,
                      best_chain=None, highlight_color="#e74c3c"):
    """Rank-overlay plot à la Vehtari et al. 2021.

    For each parameter, compute each chain's rank within each iteration
    across chains, then histogram the ranks each chain visited. A
    well-mixed chain visits all rank bins uniformly (flat histogram);
    a stuck chain concentrates at low or high ranks (spike).

    Overlays all chains' rank histograms per panel — uniform spaghetti
    of flat lines → mixed; distinct shapes → unmixed.

    Args:
        chains: polars DataFrame with columns `iteration`, `chain`, params.
        params: list of param columns.
        labels: display labels.
        outpath: output path.
        n_bins: number of rank bins (default 20).
        best_chain: optional chain ID to highlight.
    """
    _apply_rc()
    chain_ids = sorted(int(c) for c in chains["chain"].unique().to_list())
    n_chains = len(chain_ids)

    n = len(params)
    n_rows = (n + n_cols - 1) // n_cols
    if figsize is None:
        figsize = (6 * n_cols, 2.8 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()

    # Pre-pivot: (iteration × chain) matrix of values per parameter
    for ax, p, lab in zip(axes, params, labels):
        wide = chains.select(["iteration", "chain", p]).pivot(
            values=p, index="iteration", on="chain").sort("iteration")
        cols = [c for c in wide.columns if c != "iteration"]
        mat = wide.select(cols).to_numpy()   # shape (n_iter, n_chains)
        mask = np.isfinite(mat)
        # For each iteration compute ranks across chains (break ties by mean)
        ranks = np.full_like(mat, np.nan, dtype=float)
        for i, row in enumerate(mat):
            finite = np.isfinite(row)
            if finite.sum() < 2: continue
            r = np.empty_like(row, dtype=float)
            r[finite] = _rankdata(row[finite])
            r[~finite] = np.nan
            ranks[i] = r

        # For each chain, histogram of its ranks across iterations.
        bin_edges = np.linspace(0, n_chains + 1, n_bins + 1)
        for j, cid in enumerate([int(c) for c in cols]):
            r = ranks[:, j]
            r = r[np.isfinite(r)]
            if r.size == 0: continue
            h, _ = np.histogram(r, bins=bin_edges, density=True)
            # plot as stepped line
            x = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            is_best = (best_chain is not None and cid == best_chain)
            col = highlight_color if is_best else "#2c3e50"
            alpha = 0.95 if is_best else 0.15
            lw    = 1.4  if is_best else 0.6
            ax.plot(x, h, color=col, alpha=alpha, linewidth=lw,
                    zorder=5 if is_best else 3)
        # expected-uniform reference line
        exp_density = 1.0 / (n_chains + 1)
        ax.axhline(exp_density, color="#e67e22", linestyle="--", linewidth=1,
                   alpha=0.8, label="uniform (mixed)")
        ax.set_xlabel("within-iteration rank")
        ax.set_ylabel("density")
        ax.set_title(lab, fontsize=10, pad=4)
        ax.spines[["top","right"]].set_visible(False)

    axes[0].legend(frameon=False, fontsize=8, loc="upper right")
    for ax in axes[n:]:
        ax.axis("off")

    if title:
        fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath


def plot_pair_grid(chains, params, labels, *, outpath,
                    title="", ll_col="if2_perturbed_loglik",
                    dll_filter=None, best_chain=None,
                    highlight_color="#e74c3c", n_bins=14,
                    ll_axis_dll=None, figsize=None):
    """Lower-triangle pair plot of IF2 chain trajectories + endpoints,
    with marginal histograms on the diagonal.

    Optionally append log-likelihood as the last "parameter" — a
    standard exploratory-MCMC convention (pair plots / corner plots
    predate Stan's `lp__` column; the idea is the same: visualise how
    likelihood co-varies with each estimated parameter).

    Args:
        chains: polars DataFrame with columns `iteration`, `chain`, the
            params, and (if included) `ll_col`.
        params: list of column names to plot. If `ll_col` is in the list,
            its axis can be clipped via `ll_axis_dll`.
        labels: display labels.
        outpath: output path.
        ll_col: column used to colour endpoint dots and to sort chains
            by quality for filtering.
        dll_filter: if set, only show chains whose final `ll_col` value
            is within `dll_filter` of the best chain's — keeps the
            colormap informative by excluding outliers.
        best_chain: chain ID to highlight (trajectory only).
        n_bins: histogram bin count on the diagonal.
        ll_axis_dll: if `ll_col` is in `params`, clip its axes to
            [best_ll - ll_axis_dll, best_ll + 2] so the very-bad
            chains don't compress the informative range. Passing
            None uses full range.
    """
    _apply_rc()
    chain_ids = sorted(int(c) for c in chains["chain"].unique().to_list())
    final = (chains.filter(pl.col("iteration") == chains["iteration"].max())
             .sort("chain"))

    best_ll = float(final[ll_col].max())
    if dll_filter is not None:
        kept = final.filter(pl.col(ll_col) >= best_ll - dll_filter)["chain"].to_list()
        kept_set = set(kept)
        chains = chains.filter(pl.col("chain").is_in(kept))
        final  = final.filter(pl.col("chain").is_in(kept))
        print(f"plot_pair_grid: keeping {len(kept)}/{len(chain_ids)} chains "
              f"(Δ {ll_col} < {dll_filter})")
    thin = chains.filter(pl.col("iteration") % 20 == 0)

    # Tail-Rhat per param (last half of iterations) — same convention as camdl
    def tail_rhat(param):
        n_iter = int(chains["iteration"].max()) + 1
        tail = chains.filter(pl.col("iteration") >= n_iter // 2)
        g = (tail.filter(pl.col(param).is_not_null())
             .group_by("chain")
             .agg(pl.col(param).mean().alias("m"),
                  pl.col(param).var().alias("v"),
                  pl.col(param).count().alias("n")))
        if g.height < 2: return float("nan")
        m = g["m"].to_numpy(); v = g["v"].to_numpy(); nn = g["n"].to_numpy()
        n_tail = int(nn[0])
        if n_tail < 2: return float("nan")
        W = v.mean(); B = n_tail * m.var(ddof=1)
        var_hat = (n_tail - 1) / n_tail * W + B / n_tail
        return float(np.sqrt(var_hat / W)) if W > 0 else float("nan")

    rhats = {p: tail_rhat(p) for p in params}

    n = len(params)
    if figsize is None:
        figsize = (3 * n + 2, 3 * n + 1)
    fig, axes = plt.subplots(n, n, figsize=figsize)

    # log-lik axis clip range
    ll_lo = best_ll - ll_axis_dll if ll_axis_dll is not None else None
    ll_hi = best_ll + 2 if ll_axis_dll is not None else None

    def axlim(param):
        """Return (lo, hi) for axes involving param — respects ll clip."""
        if param == ll_col and ll_axis_dll is not None:
            return (ll_lo, ll_hi)
        vals = chains[param].to_numpy()
        vals = vals[np.isfinite(vals)]
        if vals.size == 0: return (None, None)
        return (float(vals.min()), float(vals.max()))

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            p_y = params[i]; p_x = params[j]
            if i == j:
                # diagonal — marginal histogram of final-iter values
                vals = final[p_y].to_numpy()
                if p_y == ll_col and ll_axis_dll is not None:
                    vals = vals[vals >= ll_lo]
                ax.hist(vals, bins=n_bins, color="#3366cc",
                        alpha=0.65, edgecolor="none")
                if best_chain is not None:
                    bval = final.filter(pl.col("chain") == best_chain)[p_y]
                    if len(bval) > 0:
                        ax.axvline(float(bval[0]), color=highlight_color,
                                   linestyle="--", linewidth=1.2)
                # Tail-Rhat annotation, coloured by gate threshold
                r = rhats.get(p_y, float("nan"))
                if np.isfinite(r):
                    rcolor = ("#27ae60" if r < 1.05 else
                              "#e67e22" if r < 1.10 else "#c0392b")
                    ax.text(0.98, 0.95, f"R̂ = {r:.2f}",
                            transform=ax.transAxes, ha="right", va="top",
                            fontsize=9, color=rcolor, fontweight="bold")
                # keep only the first (top-left) y axis visible
                if i > 0:
                    ax.set_yticks([])
                else:
                    ax.set_ylabel("count")
            elif i > j:
                # lower triangle — chain trajectories + endpoints
                for cid, grp in thin.group_by("chain", maintain_order=True):
                    c = int(cid[0])
                    if best_chain is not None and c == best_chain: continue
                    g = grp.sort("iteration")
                    ax.plot(g[p_x].to_numpy(), g[p_y].to_numpy(),
                            color="#999999", alpha=0.12, linewidth=0.4,
                            zorder=2)
                if best_chain is not None:
                    gb = thin.filter(pl.col("chain") == best_chain).sort("iteration")
                    ax.plot(gb[p_x].to_numpy(), gb[p_y].to_numpy(),
                            color=highlight_color, alpha=0.95,
                            linewidth=0.9, zorder=4)
                    bf = final.filter(pl.col("chain") == best_chain)
                    if bf.height:
                        ax.scatter(bf[p_x], bf[p_y], marker="*",
                                   s=130, color="red", edgecolor="white",
                                   linewidth=1.0, zorder=6)
                # endpoint dots coloured by final ll
                sc = ax.scatter(final[p_x].to_numpy(), final[p_y].to_numpy(),
                                c=final[ll_col].to_numpy(),
                                cmap="viridis", s=14, alpha=0.85,
                                zorder=3, edgecolor="none",
                                vmin=ll_lo, vmax=best_ll)
            else:
                ax.axis("off")
                continue

            lo_y, hi_y = axlim(p_y)
            lo_x, hi_x = axlim(p_x)
            if i == j:
                if lo_x is not None: ax.set_xlim(lo_x, hi_x)
            else:
                if lo_x is not None: ax.set_xlim(lo_x, hi_x)
                if lo_y is not None: ax.set_ylim(lo_y, hi_y)

            if i == n - 1: ax.set_xlabel(labels[j])
            else: ax.set_xticklabels([])
            if j == 0 and i != 0: ax.set_ylabel(labels[i])
            elif i == j: ax.set_ylabel("")
            else: ax.set_yticklabels([])
            ax.spines[["top","right"]].set_visible(False)

    # Colorbar + summary in one of the upper-triangle empty cells
    if n >= 2:
        cax = fig.add_axes([0.75, 0.82, 0.012, 0.12])
        fig.colorbar(sc, cax=cax, label=f"chain final {ll_col}")

    if title:
        fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath


def _rankdata(a):
    """Simple rank (1..n) with ties broken by mean rank."""
    a = np.asarray(a)
    order = np.argsort(a)
    ranks = np.empty_like(a, dtype=float)
    ranks[order] = np.arange(1, len(a) + 1)
    # handle ties by averaging
    unique_vals, inv = np.unique(a, return_inverse=True)
    for v_idx, v in enumerate(unique_vals):
        mask = inv == v_idx
        if mask.sum() > 1:
            ranks[mask] = ranks[mask].mean()
    return ranks
