"""Transmission-tree utilities: Newick parsing, plotting, and shape statistics.

camdl's ``lineage tree`` emits a **forest** of time-calibrated transmission
trees in Newick (one ``;``-terminated tree per surviving introduction). This
module parses that output into a small typed tree structure and provides the
plotting + summary statistics a phylodynamics chapter needs, with no external
phylogenetics dependency (stdlib + numpy + matplotlib only).

The Newick we consume is machine-generated and regular — tip labels look like
``ind2158`` and every node carries a ``:branch_length``. The parser is a
focused recursive-descent reader for exactly that grammar, not a general-purpose
Newick library.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


# --------------------------------------------------------------------------
# Types
# --------------------------------------------------------------------------

@dataclass
class TreeNode:
    """A node in a transmission tree.

    ``branch_length`` is the length of the branch *immediately above* this node
    (connecting it to its parent). ``deme`` is optional metadata joined in from
    the line list (see :func:`annotate_demes`).
    """

    label: str | None = None
    branch_length: float = 0.0
    children: list["TreeNode"] = field(default_factory=list)
    deme: int | None = None

    @property
    def is_tip(self) -> bool:
        return not self.children

    def tips(self) -> list["TreeNode"]:
        """All tip (leaf) nodes under this node, left-to-right."""
        if self.is_tip:
            return [self]
        out: list[TreeNode] = []
        for c in self.children:
            out.extend(c.tips())
        return out

    def n_tips(self) -> int:
        return 1 if self.is_tip else sum(c.n_tips() for c in self.children)


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def parse_newick(text: str) -> list[TreeNode]:
    """Parse a Newick string (possibly multiple ``;``-terminated trees) into a
    forest. Whitespace between trees (including newlines) is ignored.
    """
    forest: list[TreeNode] = []
    for chunk in text.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        node, idx = _parse_clade(chunk, 0)
        forest.append(node)
    return forest


def _parse_clade(s: str, i: int) -> tuple[TreeNode, int]:
    node = TreeNode()
    if i < len(s) and s[i] == "(":
        # internal node: parse a comma-separated child list until ")"
        i += 1
        while True:
            child, i = _parse_clade(s, i)
            node.children.append(child)
            if i < len(s) and s[i] == ",":
                i += 1
                continue
            if i < len(s) and s[i] == ")":
                i += 1
                break
            break
    # label (tips, and optionally internal nodes) up to : , ) end
    start = i
    while i < len(s) and s[i] not in ":,)":
        i += 1
    if i > start:
        node.label = s[start:i]
    # branch length
    if i < len(s) and s[i] == ":":
        i += 1
        start = i
        while i < len(s) and s[i] not in ",)":
            i += 1
        node.branch_length = float(s[start:i])
    return node, i


def tip_ids(forest: list[TreeNode]) -> list[int]:
    """Integer individual ids parsed from ``indNNNN`` tip labels across a forest."""
    out: list[int] = []
    for root in forest:
        for t in root.tips():
            if t.label and t.label.startswith("ind"):
                out.append(int(t.label[3:]))
    return out


def annotate_demes(forest: list[TreeNode], id_to_deme: dict[int, int]) -> None:
    """Set ``.deme`` on every tip by looking its integer id up in ``id_to_deme``."""
    for root in forest:
        for t in root.tips():
            if t.label and t.label.startswith("ind"):
                t.deme = id_to_deme.get(int(t.label[3:]))


# --------------------------------------------------------------------------
# Layout + plotting
# --------------------------------------------------------------------------

def _layout(root: TreeNode) -> tuple[dict[int, float], dict[int, float], float]:
    """Assign (x, y) to every node. x = cumulative branch length from the root;
    y = tip order for tips, mean of children for internals. Returns
    ``(x_by_id, y_by_id, n_tips)`` keyed by ``id(node)``.
    """
    x: dict[int, float] = {}
    y: dict[int, float] = {}
    counter = [0]

    def assign_x(node: TreeNode, acc: float) -> None:
        xv = acc + node.branch_length
        x[id(node)] = xv
        for c in node.children:
            assign_x(c, xv)

    def assign_y(node: TreeNode) -> float:
        if node.is_tip:
            yv = float(counter[0])
            counter[0] += 1
        else:
            ys = [assign_y(c) for c in node.children]
            yv = sum(ys) / len(ys)
        y[id(node)] = yv
        return yv

    assign_x(root, 0.0)
    assign_y(root)
    return x, y, float(counter[0])


def plot_forest(
    forest: list[TreeNode],
    ax,
    *,
    color_by_deme: bool = False,
    palette: dict[int, object] | None = None,
    default_color: str = "#555555",
    linewidth: float = 0.6,
    tip_size: float = 6.0,
    gap: float = 2.0,
    max_components: int | None = None,
):
    """Draw a rectangular phylogram for each tree, stacked vertically.

    Branch lengths are time-calibrated, so the x-axis is time. Tip marker color
    encodes ``deme`` when ``color_by_deme`` and a ``palette`` are supplied.
    Returns the total y-extent used (useful for sizing the figure).
    """
    components = forest if max_components is None else forest[:max_components]
    y_offset = 0.0
    for root in components:
        x, y, ntips = _layout(root)

        def draw(node: TreeNode) -> None:
            xn, yn = x[id(node)], y[id(node)] + y_offset
            for c in node.children:
                xc, yc = x[id(c)], y[id(c)] + y_offset
                ax.plot([xn, xn], [yn, yc], color=default_color, lw=linewidth, zorder=1)
                ax.plot([xn, xc], [yc, yc], color=default_color, lw=linewidth, zorder=1)
                draw(c)
            if node.is_tip and color_by_deme and palette is not None:
                col = palette.get(node.deme, default_color)
                ax.scatter([xn], [yn], s=tip_size, color=col, zorder=2, edgecolors="none")

        draw(root)
        y_offset += ntips + gap

    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("time (days)")
    return y_offset


# --------------------------------------------------------------------------
# Tree-shape statistics
# --------------------------------------------------------------------------

def sackin_index(root: TreeNode, normalize: bool = False) -> float:
    """Sackin's index: sum of root-to-tip depths (number of internal ancestors).

    Larger = more imbalanced (caterpillar-like). With ``normalize=True``, divide
    by the tip count to get mean tip depth (comparable across tree sizes).
    """
    total = [0]

    def walk(node: TreeNode, depth: int) -> None:
        if node.is_tip:
            total[0] += depth
        else:
            for c in node.children:
                walk(c, depth + 1)

    walk(root, 0)
    n = root.n_tips()
    return total[0] / n if (normalize and n) else float(total[0])


def colless_index(root: TreeNode, normalize: bool = False) -> float:
    """Colless imbalance: sum over internal nodes of |L - R| tip counts.

    Defined for strictly bifurcating trees; for multifurcations we use the
    spread (max - min) of child subtree sizes as the local imbalance. With
    ``normalize=True``, divide by the max possible ((n-1)(n-2)/2) for n tips.
    """
    total = [0.0]

    def walk(node: TreeNode) -> int:
        if node.is_tip:
            return 1
        sizes = [walk(c) for c in node.children]
        total[0] += max(sizes) - min(sizes)
        return sum(sizes)

    walk(root)
    if normalize:
        n = root.n_tips()
        denom = (n - 1) * (n - 2) / 2 if n > 2 else 1
        return total[0] / denom
    return total[0]


def ltt(forest: list[TreeNode]) -> tuple[np.ndarray, np.ndarray]:
    """Lineages-through-time across a forest.

    Each branching event in absolute time adds (children - 1) lineages; each tip
    removes one. Returns ``(times, n_lineages)`` as a step function, summed over
    all forest components on a shared absolute-time axis.
    """
    events: list[tuple[float, int]] = []
    for root in forest:
        x, _, _ = _layout(root)
        # root's branch starts the clock for its own introduction
        events.append((x[id(root)] - root.branch_length, 1))

        def walk(node: TreeNode) -> None:
            xn = x[id(node)]
            if node.is_tip:
                events.append((xn, -1))
            else:
                events.append((xn, len(node.children) - 1))
                for c in node.children:
                    walk(c)

        walk(root)

    events.sort()
    times = np.array([t for t, _ in events])
    counts = np.cumsum([d for _, d in events])
    return times, counts


def skyline_from_tree(root: TreeNode, window: int = 15) -> tuple[np.ndarray, np.ndarray]:
    """Generalized-skyline estimate of relative prevalence through time from a
    single time-calibrated tree.

    In each interval between consecutive coalescent (internal-node) events, with
    ``A`` lineages present, the classic skyline estimator sets the coalescent
    ``N_e`` to ``Δt · C(A, 2)``. For an epidemic genealogy the pairwise
    coalescent rate is ``∝ 1/I(t)`` (the SIR rate ``2βS/(NI)``), so this ``N_e``
    is proportional to prevalence ``I(t)`` up to a sampling/rate-dependent
    constant — the *shape* is recovered, not the absolute level. Returns
    ``(times, Ne_relative)``, median-smoothed over ``window`` intervals.
    """
    x: dict[int, float] = {}

    def assign_x(node: TreeNode, acc: float) -> None:
        xv = acc + node.branch_length
        x[id(node)] = xv
        for c in node.children:
            assign_x(c, xv)

    assign_x(root, 0.0)
    nodes: list[TreeNode] = []

    def collect(n: TreeNode) -> None:
        nodes.append(n)
        for c in n.children:
            collect(c)

    collect(root)
    coal = sorted(x[id(n)] for n in nodes if n.children)
    t_ltt, n_ltt = ltt([root])

    st, sne = [], []
    for t0, t1 in zip(coal, coal[1:]):
        dt = t1 - t0
        if dt <= 0:
            continue
        A = n_ltt[np.searchsorted(t_ltt, 0.5 * (t0 + t1), side="right") - 1]
        if A >= 2:
            st.append(0.5 * (t0 + t1))
            sne.append(dt * A * (A - 1) / 2)
    st, sne = np.array(st), np.array(sne)
    if window and len(sne) >= window:
        from numpy.lib.stride_tricks import sliding_window_view as swv
        sm = np.median(swv(sne, window), axis=1)
        st = st[window // 2: window // 2 + len(sm)]
        sne = sm
    return st, sne


def plot_tree_with_skyline(
    root: TreeNode,
    prevalence_t=None,
    prevalence=None,
    *,
    figsize: tuple[float, float] = (8.0, 7.0),
    window: int = 15,
    skyline_color: str = "#d62728",
    prevalence_color: str = "#333333",
    **tree_kwargs,
):
    """Stacked two-panel figure sharing a time axis: **top** the tree-derived
    skyline (rescaled to overlay) and, if supplied, the true prevalence
    ``I(t)``; **bottom** the tree itself (drawn by :func:`plot_forest`). Returns
    ``(fig, (ax_top, ax_bottom))``.
    """
    import matplotlib.pyplot as plt

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=figsize, sharex=True, gridspec_kw={"height_ratios": [1, 2.2]}
    )
    st, sne = skyline_from_tree(root, window=window)
    scale = 1.0
    if prevalence is not None and prevalence_t is not None and len(sne):
        # least-squares scale (through origin) over the tree's coalescent window,
        # so the skyline sits on the prevalence curve where it has signal
        I_at = np.interp(st, np.asarray(prevalence_t), np.asarray(prevalence))
        denom = float(np.sum(sne * sne))
        scale = float(np.sum(I_at * sne) / denom) if denom > 0 else 1.0
    if prevalence_t is not None and prevalence is not None:
        ax_top.plot(prevalence_t, prevalence, color=prevalence_color, lw=2,
                    label="true prevalence $I(t)$")
    ax_top.plot(st, sne * scale, color=skyline_color, lw=1.4,
                label="skyline (from tree, rescaled)")
    ax_top.set_ylabel("prevalence")
    ax_top.legend(loc="upper right", frameon=False, fontsize=8)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    plot_forest([root], ax_bot, **tree_kwargs)
    fig.tight_layout()
    return fig, (ax_top, ax_bot)


def offspring_counts(parent_ids, infected_ids) -> np.ndarray:
    """Per-individual secondary-case counts over the **full** set of infected
    individuals — the offspring distribution in the Lloyd-Smith et al. (2005)
    sense, which *must* include individuals who infected nobody (the zeros).

    ``parent_ids``: the infector id for each transmission event (a line list's
    ``parent_id`` column, filtered to ``parent_kind == "individual"``).
    ``infected_ids``: the universe of individuals who were ever infected (each
    contributes exactly one offspring count, 0 if it never transmitted).

    Omitting the zeros (counting only realized infectors) inflates the mean and
    deflates the variance/mean ratio, which biases :func:`dispersion_k` upward
    and can hide genuine superspreading.
    """
    from collections import Counter

    c = Counter(int(p) for p in parent_ids if p is not None and p >= 0)
    return np.array([c.get(int(i), 0) for i in infected_ids])


def dispersion_k(offspring_counts: np.ndarray) -> float:
    """Negative-binomial dispersion ``k`` (method-of-moments) for an offspring
    distribution. Small k (<1) = overdispersion / superspreading; k -> inf is
    Poisson. Returns ``inf`` when the variance does not exceed the mean.
    """
    m = float(np.mean(offspring_counts))
    v = float(np.var(offspring_counts))
    if v <= m or m == 0:
        return float("inf")
    return m * m / (v - m)
