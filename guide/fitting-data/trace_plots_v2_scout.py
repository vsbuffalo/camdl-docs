"""Scout trace plots for the v2 free-IC fit (narrower bounds, longer scout).

128 scout chains × 1000 iterations — too many to distinguish individually,
so all chains drawn in light gray, best chain (114) highlighted.
"""
from pathlib import Path
import polars as pl, numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

fit_dir = next((ROOT/"results/fits").glob("fit_free_init_v2-*/real/fit_42"))
scout_dir = fit_dir/"scout"

dfs = []
for c in sorted(scout_dir.glob("chain_*")):
    if not c.is_dir(): continue
    d = pl.read_csv(c/"parameter_traces.tsv", separator="\t",
                    comment_prefix="#").sort("iteration")
    cid = int(c.name.split("_")[1])
    dfs.append(d.with_columns(pl.lit(cid).alias("chain")))
scout = pl.concat(dfs)
chains = sorted(int(c) for c in scout["chain"].unique().to_list())
print(f"loaded {len(chains)} scout chains, {scout.height} total rows")

# Thin to every 5 iter so 128×200 = 25600 line segments per panel, readable
scout_thin = scout.filter(pl.col("iteration") % 5 == 0)

# best chain by last-iteration if2_perturbed_loglik
final = scout.filter(pl.col("iteration") == scout["iteration"].max()).sort("chain")
best_chain = int(final.sort("if2_perturbed_loglik", descending=True)["chain"][0])
print(f"best chain (by final if2_perturbed_loglik): {best_chain}")

params = ["beta", "gamma", "I0", "R_init", "if2_perturbed_loglik"]
labels = ["β", "γ", "I(0)", "R(0)", "log-lik (perturbed)"]
fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
ax_list = axes.flat

for ax, p, lab in zip(ax_list, params, labels):
    # Non-best chains: thin gray lines
    for cid, grp in scout_thin.group_by("chain", maintain_order=True):
        c = int(cid[0])
        if c == best_chain: continue
        g = grp.sort("iteration")
        ax.plot(g["iteration"], g[p], color="#999999", alpha=0.18,
                linewidth=0.45, zorder=3)
    # Best chain: highlight in red
    gb = scout_thin.filter(pl.col("chain") == best_chain).sort("iteration")
    ax.plot(gb["iteration"], gb[p], color="#e74c3c", alpha=0.95,
            linewidth=1.2, zorder=5, label=f"chain {best_chain} (best)")
    ax.set_ylabel(lab)
    ax.spines[["top","right"]].set_visible(False)

ax_list[0].legend(frameon=False, fontsize=9, loc="upper right")

# summary panel
ax_last = ax_list[-1]
ax_last.axis("off")
stats = [
    f"scout endpoints across {len(chains)} chains:",
    f"  β      {float(final['beta'].min()):.3f} – {float(final['beta'].max()):.3f}",
    f"  γ      {float(final['gamma'].min()):.3f} – {float(final['gamma'].max()):.3f}",
    f"  I(0)   {float(final['I0'].min()):.1f} – {float(final['I0'].max()):.1f}",
    f"  R(0)   {float(final['R_init'].min()):.0f} – {float(final['R_init'].max()):.0f}",
    f"  ll     {float(final['if2_perturbed_loglik'].min()):.1f} – "
        f"{float(final['if2_perturbed_loglik'].max()):.1f}",
    f"  spread {float(final['if2_perturbed_loglik'].max()-final['if2_perturbed_loglik'].min()):.1f} nats",
    "",
    "scout tail-Rhat:",
    "  β  = 2.237   ✗ (gated)",
    "  γ  = 1.928   ✗ (gated)",
    "  I0 = 8.831   — ungated (ivp)",
    "  R0 = 3.321   — ungated (ivp)",
    "",
    "27/128 chains diverged",
    "refine blocked by scout-Rhat gate",
]
for i, ln in enumerate(stats):
    ax_last.text(0.02, 0.97 - i*0.060, ln, transform=ax_last.transAxes,
                 fontfamily="monospace", fontsize=9)

axes[2,0].set_xlabel("iteration")
fig.suptitle("v2 scout traces (narrower bounds, 128 chains × 1000 iter) — multi-modal",
             fontsize=11, y=0.995)
fig.tight_layout()
out = FIG/"traces_v2_scout.png"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"wrote {out}")
