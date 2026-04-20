"""Classic IF2 diagnostic trace plots — each estimated parameter and
log-lik vs iteration, all 12 refine chains overlaid, chain 12 highlighted.

Tells us whether chains have settled or are still drifting."""
from pathlib import Path
import polars as pl, numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

fit_dir = next((ROOT/"results/fits").glob("fit_negbin-*/real/fit_42"))
def load_chains(stage):
    dfs = []
    for cdir in sorted((fit_dir/stage).glob("chain_*")):
        if not cdir.is_dir():
            continue
        d = pl.read_csv(cdir/"parameter_traces.tsv", separator="\t",
                        comment_prefix="#").sort("iteration")
        cid = int(cdir.name.split("_")[1])
        dfs.append(d.with_columns(pl.lit(cid).alias("chain")))
    return pl.concat(dfs) if dfs else pl.DataFrame()

refine = load_chains("refine")
# best chain by final if2_perturbed_loglik
final_iter = refine["iteration"].max()
best_chain = int(refine.filter(pl.col("iteration") == final_iter)
                 .sort("if2_perturbed_loglik", descending=True)["chain"][0])

params = ["beta", "gamma", "k", "if2_perturbed_loglik"]
labels = ["β", "γ", "k", "log-lik (if2-perturbed)"]
fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True)

chain_ids = sorted(int(c) for c in refine["chain"].unique().to_list())
# tab20 gives ~20 perceptually distinct colors; we use the first 12
palette   = plt.cm.tab20(np.linspace(0, 1, 20))[:len(chain_ids)]
cmap      = {c: palette[i] for i, c in enumerate(chain_ids)}

for ax, p, lab in zip(axes.flat, params, labels):
    for cid, grp in refine.group_by("chain", maintain_order=True):
        c = int(cid[0])
        g = grp.sort("iteration")
        ax.plot(g["iteration"], g[p], color=cmap[c], alpha=0.85,
                linewidth=0.9, zorder=3, label=f"chain {c}" if p=="beta" else None)
    ax.set_ylabel(lab)
    ax.spines[["top","right"]].set_visible(False)
    if p == "k":
        ax.set_yscale("log")

# legend in the first panel only
axes[0,0].legend(frameon=False, fontsize=7, ncol=2, loc="lower right")

axes[1,0].set_xlabel("iteration")
axes[1,1].set_xlabel("iteration")
fig.suptitle(f"Refine chain traces (NegBin SIR) — all {len(chain_ids)} chains, one colour each",
             y=0.995, fontsize=11)
fig.tight_layout()
fig.savefig(FIG/"refine_traces.png", bbox_inches="tight")
plt.close(fig)
print(f"wrote {FIG/'refine_traces.png'}")

# print a per-chain summary of final k
final = refine.filter(pl.col("iteration") == final_iter).sort("chain")
print("\nrefine chain endpoints:")
print(final.select(["chain", "beta", "gamma", "k", "if2_perturbed_loglik"]))
