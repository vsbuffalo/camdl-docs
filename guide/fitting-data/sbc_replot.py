"""Regenerate SBC plot + summary from existing raw TSV."""
from pathlib import Path
import polars as pl, numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG, OUT = ROOT/"figures", ROOT/"output"
GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

# truth = real-data MLE
mle_path = next((ROOT/"results/fits").glob("fit-*/real/fit_42/refine/mle_params.toml"))
truth = {}
for ln in mle_path.read_text().splitlines():
    s = ln.split("#",1)[0].strip()
    if s.startswith("[") or "=" not in s: continue
    k,v = s.split("=",1)
    try: truth[k.strip()] = float(v.strip())
    except: pass

df = pl.read_csv(OUT/"sbc_raw.tsv", separator="\t")
ok   = df.filter(pl.col("status") == "ok")
conv = ok.filter(pl.col("converged"))
conv_rate = conv.height / df.height

CONV_WARN_THRESHOLD = 0.80

def summarise(subset, tag):
    rows = []
    for p, true in [("beta", truth["beta"]), ("gamma", truth["gamma"])]:
        vals = subset[p].to_numpy()
        if len(vals) == 0: continue
        rows.append({"subset": tag, "param": p, "n": len(vals), "truth": true,
                     "mean_est": float(vals.mean()),
                     "bias": float(vals.mean() - true),
                     "sd":   float(vals.std()),
                     "rmse": float(np.sqrt(np.mean((vals - true)**2))),
                     "q025": float(np.quantile(vals, 0.025)),
                     "q975": float(np.quantile(vals, 0.975))})
    return rows

srows = summarise(conv, "converged_only") + summarise(ok, "all_ok")
summary = pl.DataFrame(srows)
summary.write_csv(OUT/"sbc_summary.tsv", separator="\t")
print(summary)

fig, axes = plt.subplots(2, 2, figsize=(11, 7))

for ax, p, c in [(axes[0,0], "beta", "#3366cc"), (axes[0,1], "gamma", "#9b59b6")]:
    vals_all  = ok[p].to_numpy()
    vals_conv = conv[p].to_numpy()
    bins = 6
    ax.hist(vals_all,  bins=bins, color=c, alpha=0.4, edgecolor="none",
            label=f"all ok (n={len(vals_all)})")
    ax.hist(vals_conv, bins=bins, color=c, alpha=0.9, edgecolor="none",
            label=f"converged (n={len(vals_conv)})")
    ax.axvline(truth[p], color="black", linestyle="--",
               label=f"truth = {truth[p]:.3f}")
    ax.axvline(vals_conv.mean(), color="red", linewidth=1.2,
               label=f"mean = {vals_conv.mean():.3f}")
    ax.set_xlabel("β̂" if p=="beta" else "γ̂")
    ax.set_title(f"{p}  bias={vals_conv.mean()-truth[p]:+.3f}  "
                 f"sd={vals_conv.std():.3f}  "
                 f"95% CI [{np.quantile(vals_conv,0.025):.3f}, "
                 f"{np.quantile(vals_conv,0.975):.3f}]", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top","right"]].set_visible(False)

ax = axes[1,0]
# Shared bin edges so β and γ Rhat histograms are directly comparable
rh_all = np.concatenate([ok["rhat_beta"].to_numpy(),
                         ok["rhat_gamma"].to_numpy()])
rh_bins = np.linspace(rh_all.min(), rh_all.max(), 16)
ax.hist(ok["rhat_beta"].to_numpy(),  bins=rh_bins, color="#3366cc", alpha=0.5,
        edgecolor="none", label="β Rhat")
ax.hist(ok["rhat_gamma"].to_numpy(), bins=rh_bins, color="#9b59b6", alpha=0.5,
        edgecolor="none", label="γ Rhat")
ax.axvline(1.10, color="#c0392b", linestyle="--", label="gate = 1.10")
ax.axvline(1.05, color="#e67e22", linestyle=":",  label="target = 1.05")
ax.set_xlabel("refine Rhat")
ax.set_title(f"Refine Rhat per synthetic fit "
             f"({conv.height}/{ok.height} pass gate)")
ax.legend(frameon=False, fontsize=8)
ax.spines[["top","right"]].set_visible(False)

ax = axes[1,1]; ax.axis("off")
lines = [
    "Convergence summary",
    f"    {df.height} synthetic datasets",
    f"    {ok.height} produced an MLE ({100*ok.height/df.height:.0f}%)",
    f"    {conv.height} converged at refine ({100*conv_rate:.0f}%)",
    "",
    "    β Rhat   median {:.3f}   max {:.3f}".format(
        float(ok['rhat_beta'].median()), float(ok['rhat_beta'].max())),
    "    γ Rhat   median {:.3f}   max {:.3f}".format(
        float(ok['rhat_gamma'].median()), float(ok['rhat_gamma'].max())),
    "",
    "    threshold: all refine Rhat < 1.10 must hold",
    f"    warn if converged rate < {100*CONV_WARN_THRESHOLD:.0f}%",
    "",
    "    " + ("✓ PASS" if conv_rate >= CONV_WARN_THRESHOLD else "⚠ WARN"),
]
for i, ln in enumerate(lines):
    ax.text(0.02, 0.95 - i*0.08, ln, transform=ax.transAxes,
            fontfamily="monospace", fontsize=10)

fig.suptitle(f"Synthetic calibration at real-data MLE "
             f"({df.height} synthetic datasets from θ* = "
             f"β={truth['beta']:.3f}, γ={truth['gamma']:.3f})")
fig.tight_layout()
fig.savefig(FIG/"sbc.png", bbox_inches="tight"); plt.close(fig)
print(f"wrote {FIG/'sbc.png'}")
