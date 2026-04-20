"""Simulation-based calibration at the real-data MLE.

1. Read the refine MLE as the synthetic truth θ*.
2. Generate N synthetic datasets from chain_binomial at θ*.
3. Fit each synthetic dataset with the same fit.toml pipeline.
4. Report bias, SD, and RMSE on (β, γ).
"""
import subprocess, shlex, tempfile, json
from pathlib import Path
import polars as pl, numpy as np, matplotlib.pyplot as plt, matplotlib as mpl
import time as time_mod

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"; OUT.mkdir(exist_ok=True)
SYNTH = OUT/"synth_data"; SYNTH.mkdir(exist_ok=True)
SYNTH_FITS = OUT/"synth_fits"; SYNTH_FITS.mkdir(exist_ok=True)

GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

N_SYNTH = 30
BASE_SEED = 5001

def read_toml(p):
    d = {}
    for line in Path(p).read_text().splitlines():
        line = line.split("#",1)[0].strip()
        if "=" in line:
            k,v = line.split("=",1)
            try: d[k.strip()] = float(v.strip())
            except: pass
    return d

mle_candidates = sorted((ROOT/"results/fits").glob("fit-*/real/fit_42/refine/mle_params.toml"))
if not mle_candidates:
    raise SystemExit("No fit found. Run: camdl fit run fit.toml --seed 42")
MLE = mle_candidates[-1]
truth = read_toml(MLE)
print(f"Truth θ*: β={truth['beta']:.4f}, γ={truth['gamma']:.4f}, I(0)={int(truth['I0'])}, N={int(truth['N0'])}")

# ── 1. Generate synthetic datasets (one simulate call, N reps, split by replicate) ──
t0 = time_mod.time()
batch_tsv = SYNTH/"batch.tsv"
subprocess.run(shlex.split(
    f"camdl simulate boarding_school_sir.ir.json --params {MLE} "
    f"--replicates {N_SYNTH} --seed {BASE_SEED} --scenario baseline "
    f"--backend chain_binomial --dt 1.0 --obs {batch_tsv}"),
    check=True, capture_output=True, cwd=str(ROOT))

batch = pl.read_csv(batch_tsv, separator="\t")
for rep in sorted(set(batch["replicate"].to_list())):
    sub = batch.filter(pl.col("replicate") == rep).drop("replicate").sort("time")
    sub.write_csv(SYNTH/f"synth_{rep:03d}.tsv", separator="\t")
print(f"generated {N_SYNTH} synthetic datasets in {time_mod.time()-t0:.1f}s")

# ── 2. Template + run fit for each ──
fit_template = (ROOT/"fit.toml").read_text()

rows = []
for rep in range(1, N_SYNTH+1):
    t0 = time_mod.time()
    synth_tsv = SYNTH/f"synth_{rep:03d}.tsv"
    fit_name = f"synth_{rep:03d}"
    fit_toml = SYNTH_FITS/f"{fit_name}.toml"
    toml_text = fit_template.replace('in_bed = "data/in_bed.tsv"',
                                     f'in_bed = "{synth_tsv.resolve()}"')
    toml_text = toml_text.replace('output_dir = "results/sir"',
                                  f'output_dir = "sbc/{fit_name}"')
    fit_toml.write_text(toml_text)

    log_path = SYNTH_FITS/f"{fit_name}.log"
    with log_path.open("w") as f:
        res = subprocess.run(
            ["camdl", "fit", "run", str(fit_toml), "--seed", "42"],
            stdout=f, stderr=subprocess.STDOUT, cwd=str(ROOT))
    elapsed = time_mod.time() - t0

    refine = ROOT/f"results/fits/{fit_name}-*/real/fit_42/refine"
    refine_glob = list((ROOT/"results/fits").glob(f"{fit_name}-*/real/fit_42/refine"))
    if not refine_glob:
        rows.append({"rep": rep, "status": "no refine dir", "elapsed_s": elapsed})
        continue
    refine_dir = refine_glob[0]
    try:
        hat = read_toml(refine_dir/"mle_params.toml")
        # parse Rhat from log
        import re
        rhats = {}
        for ln in log_path.read_text().splitlines():
            m = re.match(r"\s*(\w+)\s+Rhat=(\d+\.\d+)", ln)
            if m: rhats[m.group(1)] = float(m.group(2))
        # only the final (refine) Rhats — they overwrite scout's
        converged = all(r < 1.10 for r in rhats.values()) if rhats else False
        rows.append({"rep": rep, "status": "ok", "elapsed_s": elapsed,
                     "beta": hat["beta"], "gamma": hat["gamma"],
                     "rhat_beta": rhats.get("beta", float("nan")),
                     "rhat_gamma": rhats.get("gamma", float("nan")),
                     "converged": converged})
        print(f"  rep {rep:2d}: β̂={hat['beta']:.3f} γ̂={hat['gamma']:.3f} "
              f"Rhat={rhats} conv={converged} ({elapsed:.1f}s)")
    except Exception as e:
        rows.append({"rep": rep, "status": f"parse err: {e}", "elapsed_s": elapsed})

df = pl.DataFrame(rows)
df.write_csv(OUT/"sbc_raw.tsv", separator="\t")

# ── 3. Convergence audit ──
ok   = df.filter(pl.col("status") == "ok")
conv = ok.filter(pl.col("converged"))
conv_rate = conv.height / df.height if df.height else 0.0
print(f"\nconvergence audit:")
print(f"  {df.height} total fits")
print(f"  {ok.height} produced an MLE ({100*ok.height/df.height:.0f}%)")
print(f"  {conv.height} converged with all refine Rhat < 1.10 ({100*conv_rate:.0f}%)")
for p in ("rhat_beta", "rhat_gamma"):
    s = ok[p]
    print(f"  {p}: median {s.median():.3f}, max {s.max():.3f}")

CONV_WARN_THRESHOLD = 0.80
if conv_rate < CONV_WARN_THRESHOLD:
    print(f"\n⚠ WARNING: only {100*conv_rate:.0f}% of synthetic fits converged "
          f"(threshold {100*CONV_WARN_THRESHOLD:.0f}%). Interpret bias/sd with "
          f"caution — unconverged fits contaminate the summary.")

# ── 4. Summary tables: converged-only AND all-ok ──
def summarise(subset: pl.DataFrame, tag: str):
    rows = []
    for p, true in [("beta", truth["beta"]), ("gamma", truth["gamma"])]:
        vals = subset[p].to_numpy()
        if len(vals) == 0:
            continue
        rows.append({
            "subset": tag, "param": p, "n": len(vals), "truth": true,
            "mean_est": float(vals.mean()),
            "bias": float(vals.mean() - true),
            "sd": float(vals.std()),
            "rmse": float(np.sqrt(np.mean((vals - true)**2))),
            "q025": float(np.quantile(vals, 0.025)),
            "q975": float(np.quantile(vals, 0.975)),
        })
    return rows

srows = summarise(conv, "converged_only") + summarise(ok, "all_ok")
summary_df = pl.DataFrame(srows)
summary_df.write_csv(OUT/"sbc_summary.tsv", separator="\t")
print("\nSBC summary:")
print(summary_df)

# ── 5. Plot: param histograms + Rhat distributions ──
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
               label=f"mean(converged) = {vals_conv.mean():.3f}")
    ax.set_xlabel(f"{'β̂' if p=='beta' else 'γ̂'}")
    ax.set_title(f"{p}  "
                 f"(bias = {vals_conv.mean()-truth[p]:+.3f}, "
                 f"sd = {vals_conv.std():.3f}, "
                 f"95% CI: [{np.quantile(vals_conv,0.025):.3f}, "
                 f"{np.quantile(vals_conv,0.975):.3f}])",
                 fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top","right"]].set_visible(False)

# Rhat distribution (both params on same axis)
ax = axes[1,0]
rh_all = np.concatenate([ok["rhat_beta"].to_numpy(), ok["rhat_gamma"].to_numpy()])
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

# Convergence summary panel (text)
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
    f"    threshold: all refine Rhat < 1.10 must hold",
    f"    warn if converged rate < {100*CONV_WARN_THRESHOLD:.0f}%",
    "",
    "    " + ("✓ PASS" if conv_rate >= CONV_WARN_THRESHOLD else "⚠ WARN"),
]
for i, ln in enumerate(lines):
    ax.text(0.02, 0.95 - i*0.08, ln, transform=ax.transAxes,
            fontfamily="monospace", fontsize=10)

fig.suptitle(f"Synthetic calibration at real-data MLE "
             f"({N_SYNTH} synthetic datasets from θ* = "
             f"β={truth['beta']:.3f}, γ={truth['gamma']:.3f})")
fig.tight_layout()
fig.savefig(FIG/"sbc.png", bbox_inches="tight"); plt.close(fig)
print(f"wrote {FIG/'sbc.png'}")
