"""Pre-fit pipeline validation: simulate from realistic-flu guesstimate
parameters, fit with the same IF2 pipeline, check recovery.

Truth:
  β = 1.20 /day        # transmission rate, per day per infective
  γ = 0.33 /day        # recovery rate; mean infectious period ≈ 3 days
  I(0) = 2             # index cases
  N0 = 763             # known school size
  R₀ = β/γ ≈ 3.6       # pandemic-flu-in-children ballpark (literature range 2–4)

These differ substantially from the real-data MLE (β=1.91, γ=0.66). If the
IF2 pipeline recovers values close to truth here, the pipeline is sound;
if it drifts to the real-data basin, it's broken or data is unidentified.
"""
import subprocess, shlex, time
from pathlib import Path
import polars as pl, numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"; FIG.mkdir(exist_ok=True)
OUT  = ROOT/"output";  OUT.mkdir(exist_ok=True)

GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

TRUTH = {"beta": 1.20, "gamma": 0.33, "I0": 5, "N0": 763}  # I0 matches fit's fixed value
TRUTH_SEED = 7001

# 1. Write a truth params file
truth_toml = OUT/"synthetic_truth.toml"
truth_toml.write_text("\n".join(f"{k} = {v}" for k, v in TRUTH.items()) + "\n")

# 2. Generate one synthetic dataset under truth using the "synth" model:
# observations every 2 days over 28 days → 15 points spanning the whole
# outbreak (R₀ = 3.6 gives a ~25-day outbreak; daily obs over only 14 days
# would miss the decline).
synth_tsv = OUT/"synthetic_data.tsv"
subprocess.run(shlex.split(
    f"camdl simulate boarding_school_sir_synth.ir.json --params {truth_toml} "
    f"--replicates 1 --seed {TRUTH_SEED} --scenario baseline "
    f"--backend chain_binomial --dt 1.0 --obs {synth_tsv}"),
    check=True, capture_output=True, cwd=str(ROOT))
synth = pl.read_csv(synth_tsv, separator="\t")
# single-replicate output omits the replicate column already; rewrite as-is
synth.write_csv(synth_tsv, separator="\t")
print(f"synthetic peak: day {int(synth.filter(pl.col('in_bed') == synth['in_bed'].max())['time'][0])} "
      f"at {int(synth['in_bed'].max())} boys; total bed-days {int(synth['in_bed'].sum())}")

# 3. Template a fit.toml that points at the synthetic dataset and synth model
fit_template = (ROOT/"fit.toml").read_text()
fit_synth = fit_template.replace('camdl = "boarding_school_sir.ir.json"',
                                 'camdl = "boarding_school_sir_synth.ir.json"')
fit_synth = fit_synth.replace('in_bed = "data/in_bed.tsv"',
                              f'in_bed = "{synth_tsv.resolve()}"')
fit_synth = fit_synth.replace('output_dir = "results/sir"',
                              'output_dir = "synthetic_recovery"')
(ROOT/"fit_synthetic.toml").write_text(fit_synth)

# 4. Run fit (--force so cache doesn't hide the log we need for Rhat)
print("running IF2 fit on synthetic data...")
t0 = time.time()
log_path = OUT/"synthetic_fit.log"
with log_path.open("w") as f:
    subprocess.run(["camdl", "fit", "run", "fit_synthetic.toml", "--seed", "42", "--force"],
                   stdout=f, stderr=subprocess.STDOUT, cwd=str(ROOT), check=True)
print(f"fit done in {time.time()-t0:.1f}s")

# 5. Read MLE
fit_dir = next((ROOT/"results/fits").glob("fit_synthetic-*/real/fit_42/refine"))
mle = {}
for ln in (fit_dir/"mle_params.toml").read_text().splitlines():
    line = ln.split("#",1)[0].strip()
    if line.startswith("[") or "=" not in line: continue
    k, v = line.split("=", 1)
    try: mle[k.strip()] = float(v.strip())
    except: pass

# 6. Parse Rhat from log — only the final (refine) Rhat block, so scout's
# loose Rhats don't overwrite refine's.
import re
log_text = log_path.read_text()
# take everything after the last "── stage: refine" header
refine_section = log_text.split("── stage: refine")[-1]
rhats = {}
for ln in refine_section.splitlines():
    m = re.search(r"^\s*(\w+)\s+Rhat=(\d+\.\d+)", ln)
    if m: rhats[m.group(1)] = float(m.group(2))

print(f"\n=== Recovery check ===")
print(f"{'param':>8} {'truth':>10} {'estimate':>10} {'abs err':>10} {'% err':>8} {'refine Rhat':>12}")
for p in ["beta", "gamma"]:
    t, e = TRUTH[p], mle[p]
    print(f"{p:>8} {t:>10.3f} {e:>10.3f} {abs(e-t):>10.3f} {100*abs(e-t)/t:>7.1f}% {rhats.get(p, float('nan')):>12.3f}")

# 7. PPC at recovered MLE against synthetic data — use the synth model
# so the 28-day outbreak plays out fully and obs grid matches the data.
ppc_tsv = OUT/"synthetic_ppc.tsv"
subprocess.run(shlex.split(
    f"camdl simulate boarding_school_sir_synth.ir.json "
    f"--params {fit_dir}/mle_params.toml "
    f"--replicates 200 --seed 12345 --scenario baseline "
    f"--backend chain_binomial --dt 1.0 --obs {ppc_tsv}"),
    check=True, capture_output=True, cwd=str(ROOT))
ppc = pl.read_csv(ppc_tsv, separator="\t")

# 8. Plot: truth trajectory, synthetic data points, recovered PPC envelope
fig, ax = plt.subplots(figsize=(8, 4.5))
q = ppc.group_by("time").agg(
    pl.col("in_bed").quantile(0.05).alias("lo"),
    pl.col("in_bed").quantile(0.5).alias("med"),
    pl.col("in_bed").quantile(0.95).alias("hi"),
).sort("time")
ax.fill_between(q["time"], q["lo"], q["hi"], color="#2980b9", alpha=0.22, linewidth=0, edgecolor="none",
                label=f"PPC at recovered MLE\n(β̂={mle['beta']:.3f}, γ̂={mle['gamma']:.3f})")
ax.plot(q["time"], q["med"], color="#2980b9", linewidth=1.6)
ax.scatter(synth["time"], synth["in_bed"], color="black", s=36, zorder=5,
           label=f"synthetic data (truth β={TRUTH['beta']}, γ={TRUTH['gamma']}, I(0)={TRUTH['I0']})")
ax.set_xlabel("Day"); ax.set_ylabel("Boys in bed")
ax.set_title("Pipeline validation — recover guesstimate flu parameters from synthetic data")
ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.spines[["top","right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG/"synthetic_recovery.png", bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {FIG/'synthetic_recovery.png'}")

# Save a summary row for the report
summary = pl.DataFrame([
    {"param": p, "truth": TRUTH[p], "estimate": mle[p],
     "abs_err": abs(mle[p]-TRUTH[p]),
     "pct_err": 100*abs(mle[p]-TRUTH[p])/TRUTH[p],
     "rhat": rhats.get(p, float("nan"))}
    for p in ["beta", "gamma"]
])
summary.write_csv(OUT/"synthetic_recovery.tsv", separator="\t")
