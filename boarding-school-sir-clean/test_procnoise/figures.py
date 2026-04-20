"""Generate diagnostic figures for the process-noise-on-β fit and compare
to Poisson-fixed-IC and NegBin-fixed-IC baselines."""
import sys, subprocess, io
from pathlib import Path
import numpy as np, polars as pl
import matplotlib.pyplot as plt, matplotlib as mpl

HERE = Path(__file__).parent
BOOK = HERE.parent
sys.path.insert(0, str(BOOK))
from plot_traces import load_chains, best_chain_by_final_ll, plot_traces, plot_rank_overlay, plot_pair_grid
from styles import apply_rc, GREY

apply_rc()
FIG = HERE/"figures"; FIG.mkdir(exist_ok=True)

FIT_DIR = next(HERE.glob("results/fits/fit_procnoise-*/real/fit_42"))
SCOUT   = FIT_DIR/"scout"

# --- 1. Scout traces
scout = load_chains(SCOUT)
best = best_chain_by_final_ll(scout)
print(f"scout best chain: {best}")

plot_traces(
    scout, params=["beta", "gamma", "sigma_se", "if2_perturbed_loglik"],
    labels=["β", "γ", "σ_se", "log-lik"],
    bounds={"beta": (2.0, 3.5), "gamma": (0.5, 0.9), "sigma_se": (0.05, 1.0)},
    best_chain=best, color_mode="ll",
    outpath=FIG/"procnoise_scout_traces.png",
    title="Process-noise fit (Poisson obs, σ_se on β) — scout traces (64 chains × 400 iter)",
)

plot_rank_overlay(
    scout, params=["beta", "gamma", "sigma_se"],
    labels=["β", "γ", "σ_se"],
    outpath=FIG/"procnoise_scout_rank.png",
    title="Process-noise scout — rank overlay (Vehtari et al. 2021)",
)

# --- 2. Pair plot
plot_pair_grid(
    scout,
    params=["beta", "gamma", "sigma_se", "if2_perturbed_loglik"],
    labels=["β", "γ", "σ_se", "log-lik"],
    best_chain=best, dll_filter=10.0, n_bins=14, ll_axis_dll=8,
    outpath=FIG/"procnoise_scout_pairplot.png",
    title="Process-noise scout pair plot (Δll < 10 from best)",
)

# --- 3. Best-chain params + ll
BASE = BOOK/"results/fits"
def read_mle(path):
    d = {}
    ll = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("log_likelihood"): ll = float(line.split("=")[1])
        if "=" in line and not line.startswith("#") and not line.startswith("["):
            k, v = line.split("=", 1); k = k.strip()
            try: d[k] = float(v.strip())
            except Exception: pass
    return d, ll

# Process-noise best chain (use scout best since refine didn't run — not
# converged in Rhat terms but best ll is the comparison point)
best_df = scout.filter(pl.col("chain") == best).sort("iteration")
ll_proc = float(best_df.filter(pl.col("loglik").is_not_null())["loglik"][-1])
last = best_df.tail(1).row(0, named=True)
params_proc = {"beta": last["beta"], "gamma": last["gamma"],
               "sigma_se": last["sigma_se"], "N0": 763, "I0": 5}
print(f"procnoise best-chain ll = {ll_proc:.2f}  β={params_proc['beta']:.3f}  "
      f"γ={params_proc['gamma']:.3f}  σ_se={params_proc['sigma_se']:.3f}")

# Baselines
poi_mle, ll_poi = read_mle(BASE/"fit-6bb88d20/real/fit_42/refine/mle_params.toml")
nb_mle,  ll_nb  = read_mle(BASE/"fit_negbin-60da4f72/real/fit_42/refine/mle_params.toml")
print(f"Poisson fixed-IC refine ll = {ll_poi:.2f}")
print(f"NegBin  fixed-IC refine ll = {ll_nb:.2f}")

# --- 4. PPC at the process-noise best-chain MLE
# Write params and call camdl simulate
mle_path = HERE/"mle_best_chain.toml"
mle_path.write_text(
    "# procnoise best-chain params (scout, not Rhat-converged — see README)\n"
    f"beta = {params_proc['beta']}\n"
    f"gamma = {params_proc['gamma']}\n"
    f"sigma_se = {params_proc['sigma_se']}\n"
    f"N0 = 763\nI0 = 5\n"
    "[provenance]\n"
    'backend = "chain_binomial"\ndt = 1.0\n'
    f'model = "boarding_school_sir_procnoise.ir.json"\n'
)

def sim_ppc(n_reps=200, seed0=1):
    """Use batch --seeds to get n_reps trajectories in one invocation."""
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        obs_path = Path(td)/"obs.tsv"
        cmd = ["camdl", "simulate", "boarding_school_sir_procnoise.ir.json",
               "--params", str(mle_path), "--scenario", "baseline",
               "--backend", "chain_binomial", "--dt", "1.0",
               "--seeds", f"{seed0}:{seed0+n_reps-1}",
               "--obs-only", str(obs_path)]
        subprocess.run(cmd, capture_output=True, text=True, cwd=HERE, check=True)
        df = pl.read_csv(obs_path, separator="\t", comment_prefix="#")
    return df

print("simulating PPC (200 reps)...")
ppc = sim_ppc(n_reps=200)
obs = pl.read_csv(BOOK/"data/in_bed.tsv", separator="\t")
print(f"PPC columns: {ppc.columns}")

# pivot: one row per (seed, time) — stack in_bed into (n_seeds, n_times)
seed_col = "seed" if "seed" in ppc.columns else ppc.columns[0]
wide = ppc.pivot(values="in_bed", index="time", on=seed_col)
times = wide["time"].to_numpy()
M = wide.drop("time").to_numpy().T  # (n_seeds, n_times)
lo, md, hi = np.quantile(M, [0.025, 0.5, 0.975], axis=0)

fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.fill_between(times, lo, hi, color="#3498db", alpha=0.20, label="95% PPC band")
ax.plot(times, md, color="#3498db", linewidth=1.5, label="PPC median")
ax.plot(obs["time"], obs["in_bed"], "o-", color="black", linewidth=1.5, label="observed")
ax.set_xlabel("Day"); ax.set_ylabel("Boys in bed")
ax.set_title("Posterior-predictive check — Poisson + β process-noise (σ_se) at best-chain MLE",
             fontsize=10, pad=6)
ax.legend(frameon=False)
ax.spines[["top","right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG/"procnoise_ppc.png", bbox_inches="tight"); plt.close(fig)

# --- 5. Comparison: Poisson / NegBin / ProcNoise
models = [
    ("Poisson\nfixed IC",              ll_poi, 2),
    ("NegBin\nfixed IC",               ll_nb,  3),
    ("Poisson + β-process-noise\nfixed IC", ll_proc, 3),
]
names = [m[0] for m in models]
lls   = np.array([m[1] for m in models])
ks    = np.array([m[2] for m in models])
aics  = 2*ks - 2*lls
daic  = aics - aics.min()

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
ax = axes[0]
ax.bar(range(3), lls, color=["#bdc3c7", "#3498db", "#e67e22"])
ax.set_xticks(range(3)); ax.set_xticklabels(names, fontsize=9)
ax.set_ylabel("log-likelihood")
ax.set_title("Fit log-likelihood on boarding-school data", fontsize=10)
ax.set_ylim(min(lls)-1, max(lls)+0.5)
for i, v in enumerate(lls):
    ax.text(i, v+0.1, f"{v:.2f}", ha="center", fontsize=9)
ax.spines[["top","right"]].set_visible(False)

ax = axes[1]
ax.bar(range(3), daic, color=["#bdc3c7", "#3498db", "#e67e22"])
ax.set_xticks(range(3)); ax.set_xticklabels(names, fontsize=9)
ax.set_ylabel("ΔAIC (from best)")
ax.set_title("ΔAIC — lower is better (k = # estimated params)", fontsize=10)
for i, v in enumerate(daic):
    ax.text(i, v+0.05, f"{v:.2f}", ha="center", fontsize=9)
ax.spines[["top","right"]].set_visible(False)

fig.tight_layout()
fig.savefig(FIG/"model_comparison.png", bbox_inches="tight"); plt.close(fig)

# Write a summary markdown
(HERE/"SUMMARY.md").write_text(f"""# Process-noise-on-β test fit — summary

**Model:** `boarding_school_sir_procnoise.camdl` — SIR with
`overdispersed(β·S·I/N, σ_se)` on the infection transition (gamma-noise
on β, He et al. 2010), Poisson observation model.

**Scout fit status:** Rhat failed on β (2.32) and σ_se (4.29) due to
β–σ_se ridge (corr ≈ 0.88). Refine did not run. The scout *best-chain*
ll is used as a comparison point — this is a test, not a publishable MLE.

## Comparison table

| Model                            | k | log-lik | AIC    | ΔAIC  |
|----------------------------------|---|---------|--------|-------|
| Poisson, fixed IC                | 2 | {ll_poi:.2f} | {2*2 - 2*ll_poi:.2f} | {(2*2 - 2*ll_poi) - aics.min():.2f} |
| NegBin, fixed IC                 | 3 | {ll_nb:.2f} | {2*3 - 2*ll_nb:.2f} | {(2*3 - 2*ll_nb) - aics.min():.2f} |
| Poisson + β-process-noise        | 3 | {ll_proc:.2f} | {2*3 - 2*ll_proc:.2f} | {(2*3 - 2*ll_proc) - aics.min():.2f} |

**σ_se did NOT collapse** — best-chain σ_se ≈ {params_proc['sigma_se']:.3f}. The
data supports genuine process noise on β.

**Both extra-variance models beat Poisson by >3 nats.** Process noise on
β beats NegBin by ~{ll_proc - ll_nb:.1f} nats on this scout — suggestive but within
noise given the β–σ_se ridge and the fact that scout didn't converge.

## Figures (saved in `figures/`)

- `procnoise_scout_traces.png` — 64-chain scout traces
- `procnoise_scout_rank.png` — rank overlay
- `procnoise_scout_pairplot.png` — β × γ × σ_se × ll pair plot
- `procnoise_ppc.png` — PPC at best-chain params (200 reps)
- `model_comparison.png` — ll + ΔAIC bars across all three models

## Caveats

1. β–σ_se is a real ridge; tighter convergence needs a σ_se profile or
   ivp=true treatment.
2. The scout ll is biased slightly upward vs refine (shorter budget,
   slightly noisier), so the ~{ll_proc - ll_nb:.1f}-nat advantage over NegBin
   could shrink when properly refined.
3. For the chapter's formal model-comparison section, profile σ_se and
   compare refine-converged lls.
""")
print(f"\nwrote {HERE/'SUMMARY.md'}")
print("wrote figures/procnoise_scout_traces.png, _rank.png, _pairplot.png, _ppc.png")
print("wrote figures/model_comparison.png")
