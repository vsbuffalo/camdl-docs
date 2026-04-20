"""Refresh model-comparison figure + SUMMARY using the 1D-profile MLE."""
import sys
from pathlib import Path
import numpy as np, matplotlib.pyplot as plt, polars as pl

HERE = Path(__file__).parent
BOOK = HERE.parent
sys.path.insert(0, str(BOOK))
from styles import apply_rc

apply_rc()

# Profile-based MLE (converged: each cell is a fit where σ_se is fixed
# and β,γ are IF2-optimised; refine convergence is achieved per cell)
df1d = pl.read_csv(HERE/"output/profile_sigma_1d.tsv",
                   separator="\t", comment_prefix="#")
ll_proc = float(df1d["max_loglik"].max())
row = df1d.sort("max_loglik", descending=True).row(0, named=True)
sigma_star = row["sigma_se"]; beta_star = row["beta"]; gamma_star = row["gamma"]

ll_poi = -61.1391851
ll_nb  = -58.0857824

models = [
    ("Poisson\nfixed IC",              ll_poi, 2),
    ("NegBin\nfixed IC",               ll_nb,  3),
    ("Poisson + β process-noise\nσ_se profile MLE", ll_proc, 3),
]
names = [m[0] for m in models]
lls   = np.array([m[1] for m in models])
ks    = np.array([m[2] for m in models])
aics  = 2*ks - 2*lls
daic  = aics - aics.min()

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
colors = ["#bdc3c7", "#3498db", "#e67e22"]

ax = axes[0]
ax.bar(range(3), lls, color=colors)
ax.set_xticks(range(3)); ax.set_xticklabels(names, fontsize=9)
ax.set_ylabel("log-likelihood")
ax.set_title("Fit log-likelihood on boarding-school data", fontsize=10)
ax.set_ylim(min(lls)-1, max(lls)+0.6)
for i, v in enumerate(lls):
    ax.text(i, v+0.1, f"{v:.2f}", ha="center", fontsize=9)
ax.spines[["top","right"]].set_visible(False)

ax = axes[1]
ax.bar(range(3), daic, color=colors)
ax.set_xticks(range(3)); ax.set_xticklabels(names, fontsize=9)
ax.set_ylabel("ΔAIC (from best)")
ax.set_title("ΔAIC — lower is better", fontsize=10)
for i, v in enumerate(daic):
    ax.text(i, v+0.08, f"{v:.2f}", ha="center", fontsize=9)
ax.spines[["top","right"]].set_visible(False)

fig.tight_layout()
fig.savefig(HERE/"figures/model_comparison.png", bbox_inches="tight")
plt.close(fig)

# Rewrite SUMMARY
(HERE/"SUMMARY.md").write_text(f"""# Process-noise-on-β test fit — summary

**Model:** `boarding_school_sir_procnoise.camdl` — SIR with
`overdispersed(β·S·I/N, σ_se)` on the infection transition (gamma-noise
on β, He et al. 2010), Poisson observation model.

**Fit status (standards-compliant):**
- Joint IF2 fit fails the built-in Rhat ≤ 1.10 gate because of a real
  β–σ_se ridge (corr ≈ 0.88; Rhat σ_se=4.29, β=2.32).
- The honest MLE comes from a **1D σ_se profile** (σ_se fixed per cell,
  β and γ IF2-optimised, refine converges at each cell).

## Converged comparison table

| Model                              | k | log-lik  | AIC    | ΔAIC |
|------------------------------------|---|----------|--------|------|
| Poisson, fixed IC                  | 2 | {ll_poi:.2f} | {2*2 - 2*ll_poi:.2f} | {(2*2 - 2*ll_poi) - aics.min():.2f} |
| NegBin, fixed IC                   | 3 | {ll_nb:.2f} | {2*3 - 2*ll_nb:.2f} | {(2*3 - 2*ll_nb) - aics.min():.2f} |
| Poisson + β process-noise (σ_se)   | 3 | **{ll_proc:.2f}** | {2*3 - 2*ll_proc:.2f} | **{(2*3 - 2*ll_proc) - aics.min():.2f}** |

Profile MLE:  σ_se = {sigma_star:.3f},  β = {beta_star:.3f},  γ = {gamma_star:.3f}

**95% Wilks CI on σ_se:** [0.050, 0.550] — a factor of 10, but bounded
*away from 0*. Process noise is identified but weakly constrained.

## Headline

Both extra-variance models beat plain Poisson by ~3–5 nats. Process
noise on β beats NegBin by Δll ≈ {ll_proc - ll_nb:.2f} / ΔAIC ≈ {2*3 - 2*ll_proc - (2*3 - 2*ll_nb):.2f}. Neither is a
landslide — within the 95% Wilks band the two models overlap — but the
data prefers process noise.

## Figures (saved in `figures/`)

- `profile_sigma_1d.png` — 1D σ_se profile with Wilks CI + β̂, γ̂ traces
- `profile_beta_sigma.png` — 2D β × σ_se profile (the ridge)
- `model_comparison.png` — ll + ΔAIC bars
- `procnoise_scout_traces.png` — scout traces (non-converged, for
  context — the ridge is visible in σ_se drift)
- `procnoise_scout_rank.png`, `procnoise_scout_pairplot.png`,
  `procnoise_ppc.png` — earlier scout-based diagnostics

## What this means for the chapter

The obs-model-vs-process-noise comparison is a **genuine** comparison on
this data. Neither collapses; neither wins by a landslide. That's the
pedagogical gold: the reader sees two mechanistically different stories
(observer noisy / underlying dynamics noisy) fitting nearly equally well
and has to reason about which story the science supports.
""")
print(f"procnoise profile MLE: ll = {ll_proc:.2f}")
print(f"  σ_se = {sigma_star:.3f}, β = {beta_star:.3f}, γ = {gamma_star:.3f}")
print(f"  ΔAIC vs NegBin:  {(2*3 - 2*ll_nb) - aics.min():.2f}")
print(f"  ΔAIC vs Poisson: {(2*2 - 2*ll_poi) - aics.min():.2f}")
print(f"\nwrote figures/model_comparison.png + SUMMARY.md")
