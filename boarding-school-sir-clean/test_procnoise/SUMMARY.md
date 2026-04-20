# Process-noise-on-β test fit — summary

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
| Poisson, fixed IC                  | 2 | -61.14 | 126.28 | 7.50 |
| NegBin, fixed IC                   | 3 | -58.09 | 122.17 | 3.39 |
| Poisson + β process-noise (σ_se)   | 3 | **-56.39** | 118.78 | **0.00** |

Profile MLE:  σ_se = 0.120,  β = 2.294,  γ = 0.698

**95% Wilks CI on σ_se:** [0.050, 0.550] — a factor of 10, but bounded
*away from 0*. Process noise is identified but weakly constrained.

## Headline

Both extra-variance models beat plain Poisson by ~3–5 nats. Process
noise on β beats NegBin by Δll ≈ 1.70 / ΔAIC ≈ -3.39. Neither is a
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
