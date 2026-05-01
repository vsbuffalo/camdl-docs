# He et al. (2010) Measles Vignette — Chapter Plan

## The pedagogical gap

The reader has just fit a 2-parameter SIR to 14 observations in a closed
population (boarding school). He et al. is a 12-parameter SEIR with
time-varying forcing, open population dynamics, cohort pulses,
overdispersed transmission, and 1096 weekly observations over 21 years.
That's at least 6 new concepts. We introduce them one at a time.

---

## Chapter 1: "The Model" — from SIR to SEIR with forcing

Build up from the boarding-school SIR the reader already knows. Each
section adds one mechanism and shows why it's needed.

1. **Why the boarding-school SIR can't do measles.** 14 days vs 21
   years. Closed population vs births/deaths. Single outbreak vs
   recurrent biennial epidemics. Motivate each extension before
   introducing it.

2. **Open the population — SIR + demography.** Add births and deaths
   (μ parameter). Show in camdl DSL side-by-side with the closed SIR.
   Simulate: the epidemic now recurs because the susceptible pool
   refills. This is the single biggest conceptual step — from
   "outbreak" to "endemic dynamics."

3. **Add a latent period — SEIR + demography.** Introduce the E
   compartment and σ parameter. Why it matters for measles (8–13 day
   incubation). Show the effect on epidemic timing and shape vs SIR.

4. **School-term forcing.** The central mechanism of He et al.:
   measles transmission tracks the school calendar. The `periodic`
   forcing block. Show what the forcing function looks like. Simulate:
   biennial epidemics emerge from the interaction of seasonality and
   susceptible depletion.

5. **Time-varying covariates.** `interpolated` forcing for pop(t) and
   birthrate(t). The `forcing {}` block. Why a 21-year model needs
   demography that tracks real population changes.

6. **Cohort births.** The `events {}` block. Children enter school
   once per year on day 251. The `balance {}` block to keep
   S+E+I+R = pop(t).

7. **Overdispersed transmission.** `overdispersed()` primitive (Gamma
   noise on the force of infection). Why deterministic rates can't
   explain the observed variance — connect back to the process-noise
   analysis from the model comparison chapter.

8. **The observation model.** Incidence vs prevalence (observe I→R
   recoveries, not I). Reporting fraction ρ. The discretized Normal
   (He et al. eq. 6). Weekly aggregation.

9. **The full model.** Show the complete `he2010_london.camdl` — by
   now every line is recognizable. Forward-simulate at published
   parameters, overlay on London measles data. The reader sees the
   model generates realistic biennial epidemics.

**Source material:** `he2010_london.camdl`, forward simulation from
`he2010-forward/`, data files in `data/`.

---

## Chapter 2: "Fitting" — IF2 on a real-world model

The reader already knows IF2 from the boarding school. Now show what
changes at scale.

1. **The fit config** — 12 parameters, wider bounds, more particles,
   longer scouts. Walk through `fit_he2010.toml`.

2. **IF2 scout** — convergence traces. 8 chains, 1096-week traces.
   Which parameters converge? Which don't? (s0 multimodal, iota
   plateaus.)

3. **Convergence diagnostics** — Â per parameter. Real
   non-convergence on s0 (Â=3.2 on full 21-year).

4. **Profile likelihoods** — the s0 ridge, the iota plateau. Connect
   back to sloppy-axis analysis on k from the boarding school.

5. **The MLE** — compare to He et al.'s published values.

6. **Observation model comparison** — overdispersion ablation study.
   Connect back to Poisson vs NegBin analysis from boarding school.

**Source material:** `he2010-inference/` results,
`he2010-diagnostics/` outputs, `he2010-model-comparison/index.qmd`.

---

## Chapter 3: "Posterior and Identifiability" (blocked on upstream PGAS)

1. **Why the MLE isn't enough** — R0-alpha ridge, sigma-gamma
   compensation. Profile likelihoods showing structural
   non-identifiability.

2. **PMMH and why it fails here** — PF variance problem (sd=173
   nats). Brief, honest, instructive.

3. **PGAS (when available)** — complete-data likelihood breaks the
   marginal degeneracy. Placeholder.

**Source material:** `he2010-pmmh/pmmh-report.md`, identifiability
analysis from diagnostics.
