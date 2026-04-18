# Fitting Chapter Rewrite Plan

**Date:** 2026-04-18
**Status:** proposal (revised after review)
**Scope:** Complete rewrite of `guide/fitting.qmd`
**Target:** ~8 pages of content + referenced vignettes

## Thesis

At K=14 daily observations from a stochastic SIR epidemic, finite-sample
MLE bias doesn't disappear — it redistributes across parameters depending
on the observation model. You choose which parameter to trust.

## Revised structure (8 pages)

### 1. Opening: the problem (1 page)

Latent states, intractable likelihood, particle filter, IF2. One paragraph
each.

End with the methodological frame: "Before fitting real data, we
generated synthetic datasets at plausible parameters and fit those —
partly to sanity-check the pipeline, partly because we wanted to know
what 14 observations could tell us before we asked them to tell us
anything."

### 2. The model and the first fit (1.5 pages)

Present the boarding school SIR. Show the `.camdl` model file with
NegBin observation model — the "obvious" choice. Show `fit.toml`.

**Do NOT introduce Poisson/Binomial yet.** The reader has no reason to
care about alternatives until they've seen the problem. Introduce only
NegBin: "the standard choice for overdispersed count data."

### 3. Synthetic validation: discovering the bias (2.5 pages)

#### 3a. Generate, fit, discover

Generate 30 synthetic datasets at known (β=2.20, γ=0.58, k=50). Fit each.

Before showing results: "If this were a normal stats problem, we'd
expect IF2 to converge to the truth on average. It doesn't — and
that's not a bug."

Show: parameter scatter (30 dots, truth marked), R₀ scatter showing
~+50% bias. Use scatter-plot R₀ (compute per-dataset, take mean) not
ratio-of-means — more honest and consistent throughout.

#### 3b. Does more data fix it?

R₀ bias vs observation density (K=14, 28, 56, 140). Two-panel: β and R₀.
Bias collapses by K=28. This is a small-data problem, not a pipeline bug.

#### 3c. The bias budget

Table: for each parameter at K=14, report bias, variance, MSE, and
cross-correlations Cor(β̂, k̂), Cor(β̂, γ̂). The Cor(β̂, k̂) is the
smoking gun that k is compensating for β.

### 4. The observation model tradeoff (2 pages)

"Maybe it's the k parameter." Now introduce Poisson and Binomial as
alternatives — the reader has a reason to want them.

#### 4a. Switch to Poisson, rerun SBC

β bias drops from +59% to +3%. But γ bias jumps from +9% to +41%.
"The bias didn't disappear. It moved."

"Notice that neither obs model has both parameters in the right place.
There is no row in the SBC table where bias is small for everything."

#### 4b. Three-model comparison figure

Side-by-side scatter: NegBin, Poisson, Binomial. Two panels: β and γ.
The headline figure.

#### 4c. The bias-sink principle

Table: bias budget for all three models including R₀.

"The R₀ you'd report to a public health agency depends on which
observation model you picked — NegBin says ~5.7, Poisson says ~2.8."

One-line callout: "Neither estimate is 'wrong' in the usual sense —
both MLEs converged to the likelihood maximum under their respective
models. The choice of observation noise structure silently determines
the answer."

### 5. Fitting the real data (1.5 pages)

Bridge: "With the pipeline's limits characterized, we turn to the
real 1978 outbreak — bringing our calibrated skepticism with us."

Three-panel figure: NegBin, Poisson, Binomial on real data.
Parameters and R₀ differ. The reader already knows why.

"The three 'fits' are not three estimates of the same thing — they're
three different inference problems posing as one."

Note the peak-timing mismatch (model peaks day 3-4, data peaks day
5-6). Compact discussion: run 200 replicates at Poisson MLE to check
stochasticity contribution, fit SEIR-Poisson (fix 1/σ=1.5 days from
clinical flu) to check misspecification. Report in one paragraph with
one comparison figure. State the conclusion: "the peak-height gap is
primarily misspecification (missing latent period), with a minor
contribution from stochasticity."

### 6. What's next (0.5 pages)

- **Bootstrap bias correction** — one paragraph summarizing the
  result (helps at K=14, neutral at K=56), pointing to the vignette
  for details
- **Bayesian inference with informed priors** — one paragraph
  summarizing the γ-prior result, pointing to the vignette
- **SIBCR model** — the observation *mapping* question (bed count =
  post-infectious), distinct from the observation *noise* question
  explored in this chapter. Concrete forward reference to the vignette.
- **He et al. measles** — same pipeline, much harder model

## Vignettes (pulled from chapter)

### Vignette: Bootstrap bias correction
Sections 9 from old plan. Full protocol, MSE comparison across K,
per-parameter analysis. Referenced from section 6.

### Vignette: Bayesian inference on small stochastic data
Sections 8 from old plan. PGAS with informed priors, prior vs
posterior plots, pairplots, when-is-an-informative-prior-justified
callout. Referenced from section 6.

### Vignette: Model comparison (SIBCR and observation mapping)
The SIBCR/SEIBCR analysis. Avilov/Tverskoi/Wearing papers. The
observation *mapping* question vs the observation *noise* question.
Referenced from section 6.

## Experiments

### Already done (reuse)

- [x] NegBin SBC at K=14 (30 datasets, well-specified)
- [x] NegBin SBC at K=28, K=56, K=140
- [x] Poisson SBC at K=14 (20 datasets, well-specified)
- [x] Binomial SBC at K=14 (20 datasets, well-specified)
- [x] Real data fits: NegBin, Poisson, Binomial
- [x] PGAS with informed priors (→ vignette)
- [x] Bootstrap at K=14, K=28, K=56 (→ vignette)

### Must have (new)

- [ ] **Bias budget table**: Cor(β̂, k̂), Cor(β̂, γ̂) across NegBin
      SBC datasets. Also SD(k̂)/mean(k̂). Anchors sections 3c and 4c.
- [ ] **Peak decomposition**: 200 replicates at Poisson MLE on real
      data (fraction with peak ≥ 298), SEIR-Poisson fit to real data.
      Anchors section 5.
- [ ] **SEIR model**: `boarding_school_seir_poisson.camdl`, fix
      1/σ=1.5 days, estimate β, γ with Poisson obs. Fit real data.
- [ ] **R₀ reconciliation**: recompute all R₀ biases using per-dataset
      R₀ = β̂/γ̂ consistently (not ratio of means).

### Nice to have (cut if time-limited)

- [ ] Profile likelihoods for three obs models. Reinforces section 4
      but may not earn its keep in a Getting Started chapter. Consider
      for vignette instead.

## Key sentences to land

Section 3a: "If this were a normal stats problem, we'd expect IF2 to
converge to the truth on average. It doesn't — and that's not a bug."

Section 4a: "Notice that neither obs model has both parameters in the
right place. There is no row in the SBC table where bias is small for
everything."

Section 4c: "Neither estimate is 'wrong' in the usual sense — both
MLEs converged to the likelihood maximum under their respective models.
The choice of observation noise structure silently determines the
answer."

Section 5: "The three 'fits' are not three estimates of the same
thing — they're three different inference problems posing as one."

## Approach

Stop perfecting the outline. Run the bias budget table and peak
decomposition experiments, then draft. The narrative follows from
what the experiments show. If SEIR-Poisson fits the peak well,
section 5 ends with misspecification confirmed. If it doesn't, section
5 ends with "even the right obs model + right compartmental structure
leaves residual bias that's a true small-K problem." Either ending
works. Just write it.
