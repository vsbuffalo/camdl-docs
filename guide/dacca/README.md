# Dhaka cholera — model-comparison case study

This directory holds the data, source paper, and (forthcoming)
`.camdl` models for a model-comparison chapter built on the historic
Dhaka district cholera-mortality series — the canonical pomp
model-comparison dataset, from King, Ionides, Pascual & Bouma,
*"Inapparent infections and cholera dynamics,"* **Nature 454:877–880
(2008)** ([doi:10.1038/nature07084](https://doi.org/10.1038/nature07084);
PDF in `king_et_al_2008_nature.pdf`).

It will likely **replace the boarding-school flu** as the running
example in `guide/fitting/model_comparison.qmd` (the boarding-school
example ended up with *the simple model wins* — fine as a one-off,
weak as a model-comparison teaching case; that chapter is currently
`draft: true`). Whether the boarding-school flu stays as a short
"simple-can-win" aside inside the rewritten chapter, or moves out
entirely, is TBD.

## Data

| file | what | shape |
|---|---|---|
| `data/cholera_deaths.tsv` | monthly recorded cholera deaths, Dacca district, former British province of Bengal | 600 rows; `time` = fractional year on the grid `1891 + 1/12, 1891 + 2/12, …, 1941.0` (the series covers Jan 1891 – Dec 1940), `cholera_deaths` = deaths recorded in that month |
| `data/census.tsv` | decadal census population, Dacca district | 6 rows (1891, 1901, 1911, 1921, 1931, 1941; ≈ 2.42M → 4.22M); interpolate to a smooth $N(t)$ — same trick the He et al. measles model uses for London |
| `king_et_al_2008_nature.pdf` | the source paper | 4 pp. + Methods; Fig. 1 = the model schematics, Fig. 2 = the profile likelihood of immunity duration for Dacca, Fig. 3 = mechanism contrast (new SIRS vs. TSIRS) |

**Provenance.** Both TSVs were extracted from the `dacca` example in
the **pomp** R package — file
[`R/dacca.R` in `kingaa/pomp`](https://github.com/kingaa/pomp/blob/master/R/dacca.R)
— on 2026-05-11 (`gh api repos/kingaa/pomp/contents/R/dacca.R`). pomp's
`dacca()` constructs a `pomp` object holding these data plus the King
et al. (2008) model at its MLE. The pomp documentation records: *"Data
are provided courtesy of Dr. Menno J. Bouma, London School of Hygiene
& Tropical Medicine."* The `time` grid, the census years/populations,
and the `cholera_deaths` integers are reproduced verbatim from
`R/dacca.R`; nothing is rescaled or re-aggregated here. (pomp ships the
model + MLE as a built-in — `?dacca`, `examples/dacca.R` — but does
**not** ship a tutorial walking through the 2008 comparison; that
analysis lives in the paper + its supplement and the pomp JSS paper,
King, Nguyen & Ionides, *J. Stat. Softw.* 69(12), 2016.)

## The King-2008 model set (from the paper)

All models are SIRS-with-Erlang-waning at heart, fitted as SDEs
(Euler–Maruyama, **Δt = 0.05 months ≈ 1.5 days** — so γ·Δt ≈ 0.09;
they were already respecting the leap condition), MIF for the fit
(J = 10⁴ particles, M = 80, cooling α = 0.95; 3×10⁴ for the final MLE
refinement), compared by log-likelihood + AICc + nested LRT. The
measured quantity is monthly deaths; observation model is
constant-CV Normal (mean = accumulated deaths over the month, SD =
mean·τ). Background mortality δ = 0.02 /yr is fixed; births feed S
(newborns susceptible); total population $H(t)$ is the time-varying
census interpolation; the force of infection λ(t) combines a human
mass-action term (∝ β_seas(t)·I^α/H) and an environmental-reservoir
term, both potentially seasonal.

| | model | structure | what it tests | paper's verdict |
|---|---|---|---|---|
| baseline | **TSIRS** | Koelle & Pascual 2004 / Koelle et al. 2005 semi-mechanistic time-series SIRS; immunity ~3–10 yr | the "previous best fit" | new SIRS beats it by **ΔAICc > 58** |
| **A** | **SIRS** (Fig 1a) | S→I at λ(t) [human mass-action + seasonally-*constant* environmental term]; I→R₁→…→R_k→S (k-stage Erlang, mean immunity 1/ε, variance 1/(kε²)); excess death rate m on I | the new mechanism — short immunity ⇒ high asymptomatic ratio | immunity **2–12 weeks** (Dacca, Fig 2); case fatality ≈ 0.004 ⇒ huge asymptomatic ratio; R₀ ≈ 1.5 |
| **B** | **two-path** (Fig 1b) | exposure → I with prob c (severe, *infectious*, dies at m) **or** → Y with prob 1−c (inapparent, *not* infectious, negligible death, short immunity 1/ρ → S). **Nests A** at c = 1, ρ = ∞ | separate the *immunological* consequence of exposure from the *infectiousness* — is it short immunity, or "silent shedders"? | short-term immunity ≈ 9.9 wk (agrees with A); case fatality ≈ 0.34 (era-consistent); severe-infection immunity ≈ 1.5 yr. Data **equivocal** on A vs B — it's the immunology (in both), not the shedding |
| **C** | **SIRS + seasonal reservoir** | A, but the environmental FOI is also seasonal (a second periodic basis). *This is the `pomp::dacca()` MLE model.* | does seasonal forcing of the aquatic path improve the fit? | better in many districts; low-R₀ / fast-immunity / high-inapparent conclusions **unchanged** |
| **D** | **environmental-phage** (Fig 1c) | A, plus phage W builds up (shed by I) and *attenuates* transmissibility — a negative feedback (the Faruque phage–vibrio seasonality hypothesis) | does phage–vibrio interaction explain cholera seasonality? | **not** favoured by the mortality data; and the predicted phage–cholera phase relation is contradicted by the one concurrent dataset (phage lags ≈ 180°, not coincident) — model comparison *and* an out-of-model check both say no |

Methods notes worth carrying into the chapter: King et al. explicitly
flag **weak non-identifiabilities among c, γ, m, ε, ρ** but note that
the *combinations* `R₀ = c·⟨β_seas⟩/(γ+δ+m)` and case fatality
`m/(γ+δ+m)` *are* well-identified — i.e. profile the right things;
they profiled over ε and γ (Fig 2 is the ε one). The B-nests-A test is
**non-standard** (at c = 1, ρ is unidentified) — Self–Liang /
Anisimova: the χ² approximation to the LRT is conservative there; they
report the P-value on 3 df (ρ, c, Y₀). Fig 2 (profile of 1/ε —
plateau → drop, 95% CI 2–12 wk) and Fig 3 (mechanism contrast: rapid
waning lets S return to high levels each season → susceptible
depletion halts the epidemic, vs. the old "spatial effects + seasonal
transmissibility drop" view) are directly reproducible chapter
figures.

## Proposed camdl model ladder

Faithful to the paper, plus one deliberate addition (M6) for the
"indistinguishable in-sample, divergent counterfactuals" lesson.
M1's `.camdl` is the core all of M2–M6 build on, so build it first.

- **M1 = A (SIRS).** β_seas(t) as a periodic spline; constant
  environmental FOI; k-stage Erlang waning (k = 3, mean 1/ε free);
  excess death rate m; δ = 0.02/yr fixed; constant-CV Normal obs on
  monthly deaths; `overdispersed()` (gamma white noise) on β for the
  environmental process noise (King's `sd_beta`). → profile 1/ε,
  reproduce **Fig 2**.
- **M2 = SIRS with 1/ε pinned at ≈ 5 yr** (the Koelle–Pascual textbook
  value). M1 vs M2 is the headline LRT — *the previously-believed
  immunity duration is decisively rejected*. (Nested, standard.)
- **M3 = B (two-path).** Adds the inapparent Y path (prob 1−c, short
  immunity 1/ρ). M1 vs M3: data equivocal; it's the immunology not the
  shedding — and a teachable moment on **boundary nesting** (the
  c = 1 / ρ-unidentified caveat).
- **M4 = C (SIRS + seasonal reservoir).** Environmental FOI also a
  periodic spline. M1 vs M4: better fit, conclusions robust. (Likely a
  "we also checked" row rather than a headline.)
- **M5 = D (environmental-phage).** Phage W shed by I, attenuates
  transmissibility. M1 vs M5: not favoured; plus the external
  phase-lag check kills it. (Also a "we also checked" row — but a good
  one: model comparison + an out-of-model check agreeing.)
- **M6 = M1 + immune boosting** — *not in the paper.* Re-exposure of an
  immune individual resets the immunity clock (R_k → R₁ at the FOI on
  R): long *intrinsic* immunity that behaves like short *effective*
  immunity in high-transmission seasons. The epistemic-laundering
  hook: M1 (genuinely short immunity) vs M6 (long + boostable) —
  likely indistinguishable in-sample, divergent in counterfactuals
  (a clean-water intervention cutting the environmental FOI, say).

**Chapter spine:** M1 → M2 → M3 → M6 as the four headline models (the
"is the textbook value right?" LRT; the "immunology vs shedding"
equivocal comparison; the "laundering" hook), with M4 and M5 as a
short "robustness — two more hypotheses, same conclusions" subsection /
two table rows. Plus the methodology section (the c/γ/m/ε/ρ sloppy
directions → profile R₀ and case fatality; the boundary-nesting LRT
caveat; reproduce Fig 2 and Fig 3) and prequential / out-of-sample
scoring (rolling windows over the ~600 monthly obs — the thing this
dataset enables that the boarding-school flu didn't).

## camdl machinery to use

- **Seasonal forcing — solved by gh#59** (`forcing {}` block;
  rebuilt camdl onto the v2 de-Boor periodic B-spline). Use
  `periodic_spline 'ratio { ... }` for β_seas(t) (and for the
  seasonal environmental FOI in M4) — this *is* the "King-2008-style
  6-coef spline" the forcing-kinds doc names. A 2-harmonic
  `fourier 'ratio` is the fewer-params alternative if the spline turns
  out to overfit. `interpolated 'count { method = "spline" }` for the
  census → $H(t)$ (matches pomp's `smooth.spline`). `cos`/`sin`/`pi`
  primitives (gh#58) are available if anything wants writing inline.
- **Backend:** `chain_binomial` — King 2008 used a continuous-state
  SDE; chain_binomial is the discrete-state analogue, fine for
  H ≈ 3M (it tracks integer compartment *counts*, not individuals).
- **Δt:** γ ≈ 20.8 /yr ⇒ infectious period ≈ 17.5 days ⇒ the leap
  condition forces `Δt ≲ 3–4 days` (King used ≈ 1.5 days). The post-fit
  Richardson dt-check (see [The Integrator Step](../integrator-step.qmd))
  will confirm. Over a 50-year integration this means ~10⁴ internal
  steps per particle-filter evaluation × ~10⁴ particles × 600 obs —
  the fits are **expensive**: they go in a `Makefile` (`sim`, `fit`,
  `profiles`, `clean`, `help` targets, the per-vignette pattern); the
  `.qmd` loads the artifacts and produces the figures/tables/asserts
  at render time, with the freshness-check opening cell the project
  CLAUDE.md prescribes.
- **`overdispersed()`** for the environmental process noise on β
  (King's `sd_beta`); a chain of R compartments for the Erlang waning;
  two transitions out of S for the severe/inapparent split (M3); a
  `W` accumulator-style compartment for the reservoir (M4) / the phage
  (M5) — all standard.

## Open / to-pin-down before / during the build

- King's λ(t): the exact functional form of the human + environmental
  FOI terms (the mixing exponent α — pomp's MLE has α = 1; whether to
  keep α estimable), and the secular `beta_trend`. Pull from the
  Supplementary Equations if needed (not in the main PDF here).
- The Erlang stage count k (pomp uses 3) — fix at 3, or compare.
- Starting values / bounds for IF2 — anchor on the `pomp::dacca()`
  MLE (in `R/dacca.R`, reproduced in the chapter prose) but converted
  to camdl's parameterization and made into reasonable starts/bounds,
  not pinned at the pomp MLE (per the project CLAUDE.md — never start
  IF2 from a "truth" / published-MLE point on real data).
