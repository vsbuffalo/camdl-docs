# Pre-COVID-19 coronavirus epi-parameter priors

Scope: parameter ranges and informative priors for SEIR-class compartmental
models of a *novel emerging coronavirus*, using only information that would
have been available to an analyst in **January–February 2020** — before
COVID-19-specific epidemiological measurements began to accumulate.

Sources are SARS-CoV-1 (2003), MERS-CoV (2012–present), and the seasonal
human coronaviruses (229E, OC43, NL63, HKU1). Where SARS- and MERS-specific
values diverge wildly (R₀, IFR), the prior should reflect that uncertainty
honestly with wide tails. Where they agree (latent/infectious periods of
~days, not hours or weeks), the prior can be moderately informative.

This document guides the SEIR re-fit of the WA seed-timing chapter. The
chapter's first round used a plain SIR model with a COVID-19-specific γ
pinned at 0.143 (1/7 d); both choices are wrong if our methodological
target is "what could a real-time analyst defensibly have inferred." See
[[two-estimation-modes]] in memory for the framing.

---

## 1. Latent period (1/σ, mean residence time in E compartment)

Pre-COVID coronavirus literature does not cleanly separate *latent* (time
to becoming infectious) from *incubation* (time to symptom onset) — most
SARS and MERS studies report the incubation period directly. For SEIR
modelling, the latent period is typically taken as **slightly shorter**
than the incubation period because pre-symptomatic transmission is common
for respiratory viruses.

### SARS-CoV-1

[Donnelly et al. 2003, *Lancet*](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(03)13410-1/fulltext) — Hong Kong, n=1425 cases:
**mean incubation 6.4 d (95% CI 5.2–7.7)**. Lognormal-fitted distribution.

[Lessler et al. 2009, *Lancet Infect Dis*](https://www.sciencedirect.com/science/article/abs/pii/S1473309909700696) — systematic
review of nine respiratory viral incubation periods:
- SARS-CoV: median ~4 d, IQR 2–7 d
- Human coronavirus (HKU1, 229E, OC43, NL63 pooled): median ~3 d

### MERS-CoV

[Cauchemez et al. 2014, *Lancet Infect Dis*](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(13)70304-9/fulltext) —
**mean incubation 5.5 d (range 2–14 d)**. Used a gamma distribution.

### Seasonal human coronaviruses (229E, OC43, NL63, HKU1)

Limited epidemiological characterisation. Lessler et al. report ~3 d median.
Generally mild, often not separately diagnosed from other URIs.

### Recommended prior for novel emerging CoV

Latent period **shorter** than incubation by ~1–2 d (presymptomatic infectiousness):

```
1/σ ~ LogNormal(mu = log(3), sigma = 0.45)
```

Mean ≈ 3.3 d, 90% CI ≈ [1.5, 7.0]. Covers SARS, MERS, and seasonal-CoV
estimates without over-committing. Equivalent prior on σ: median **σ ≈ 0.33/d**.

---

## 2. Infectious period (1/γ, mean residence time in I compartment)

### SARS-CoV-1

[Lipsitch et al. 2003, *Science*](https://www.science.org/doi/10.1126/science.1086616) — used **1/γ = 7 d**
in transmission model (Hong Kong/Singapore).

[Donnelly et al. 2003, *Lancet*](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(03)13410-1/fulltext):
**onset → recovery 26 d** (mean), **onset → death 21 d** (mean). These
include the symptomatic phase; the *infectious* phase (for transmission
modelling) is shorter because viral shedding peaks early and isolation
typically truncates the late tail.

### MERS-CoV

[Cauchemez et al. 2014, *Lancet Infect Dis*](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(13)70304-9/fulltext) —
recovery rate for community cases **γ_{I,s} = 1/5 d**, in line with SARS
literature. So **1/γ ≈ 5 d**.

### Recommended prior

```
1/γ ~ LogNormal(mu = log(5.5), sigma = 0.35)
```

Mean ≈ 5.8 d, 90% CI ≈ [3.1, 10.5]. Covers MERS (5 d), SARS (7 d), and
some uncertainty above and below. Equivalent prior on γ: median **γ ≈ 0.18/d**.

---

## 3. Generation interval (T_g = 1/σ + 1/γ for SEIR; ≈ 1/γ for SIR)

For SEIR with one E and one I compartment, the **mean** generation interval
is roughly 1/σ + 1/γ. Using the priors above: T_g median ≈ 3 + 5.5 ≈ **8.5 d**.

### SARS-CoV-1

[Lipsitch et al. 2003, *Science*](https://www.science.org/doi/10.1126/science.1086616) — generation interval ≈ 8–12 d
(Hong Kong/Singapore).

### MERS-CoV

[Cauchemez et al. 2014](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(13)70304-9/fulltext): 7–14 d.

### Distribution shape

Both the latent and infectious distributions are right-skewed; the *generation
interval* therefore has CV ≈ 0.5–0.7 (peaked, not exponential). Pure SIR
implies a CV of 1.0, which can bias the inferred R₀ when fitting to growth
data. [Park et al. 2020, *medRxiv*](https://www.medrxiv.org/content/10.1101/2020.10.17.20214262v3.full)
derive the exact relationship between growth rate r, R₀, and the shape of
the generation interval distribution; using the wrong CV biases R₀ by
factors of ~1.2–1.5 in either direction.

**Implication for the seed-timing chapter:** the SIR fit used 1/γ = 7 d and
implicitly assumed an exponential generation interval (CV = 1). The SEIR
re-fit with two-stage (E, I) chain will give CV ≈ 0.7, closer to the SARS
truth, and R₀ inferences will differ by ~10–20 %.

---

## 4. R₀ (basic reproduction number)

### SARS-CoV-1

- [Riley et al. 2003, *Science*](https://www.science.org/doi/10.1126/science.1086478): Hong Kong, R₀ ≈ **2.7**.
- [Lipsitch et al. 2003, *Science*](https://www.science.org/doi/10.1126/science.1086616): Hong Kong/Singapore, R₀ ≈ **3** without
  control; reported uncertainty range **1.1–7.7** across settings.
- [Anderson et al. 2004, *Phil Trans R Soc B*](https://royalsocietypublishing.org/doi/10.1098/rstb.2004.1490): R₀ ≈ 2–4 across most analyses.

### MERS-CoV

- [Cauchemez et al. 2014, *Lancet Infect Dis*](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(13)70304-9/fulltext):
  R₀ ≈ **0.63** in observed data; **0.8–1.3** in counterfactual without
  infection control. Indicates limited sustained human-to-human
  transmission under routine surveillance.
- [Breban et al. 2013, *Lancet*](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(13)61492-0/fulltext): R₀ < 1
  consistent with sporadic clusters rather than self-sustaining epidemic.

### Seasonal human coronaviruses

No clean R₀ estimates in the human population (endemic, age-structured,
many subclinical). Pandemic influenza literature suggests R₀ ≈ 1.5–2.5
as a comparison class for "transmissible respiratory virus."

### Recommended prior

For a *novel emerging coronavirus* of unknown severity in Feb 2020, the
SARS analogue is the relevant reference (highly transmissible, severe).
MERS is the cautionary case (low R₀ but high CFR; might suggest the new
agent is closer to that profile).

```
R0 ~ LogNormal(mu = log(2.5), sigma = 0.5)
```

Mean ≈ 2.8, 90% CI ≈ [1.1, 5.7]. Covers SARS (≈2.7), the lower MERS
range (~1), and pandemic-influenza territory. Does not include the
implausibly-high tail of Lipsitch's 7.7 estimate but accepts it as
"possible but unlikely."

Equivalent prior on β (with γ as defined above, median 0.18): median
**β ≈ 0.50/d**.

---

## 5. Infection fatality ratio (IFR)

**This is the parameter with the widest uncertainty in pre-COVID
coronavirus literature.** Pre-2020, the three reference points span
three orders of magnitude:

| pathogen | IFR | source |
|---|---|---|
| SARS-CoV-1 | ~**10 %** | [WHO 2003 CFR](https://www.who.int/csr/sars/country/table2004_04_21/en/) (CFR ≈ IFR for severe-only ascertainment) |
| MERS-CoV | ~**34 %** | [WHO MERS dashboard](https://www.who.int/health-topics/middle-east-respiratory-syndrome-coronavirus-mers); biased by hospital-detected cases |
| Seasonal CoV | ~**0.01–0.1 %** | inferred from age-stratified excess mortality; mostly elderly |
| 1918 influenza | ~**1–3 %** | classic catastrophe reference |
| Seasonal influenza | ~**0.05–0.1 %** | typical year |

For a novel emerging coronavirus in Feb 2020 (severity unknown, first cases
mostly hospital-ascertained → upward biased), the SARS-like prior was the
common default for emergency planning, but the wide range matters.

### Recommended prior

Three tiers from most to least diffuse:

```
IFR ~ LogUniform over [0.0005, 0.5]    # diffuse — 3 orders of magnitude
IFR ~ LogNormal(log 0.02, 1.0)         # weakly informative — anchored at 2%
IFR ~ Beta(2, 100)                     # informative — SARS + pandemic flu class
```

The **weakly informative** middle option `LogNormal(log 0.02, 1.0)` is the
right default for the seed-timing chapter:

- Median 0.02 (2%), 90% CI ≈ [0.004, 0.10], P(IFR > 0.5) ≈ 0.0001.
- The 2% anchor is the log-mean of seasonal CoV (~0.1%) and SARS-CoV-1 (~10%) — the natural midpoint of the "novel respiratory CoV with severity unknown" reference class.
- σ = 1.0 keeps the SARS extreme (10%) within the 90% CI and the seasonal extreme (0.1%) within the 95% CI. The MERS extreme (34%) is in the prior tail at P ≈ 0.001 — admitted but unlikely.
- Does not commit the analysis to a specific reference pathogen, but does encode that 0.0005 (1-in-2000) and 0.5 (1-in-2) are physically implausible for a virus producing the observed clinical syndrome.

Use the diffuse `LogUniform` only when *deliberately* wanting to refuse to
constrain. Use the informative `Beta(2, 100)` only when ready to defend the
SARS+pandemic-flu reference class on the record.

**Toolchain note (camdl).** In camdl, parameters declared as `: probability`
use a Logit transform on the unconstrained scale, and prior-transform
compatibility (validated in `fit::runner::validate_prior_transform_compat`)
requires the prior family to match: `Beta` and `Uniform` on Logit, `LogNormal`
on Log. So the weakly-informative `LogNormal(log 0.02, 1.0)` middle tier
above is **not directly representable** for a `probability`-kind IFR
parameter in a `.camdl` model file. Two options:

- Use `Beta(2, 100)` (the informative tier) — same mean as the LogNormal
  middle tier, slightly tighter tails, fully native to the Logit transform.
- Declare `ifr : positive in [0, 1]` (Log transform) with the LogNormal
  prior — semantically weird (probability-as-positive) but lets the prior
  shape come through. Not recommended for chapter use.

Practical recommendation: take the toolchain at its word and use
`Beta(2, 100)`. The numerical difference vs LogNormal(log 0.02, 1.0) in
posterior IFR for a fit with informative deaths data is negligible.

---

## 5b. Onset-to-death lag (1/γ_d, mean residence time in pre-death compartment P)

In an SEIR+death model with a P (pre-death) compartment, the death pathway is

```
I --(1-ifr)·γ·I--> R                            (survives)
I -- ifr·γ·I    --> P                            (fated to die: enters pre-death)
P --   γ_d·P    --> D                            (dies after ~1/γ_d days)
```

The **total onset-to-death lag** = 1/γ (mean time in I) + 1/γ_d (mean
time in P). Pre-COVID, this lag was measured for SARS-CoV-1 and MERS-CoV:

| pathogen | symptom-onset-to-death (median or mean) | source |
|---|---|---|
| SARS-CoV-1 (Hong Kong 2003) | **median 36 days** (IQR 26–53), mean truncated by follow-up | [Donnelly et al. 2003, *Lancet*](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(03)13410-1/fulltext) Table 1 |
| SARS-CoV-1 (HK 2003, recomputed) | **mean 23.7 days** (gamma fit) | [Ghani et al. 2005, *AJE*](https://academic.oup.com/aje/article/162/5/479/138989) |
| MERS-CoV | **median ~14 days** (range 5–36) | [Assiri et al. 2013, *NEJM*](https://www.nejm.org/doi/full/10.1056/NEJMoa1306742) |
| 1918 influenza (severe pneumonia) | **~7–14 days** | classical pathology series |
| Severe community-acquired pneumonia | **~10–20 days** | hospital-based studies |

With 1/γ ≈ 6 d (typical infectious period for an emerging respiratory CoV
from §2), the implied 1/γ_d that matches a 20-day onset-to-death lag is
**~14 days** (γ_d ≈ 0.07/d). For a SARS-like 24-day lag, 1/γ_d ≈ 18 days
(γ_d ≈ 0.055/d).

### Recommended prior

```
gamma_d ~ LogNormal(log 0.077, 0.4)    # median γ_d = 0.077/d → 1/γ_d ≈ 13 d
```

- Median 0.077/d corresponds to mean P-stage duration ≈ 13 days, giving total
  symptom-onset-to-death ≈ 6 + 13 = 19 days — within the 14–24 day pre-COVID
  range above.
- σ = 0.4 (on log scale): 90% CI roughly [0.040, 0.148] /d ⇒ 1/γ_d ∈ [6.8, 25] d
  ⇒ total lag ∈ [13, 31] d. Covers MERS-fast through SARS-slow.
- Does **not** commit to a specific clinical course; an analyst willing to
  pin γ_d at the literature mean (1/γ_d = 13 d, γ_d = 0.077) would do so
  via `[fixed]` instead.

### Observation overdispersion for deaths (k_d)

For the negative-binomial observation `deaths ~ NegBin(mean = projected, r = k_d)`:

```
k_d ~ LogNormal(log 20, 1.0)    # median r = 20 → variance ≈ mean × (1 + mean/20)
```

Same shape as the cases-stream prior (§6). Pre-COVID, daily-death counts
in respiratory outbreaks were typically more overdispersed than cases due
to clustering of severe outcomes in nursing homes and hospitals; the
LogNormal(log 20, 1.0) median r=20 allows for substantial overdispersion
(variance ≈ 6× mean at mean=100) while admitting near-Poisson (r ≫ mean) in
the tail.

### Ascertainment of deaths (rho_d)

```
rho_d = 1    # fixed
```

We assume all deaths are observed. This is a **simplifying assumption**:
in reality, the earliest WA COVID-19 deaths were probably attributed to
"pneumonia of unknown cause" until retrospective re-classification. But for
a teaching-chapter joint cases+deaths fit, the IFR parameter absorbs both
true IFR and any death under-reporting in the early window. Sensitivity to
this assumption should be flagged but not estimated jointly without a
second independent informant on death ascertainment.

---

## 6. Reporting / ascertainment fraction (ρ)

There is **no useful pre-COVID prior** on case-ascertainment for a novel
emerging coronavirus. Ascertainment depends entirely on the local
surveillance system, testing capacity, and case-definition stringency at
the time of detection. SARS in Hong Kong had near-complete ascertainment
once the outbreak was recognised (~90 %); MERS in Saudi Arabia varied
from ~30 % (active surveillance) to <5 % (community circulation, never
detected). Seasonal CoV is rarely measured at all.

### Recommended priors

For *late-detection* surveillance (a few weeks of testing buildup):

```
rho_max ~ Beta(2, 8)    # mean 0.20, 90% CI [0.05, 0.45]
```

For the *ramp-onset time* t_rep — the date when the local lab/surveillance
system reaches half of its eventual capacity — anchor it from independent
information about the surveillance timeline:

```
t_rep ~ Normal(t_external, sigma = 5 d)
```

where t_external is the documented date of test availability in the region
(e.g., for WA in early 2020: CDC assay Feb 26, UW Virology Feb 29 → centre
prior on t_external ≈ Feb 28).

---

## 7. Seed size / introduction count

Genomic monophyly evidence (when available) directly constrains the
**number of founding lineages** (n_seed). In the WA seed-timing chapter,
[Bedford et al. 2020, *Science*](https://www.science.org/doi/10.1126/science.abc0523)
identified one dominant clade (~384 sampled viruses) plus one small
secondary clade (~39 viruses) → effective n_seed in the small-integer
range (1–10) for the dominant outbreak.

Without genomic data, n_seed is effectively unidentified from case-count
data (this is the central finding of the seed-timing chapter). Use a
wide prior:

```
n_seed ~ LogUniform over [1, 1000]
```

---

## 8. Summary table — recommended SEIR priors for novel emerging CoV (Feb 2020 vintage)

| param | meaning | median | 90% CI | distribution |
|---|---|---|---|---|
| 1/σ | latent period (d) | 3.3 | [1.5, 7] | LogNormal(log 3, 0.45) |
| 1/γ | infectious period (d) | 5.8 | [3.1, 10.5] | LogNormal(log 5.5, 0.35) |
| R₀ | basic reproduction | 2.8 | [1.1, 5.7] | LogNormal(log 2.5, 0.5) |
| IFR | infection fatality ratio | 0.02 | [0.003, 0.05] | Beta(2, 100) — SARS-anchored |
| | (or 0.05 with [0.001, 0.5] if class-uncommitted) | | | LogUniform |
| ρ_max | ascertainment fraction | 0.20 | [0.05, 0.45] | Beta(2, 8) |
| t_rep | surveillance ramp midpoint | t_doc | ± 5 d | Normal(t_doc, 5) |
| n_seed | founder count | (genomic) | (wide) | LogUniform(1, 1000) |
| w_rep | surveillance ramp width | 5 | [2, 12] | LogNormal(log 5, 0.4) |

**Derived sanity check:** with these priors, the generation interval has
median ≈ 8.5 d, growth-rate-implied R₀ has median ≈ 2.8, and the case-to-
death ratio at steady state is ρ_max × incidence / (IFR × incidence_lagged).

---

## 9. How these compare to what the WA seed-timing chapter v1 used

| param | v1 (COVID-specific, with hindsight) | v2 priors above (Feb-2020 vintage) | shift |
|---|---|---|---|
| 1/γ | 7.0 d (pinned) | 5.8 d median (estimated) | shorter, less certain |
| 1/σ | — (SIR, no latent) | 3.3 d median (estimated) | NEW dimension |
| R₀ | β/γ ≈ 3.1 (estimated) | 2.8 median (estimated) | similar central, wider prior |
| ρ_max | 0.10 (pinned at lit) | 0.20 median (estimated) | higher central, much wider |
| IFR | 0.007 (pinned) | 0.02 median or LogUniform | 3× higher central, much wider |
| Generation interval CV | 1.0 (SIR) | 0.7 (SEIR) | sharper, more realistic |

The v2 priors are *broader* and *more honest* about what an analyst in
Feb 2020 could reasonably have assumed. They will produce wider posterior
intervals for τ and R₀ but more defensible point estimates.

**Predicted directional shift of τ̂ from v1 to v2:** SEIR adds the latent
period as a forward delay (~3 d), so τ̂ shifts ~3 d **earlier** than the
SIR estimate, all else equal. This brings the case-data τ̂ ~3 d closer to
the Bedford genomic centroid (Feb 12) without changing the identifiability
ridge structure. Sensitivity to ρ_max, t_rep, n_seed — unchanged in
direction, possibly narrower in magnitude because the more flexible priors
let the data speak more.

---

## 10. Open questions for the next round

1. **Should we use the diffuse IFR prior or the SARS-anchored one?** Diffuse
   is more honest; SARS-anchored is what a real analyst would have used.
   For the chapter we should probably show both and discuss the difference.
2. **Should γ_d (the death-lag time) be a separate parameter or fixed?**
   With the SEIR + death-pathway model from earlier work, this is its own
   parameter and is *also* poorly known a priori. Probably pin from
   literature (~13 d post-infectious to death) but acknowledge.
3. **w (seed pulse width) vs discrete events?** For SEIR the issue is moot
   in principle (the E compartment provides a natural smoothing) — a
   "discrete pulse into E" then spreads naturally over 3 d before
   becoming infectious. So the SIR-era debate about narrow-Gaussian vs
   discrete might not arise the same way.

---

## Sources

- [Assiri A et al. 2013. Epidemiological, demographic, and clinical characteristics of 47 cases of Middle East respiratory syndrome coronavirus disease from Saudi Arabia. *NEJM* 369:407–416.](https://www.nejm.org/doi/full/10.1056/NEJMoa1306742)
- [Donnelly CA et al. 2003. Epidemiological determinants of spread of causal agent of severe acute respiratory syndrome in Hong Kong. *Lancet* 361:1761–1766.](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(03)13410-1/fulltext)
- [Ghani AC et al. 2005. Methods for estimating the case fatality ratio for a novel, emerging infectious disease. *AJE* 162:479–486.](https://academic.oup.com/aje/article/162/5/479/138989)
- [Lipsitch M et al. 2003. Transmission dynamics and control of severe acute respiratory syndrome. *Science* 300:1966–1970.](https://www.science.org/doi/10.1126/science.1086616)
- [Riley S et al. 2003. Transmission dynamics of the etiological agent of SARS in Hong Kong: impact of public health interventions. *Science* 300:1961–1966.](https://www.science.org/doi/10.1126/science.1086478)
- [Anderson RM et al. 2004. Epidemiology, transmission dynamics and control of SARS: the 2002–2003 epidemic. *Phil Trans R Soc B* 359:1091–1105.](https://royalsocietypublishing.org/doi/10.1098/rstb.2004.1490)
- [Lessler J et al. 2009. Incubation periods of acute respiratory viral infections: a systematic review. *Lancet Infect Dis* 9:291–300.](https://www.sciencedirect.com/science/article/abs/pii/S1473309909700696)
- [Cauchemez S et al. 2014. Middle East respiratory syndrome coronavirus: quantification of the extent of the epidemic, surveillance biases, and transmissibility. *Lancet Infect Dis* 14:50–56.](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(13)70304-9/fulltext)
- [Breban R, Riou J, Fontanet A. 2013. Interhuman transmissibility of Middle East respiratory syndrome coronavirus: estimation of pandemic risk. *Lancet* 382:694–699.](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(13)61492-0/fulltext)
- [Park SW et al. 2020. Exact solution of infection dynamics with gamma distribution of generation intervals. *medRxiv* 2020.10.17.20214262.](https://www.medrxiv.org/content/10.1101/2020.10.17.20214262v3.full)
- [Bedford T et al. 2020. Cryptic transmission of SARS-CoV-2 in Washington state. *Science* 370:571–575.](https://www.science.org/doi/10.1126/science.abc0523)
