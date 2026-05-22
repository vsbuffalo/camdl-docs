# Seed-timing chapter — data, models, fits

Self-contained working directory for `guide/fitting/seed-timing.qmd` ("When did
the outbreak start?"). Everything the chapter needs lives here; nothing here is
shared with the boarding-school fitting chapters (`identifiability.qmd`,
`likelihood.qmd`, `model_comparison.qmd`), which keep their own
`../{data,models,fits}`.

## Layout

```
seed-timing/
  data/      input data + fetch scripts (provenance below); git-tracked
  models/    .camdl model definitions; git-tracked
  fits/      fit.toml configs (IF2, PMMH, survey); git-tracked
  results/   regenerated artifacts (fits, draws, surveys); GIT-IGNORED
  scripts/   helper scripts; git-tracked
  archive/   superseded fits/configs/results; GIT-IGNORED, safe to delete
  Makefile   orchestrates the full rebuild (see `make help`)
  README.md  this file
```

`results/` and `archive/` are git-ignored: results are reproducible from the
configs + data, and `archive/` is a snapshot of the prior exploration kept only
so a surprising result can be traced back. **`archive/` can be deleted at any
time.**

## Data provenance

| file | source | how to regenerate |
|---|---|---|
| `data/covid_wa_daily.tsv` | NYT `covid-19-data` `us-states.csv` (cumulative WA cases, diffed to daily), anchored t=0 = 2020-01-21 | `python data/fetch_covid_wa.py` |
| `data/covid_wa_community.tsv` | derived: `daily` with the 2020-01-21 travel-import singleton dropped | (same script) |
| `data/covid_wa_growth.tsv` | derived: `community` truncated to days 0–55 (pre-lockdown growth) | (same script) |
| `data/hagelloch_daily.tsv` | `outbreaks::measles_hagelloch_1861` R package (Pfeilsticker 1863; Oesterle 1992); prodrome-onset dates → daily incidence, t=0 = 1861-10-01 | `Rscript data/fetch_hagelloch.R` (run from `data/`) |
| `data/seed_cases.tsv` | **synthetic** — `camdl simulate models/seed_timing.camdl` at the chapter's `TRUTH` params (β=0.6, γ=0.2, λ=2, w=3, N0=5000, τ=30, ρ=0.5, k=20), seed 11 | `make sim` |

The **trustworthy anchor for the COVID-WA seed time is genomic, not the case
counts**: the outbreak-clade common ancestor (tMRCA) is dated **18 Jan–9 Feb
2020** by Bedford et al., *Cryptic transmission of SARS-CoV-2 in Washington
State*, *Science* 2020 (doi:10.1126/science.abc0523). The detected case curve is
dominated by testing ramp-up in exactly the early window that carries τ — which
is the whole point of the chapter.

## Models

| file | description | used by |
|---|---|---|
| `seed_timing.camdl` | SIR + smooth-importation seed, constant ρ (synthetic experiments) | E1/E3 synthetic |
| `seed_timing_report.camdl` | + logistic reporting ramp ρ(t); likelihood only | IF2 / survey |
| `seed_timing_report_prior.camdl` | + general-coronavirus priors | PMMH (ρ(t)) |
| `seed_timing_const_rho_prior.camdl` | constant-ρ variant + priors | model comparison |
| `seed_timing_report_od_prior.camdl` | ρ(t) + overdispersion + priors | model comparison |

Priors are **general-coronavirus** (what one would have believed at an emerging
coronavirus outbreak's start), not back-fitted to SARS-CoV-2.

## Reproduce

```
make sim        # regenerate synthetic seed_cases.tsv
python data/fetch_covid_wa.py     # COVID-WA (needs network)
Rscript data/fetch_hagelloch.R    # Hagelloch (needs R + outbreaks pkg)
make survey     # likelihood-landscape survey (interactive HTML)
make if2        # IF2 fits (fail the gate — the convergence-issues section)
make pmmh       # PMMH posteriors: informative, vague, const-rho, od
make compare    # camdl compare across the 3 PMMH models
```

Then render the chapter from `guide/fitting/`:
`quarto render seed-timing.qmd`.
