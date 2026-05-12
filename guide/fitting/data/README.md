# Boarding-school influenza data

Daily prevalence of boys "confined to bed" and "convalescent" at an
English boarding school during a 14-day H1N1 (A/USSR/77-like) outbreak
in January 1978. 763 students at risk, 512 (67 %) became ill. Three
boys were in the college infirmary on day 0 (Sunday 22 January).

## Files

| file | columns | use |
|---|---|---|
| `in_bed.tsv` | `time`, `in_bed` | headline series consumed by `fits/real.toml`; equals `prevalence(I)` in the SIR model |
| `influenza_1978_school.tsv` | `time`, `date`, `in_bed`, `convalescent` | wide canonical export — kept so chapter cells can overlay the convalescent series for context |
| `source/bmj_1978_influenza_boarding_school.pdf` | — | one-page mirror of the BMJ note that originally reported the outbreak (free / non-paywalled when this README was written) |

`time` is days since 1978-01-22 (the outbreak's first observation
day). `in_bed` is the count of boys currently in the school
infirmary; `convalescent` is the count of boys recovering but no longer
acutely ill (no longer infectious — these belong in `R` of the SIR
model, not `I`). Convalescent stays at 0 for the first four days
because the BMJ note records "boys were allowed up 36 hours after
their temperatures had returned to normal" — convalescents
materialise only on day 4.

## Provenance

Both TSVs are derived from the R package
[`outbreaks`](https://cran.r-project.org/package=outbreaks) — the
`influenza_england_1978_school` dataset (curated by Salmon, Schumacher
and Höhle as part of the package). The package is itself sourced from
the original publication:

> Anonymous (1978). "Influenza in a boarding school." *British Medical
> Journal* 1(6112): 587.

A one-page PDF copy is mirrored under `source/` for convenience and
because the BMJ note is the load-bearing source for several
modelling decisions documented in `fitting.qmd` (e.g. fixed-IC
$I(0) = 3$, the number explicitly recorded in the paper for
22 January 1978).

## Reproducing

```sh
make          # regenerate the TSVs from the R outbreaks pkg
make clean    # remove the TSVs (the PDF is committed and untouched)
```

The R script (`fetch.R`) auto-installs `outbreaks` on first run if it
isn't already present. Re-running is safe — outputs are deterministic
because the source dataset is a fixed-historical record.

## Why the PDF lives in the repo

The BMJ note is one printed page, dated 4 March 1978, and is the only
contemporaneous record of the outbreak's day-by-day counts. Its
direct quotations (population at risk, day-0 infirmary count,
vaccination context, putative Hong Kong index case) are load-bearing
in the chapter prose. The Makefile deliberately does **not** fetch
this PDF: 50-year-old hosting URLs rot, and silent breakage of
load-bearing citations is the worst kind of broken. Treat the PDF as
a research-archive artefact committed alongside the data, under
fair-use (single-page note, scholarly citation, full attribution).
