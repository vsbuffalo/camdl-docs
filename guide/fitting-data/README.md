# boarding-school-sir-clean

A minimal, from-scratch SIR fit of the 1978 English boarding-school influenza
outbreak. Teaching version — shows the IF2 fit, the two-panel diagnostic, and
synthetic calibration at the MLE. Nothing else.

## Why this exists

The parent `.fit-scratch/` tree ran into a [backend-default mismatch
bug](/tmp/camdl-incidents/2026-04-19-backend-default-mismatch.md) where
`camdl simulate` and `camdl fit` default to different backends. Every PPC
figure in that tree was contaminated. This directory starts fresh with
explicit `--backend chain_binomial --dt 1.0` on every simulation call.

## Run the pipeline

```bash
# 1. Compile the model (once)
camdl compile boarding_school_sir.camdl > boarding_school_sir.ir.json

# 2. Fit
camdl fit run fit.toml --seed 42

# 3. Two-panel diagnostic at the MLE
uv run --with polars --with matplotlib python two_panel.py

# 4. Synthetic calibration at the MLE (~10 min wall)
uv run --with polars --with numpy --with matplotlib python sbc.py

# 5. Render report
uv run --with polars --with jupyter --with matplotlib --with numpy \
    quarto render report.qmd
```

## Layout

```
boarding_school_sir.camdl   # SIR + Poisson-obs model source
boarding_school_sir.ir.json # compiled IR (camdl compile output)
fit.toml                    # IF2 fit config (chain-binomial, dt=1, explicit)
data/in_bed.tsv             # 14 daily observations
two_panel.py                # Panel A (simulate) + Panel B (pfilter --save-paths)
sbc.py                      # 30-dataset SBC at the MLE
report.qmd                  # Quarto report wiring it all together
output/                     # camdl outputs (content-addressable)
figures/                    # rendered PNGs
```

## Non-goals

- **No NegBin observation.** Poisson only. The matron counts are direct
  observations; Poisson is the natural model and avoids an unnecessary
  dispersion parameter.
- **No tvbeta, Erlang, Γ-noise on β, IC-free, or SEIR/SIBR variants.**
  Those live in the extended analysis if/when we bring them back. This
  directory is the teaching baseline.
- **No multi-truth coverage experiment.** SBC at one truth (the MLE) is
  a calibration check, not a coverage test — less demanding, enough for
  a teaching chapter.
