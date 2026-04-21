# Boarding-school SIR chapter — reorganization plan

Captured from the planning discussion, 2026-04-19.

## Module layout — `camdl_diag/` (sibling to `report.qmd`, self-contained)

All shared plotting / CLI / Rhat helpers move into one subdirectory at the
same level as `report.qmd`. Nothing outside this worktree is touched. The
module is designed to later extract cleanly into a standalone
`camdl-diag` package.

```
camdl_diag/
  __init__.py       # re-exports the public surface
  styles.py         # (move existing) GREY, CMAP_DLL, CMAP_NUIS,
                    #   CI_LEVELS_1D/2D, marker dicts, apply_rc()
  cli.py            # NEW — run_cli, _ansi_to_html
                    #   (lifted verbatim from guide/fitting.qmd:48)
  rhat.py           # NEW — parse_rhat + render_rhat_table
                    #   (HTML table, traffic-light colored by threshold)
  traces.py         # (rename plot_traces.py) load_chains,
                    #   plot_traces, plot_rank_overlay, plot_pair_grid
  profile_plot.py   # (move existing) plot_profile_2d + panel helpers
  diagnostics.py    # NEW — render_fit_diagnostics(fit_dir, ...)
                    #   → traces + rank + pair + Rhat table in one call
```

Report-side imports collapse to:

```python
from camdl_diag import (run_cli, render_fit_diagnostics,
                         plot_profile_2d, plot_pair_grid, apply_rc)
```

## Ordering decision — obs model before free-IC

**Rule:** fix the observation likelihood first, then push on the
dynamics. Reasoning:

- Modeler's workflow matches: pick the obs model before inference on
  parameters. Otherwise free-IC has to be run twice (Poisson then
  NegBin), and the geometric ridge (β–γ vs I₀–R₀) is the same under
  either likelihood, so the first run adds narrative cost without
  payoff.
- NegBin's dispersion `k` is its own sloppy axis that only makes sense
  to explore once NegBin is adopted. Lumping k-sloppiness and
  IC-identifiability into one obs model reads cleanly; interleaving
  doesn't.
- Small cost: the current "Poisson free-IC fails" section is good
  motivation for profile likelihoods. Under the new order the
  motivation shifts to NegBin — arguably cleaner, because NegBin
  already absorbs noise so residual non-convergence points directly at
  identifiability rather than misspecification.

The only case for free-IC-on-Poisson-first is pedagogical simplicity
(fewer moving parts). Not worth it — the reader sees the ridge under
NegBin just as clearly.

## Outline — what to keep, in the new order

### Page 1 — Building and calibrating the fit

1. The data
2. The model (`camdl` source)
3. Synthetic pipeline validation — SBC / recovery (Rhat floor-check
   for the whole pipeline; concise, fold details)
4. Baseline fit: Poisson, fixed IC — `camdl fit run fit.toml`
   - Rhat table; scout + refine traces (collapsible); two-panel PPC
5. Choosing the observation model: Poisson → NegBin
   - Side-by-side Poisson vs NegBin fits (fixed IC)
   - Visual + Rhat + raw ll comparison
   - Lands on: NegBin fits better; subsequent analyses use NegBin
6. Backend-consistency callout (keep as-is)

### Page 2 — Identifiability under NegBin

1. Freeing initial conditions — 4-param NegBin fit, fails to converge
   - Scout + refine trace plots, Rhat table, Δll across chains
2. 1D profile: is I(0) really unidentified, or is it joint?
3. 2D profile likelihoods — the identifiability ridge
   - (I₀, R₀) profile with (β, γ) optimized per cell
   - (β, γ) profile with (I₀, R₀) optimized per cell
   - Chain trajectories on the profile surfaces
   - Scout pair plot
4. The `k` sloppy axis — NegBin-specific, predictions along it

### Page 3 (short) — Model comparison, formally

Deserves its own page: *introducing* AIC / ΔAIC / BIC / LOO is a
reference section other chapters will link to.

1. Why "NegBin fits better" needs a proper measure
   (over-parameterization concern)
2. AIC / BIC — definitions, assumptions, caveats for SSMs
   (IF2's particle-filter MLE isn't regular MLE in the textbook
   sense; worth a sentence)
3. Leave-one-out log-likelihood — time-series analogue
   (leave-one-day-out or leave-last-k for forecast validation) —
   independent confirmation that NegBin's advantage isn't just
   extra-parameter flexibility
4. Final comparison table: Poisson vs NegBin on ΔAIC, ΔBIC, LOO-LL
5. Conclusion + reading guide

## Phases (execution order)

**Phase 0 — create `camdl_diag/`**
- Move `styles.py`, `profile_plot.py`, `plot_traces.py → traces.py`
- Add `cli.py` (lift from guide/fitting.qmd), `rhat.py`, `diagnostics.py`
- Update driver-script imports. Nothing visible changes.

**Phase 1 — inline CLI in report.qmd**
- Replace ```bash fences with `run_cli(...)` calls, one fit section at
  a time. Start with baseline `fit.toml`.

**Phase 2 — diagnostic coverage parity**
- Every major fit section gets: scout + refine traces (collapsible),
  Rhat table. Missing currently on NegBin and long-NegBin-refine.

**Phase 3 — driver consolidation**
- `render_fit_diagnostics(fit_dir, *, params, labels)` collapses each
  fit section to ~5 lines in the qmd.

**Phase 4 — split into pages + prune**
- Break `report.qmd` into the three pages above; delete duplicate /
  legacy scripts and unreferenced figures.
