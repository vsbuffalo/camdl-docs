# camdl-book — project instructions for Claude

## Fitting diagnostics: always check convergence first

When running any `camdl fit` pipeline (scout, refine, validate), **always report
and check Rhat for every estimated parameter before drawing any conclusions
from the MLE**. This is non-negotiable.

- Rhat ≥ 1.1 on any parameter → the fit is **not converged**. Do not interpret
  the MLE, do not run PPC on it, do not include it in the report. Instead:
  rerun with more iterations, tighter cooling, more chains, wider `rw_sd`, or
  narrower bounds — whichever diagnostic indicates.
- Rhat < 1.05 on all parameters → fit is converged. Safe to interpret.
- Rhat 1.05–1.10 → marginal; report but flag; preferably tighten and rerun.
- Also check: divergent-chain count, `best_loglik` vs `initial_loglik`
  progression, and whether any parameter hit a bound at convergence.

**Include Rhat as a column in every MLE parameter table in the report.** Do
not hide it. Readers should see the convergence status alongside the point
estimates.

**Root-cause reruns, don't paper over.** If Rhat is bad, find out why — ridge?
too-tight rw_sd? pinned at bound? multimodal? — rather than just cranking
iteration counts blindly.

**Never use `--allow-nonconverged-scout`.** The upstream scout-Rhat gate
exists to stop refine from laundering multi-modal scout output into a
chain-12-endpoint "MLE" that isn't a converged fit. If scout's tail-Rhat
fails, the remedy is one of: narrow bounds toward scout's best basin,
increase scout chains or iterations, set informed `start` values, or
declare the parameters unidentified from this data and stop. Do not
bypass the gate to produce a report-able number — you will get numbers
that look fine per-parameter (low refine Rhat) while chain log-lik
spreads 50+ nats across chains. This burns reviewer trust when the
narrative turns out to rest on laundered output. Same rule applies to
any upstream convergence gate: if a gate refuses, fix the underlying
problem, do not pass a flag to ignore it.

## Self-contained chapters

Each `.qmd` chapter should be **self-contained** — all imports, helpers, and
path setup live in the file itself (or import from `_lib/`). No implicit
dependencies on other chapters having been rendered first. If a chapter
genuinely can't be self-contained (shared state, sequential render order),
flag it to Vince rather than silently coupling files.

## Local preview

- Tailscale hostname: `thuja`
- Serve locally with `python3 -m http.server 8787 --directory _site`
- Preview URL: `http://thuja:8787/`

## General

- Upstream first — see the user's global memory.
- Never write aspirational CLI syntax — see the user's global memory.
