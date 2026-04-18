# SBC Bias Discrepancy: +59% β Bias vs ~0% Bias

**Date:** 2026-04-18
**Status:** root cause identified, upstream fix pending
**Severity:** high — affects any SBC validation workflow

## Summary

An SBC validation of IF2 on a boarding school SIR model produced
dramatically different results depending on how synthetic data was
generated. Using `camdl simulate --obs-only` (the bash-scripted
pipeline), β was biased +59% across 20 synthetic datasets. Using the
new `[synthetic]` block in `camdl fit run`, β was biased -0.4% across
50 datasets. The root cause is that the two code paths use different
RNG seeding, producing completely different stochastic realizations at
the same nominal seed. The +59% bias was real — reproducible across
binary versions — but unrepresentative of the true MLE bias, which is
near zero.

## Timeline

### Phase 1: discovering the bias (2026-04-17)

We built a boarding school SIR model (`boarding_school_sir.camdl`) with
NegBin observations (`neg_binomial(mean = projected, r = k)`), no
reporting fraction (rho fixed at 1). Generated 30 synthetic datasets
with a bash loop:

```bash
for dseed in $(seq 1 30); do
    camdl simulate boarding_school_sir.camdl \
        --params data/sbc/true_params.toml \
        --scenario baseline --seed $dseed \
        --obs-only data/sbc/syn_${dseed}.tsv
done
```

Fit each with IF2 (24 chains × 1000 particles × 100 iterations scout,
8 chains × 2000 particles × 120 iterations refine). True parameters:
β=2.1961, γ=0.5782, k=50.35.

**Result:** β mean = 3.302 (+50% bias), γ mean = 0.659 (+14%), k mean
= 22.6 (-55%). R₀ = β/γ biased +39%.

This was internally consistent: the β-k correlation was strong
(high β compensated by low k), the bias scaled with observation
frequency (collapsed by K=28), and the pattern matched theoretical
expectations for finite-sample MLE bias with a nuisance scale parameter.

### Phase 2: the observation model tradeoff (2026-04-17–18)

We tested Poisson and Binomial observation models with well-specified
SBC (each model generates and fits from its own observation likelihood).
Each used a similar bash loop with `--obs-only` data generation.

| Model    | β bias  | γ bias  | R₀ bias |
|----------|---------|---------|---------|
| NegBin   | +50%    | +14%    | +39%    |
| Poisson  | +3%     | +41%    | -26%    |
| Binomial | +7%     | +42%    | -24%    |

This produced the "bias-sink principle": finite-sample bias
redistributes across parameters depending on the observation model.
NegBin's k absorbs noise, freeing γ but letting β drift. Poisson/
Binomial force γ to absorb noise instead. The narrative was
internally coherent and theoretically motivated.

### Phase 3: the upstream `[synthetic]` feature (2026-04-18)

The upstream repo shipped a `[synthetic]` block for `camdl fit run`
(commits 251c2a6 through 331a98c). We rewrote the SBC using this
feature:

```toml
[synthetic]
true_params = "data/true_params_negbin.toml"
datasets = 30
sim_seeds = "1:30"
```

One config file, one command, auto-generated `summary.tsv` and
`coverage.tsv`.

**Result:** β mean = 2.177 (-0.9% bias), γ mean = 0.592 (+2.4%),
k mean = 51.6 (+2.4%). All parameters essentially unbiased. All
three observation models showed <3% bias.

This contradicted the earlier finding completely.

### Phase 4: investigation (2026-04-18)

**Hypothesis 1: binary/code difference.** Tested by building the old
binary (commit 76cb7a2) and running the old SBC loop. Result: +59%
bias, identical to the original finding. Then ran the same data
through the new binary: +59% bias, byte-for-byte identical MLEs.
**Rejected** — the binary doesn't matter.

**Hypothesis 2: cache contamination.** The old bash SBC used `/tmp/`
config files that could collide. Tested by deleting all cached results
(`rm -rf results/fits/fit_sbc_*`) and re-running with the new binary.
Result: +59% bias again. **Rejected** — not a cache issue.

**Hypothesis 3: scout vs refine.** The old pipeline might have been
reading scout MLEs instead of refine. Checked: the bash script
explicitly read from `results/fits/.../refine/mle_params.toml`. Also
ran a clean scout-only SBC with the `[synthetic]` feature: scout β
bias was +3.7%, refine was -0.9%. Scout heating bias is real but
modest — nowhere near +59%. **Partially confirmed** (scout is
slightly biased) but doesn't explain the magnitude.

**Hypothesis 4: model priors.** The old model file had `~ log_normal`
and `~ half_normal` prior declarations on parameters; the new model
didn't. Tested by running `[synthetic]` SBC with the old model file
(priors included). Result: identical coverage table — priors make
zero difference to IF2. **Rejected.**

**Hypothesis 5: different data.** Compared synthetic dataset #1 from
both generation methods:

```
# camdl simulate --obs-only (seed=1):
time  in_bed
0     10
1     19
2     64
3     172
4     354    ← peaks day 4

# [synthetic] feature (sim_seed=1):
time  in_bed
0     10
1     5
2     27
3     60
4     161
5     208
6     425    ← peaks day 6
```

**Same seed, completely different realizations.** The RNG state
diverges immediately after the first time step. This is the root
cause.

### Phase 5: confirmation (2026-04-18)

Fit the OLD data (from `--obs-only`) with the NEW binary:
β bias = +59%. Identical to old binary. **Confirmed: the bias is in
the data, not the software.**

Ran `[synthetic]` SBC with 50 datasets: β bias = -0.4%. **Confirmed:
the `[synthetic]` data produces unbiased results at scale.**

## Root cause

The `camdl simulate --obs-only` code path and the `[synthetic]` data
generator in `camdl fit run` use different RNG seeding strategies.
At the same nominal seed, they produce completely different stochastic
realizations. The `--obs-only` batch (seeds 1:20) happened to produce
realizations that are systematically biased when fit with IF2. The
`[synthetic]` batch (seeds 1:30, 1:50) produces representative
realizations.

Specifically:
- `--obs-only` data peaks earlier and higher (day 4, 354 boys)
- `[synthetic]` data peaks later and with more variation (day 6, 425)
- The early-peaking `--obs-only` data is consistent with higher β,
  which is what IF2 finds — the MLE is correct *for that data*
- But these realizations are not representative of the generative
  distribution at the true parameters

The question is whether this is (a) a seeding bug where the two paths
should agree, or (b) a legitimate difference in how stochasticity is
sampled that happens to produce different representative batches at
small N=20-30.

## What we got wrong

1. **We accepted the +59% bias as a fundamental statistical finding**
   without checking the data generation step. The internal consistency
   (β-k correlation, scaling with K, bias-sink tradeoff) made it look
   real. It was real — but only for that particular batch of data.

2. **We built a narrative ("bias-sink principle") on 20 datasets.**
   With N=20 stochastic realizations at K=14, the sample mean of the
   MLEs has a standard error of roughly σ/√20. If σ is large (which
   it is — individual MLEs ranged from β=2.1 to β=4.6), the sample
   mean can be far from the population mean.

3. **We assumed `--obs-only` and `[synthetic]` would produce the same
   data at the same seed.** They don't. This is either a bug or a
   documentation gap in camdl.

## What we got right

1. **The SBC methodology itself was correct.** Both batches of
   synthetic data were generated from the true parameters and fit
   with well-specified models. The fits converged. The bias we
   measured was real for the data we had.

2. **The scout vs refine investigation** correctly identified that
   scout has a +3.7% heating bias on β and +26.5% on k. This is a
   real methodological finding independent of the data issue.

3. **The observation model comparison on real data** was unaffected.
   The real-data MLEs (NegBin: β=2.06, Poisson: β=1.94, Binomial:
   β=1.89) are consistent across binary versions and data generation
   methods.

## Remaining questions

1. **Is the RNG difference a bug?** The `[synthetic]` feature and
   `--obs-only` should produce identical data at the same seed for
   reproducibility. If they don't, the seeding needs to be
   reconciled. Filed in agent-channel.

2. **Which batch is more representative?** 50 datasets from
   `[synthetic]` showing ~0% bias is more convincing than 20 from
   `--obs-only` showing +59%. But we should verify with even more
   datasets (100+) and ideally cross-validate with POMP.

3. **Is K=14 MLE truly unbiased for this SIR?** The 50-dataset
   result says yes (<1% β bias). This is theoretically surprising
   — classical O(1/K) bias should be visible at K=14. Either the
   Fisher information is higher than expected, or the chain-binomial
   backend has favorable properties for bias.

4. **Does the bias-sink principle hold at all?** It might — just not
   at the +59%/+41% magnitude we reported. With proper data
   generation, the NegBin vs Poisson difference might be a few
   percent, not 50 percentage points. Need to re-run the three-model
   comparison with the `[synthetic]` feature for all three models.

## Action items

### Upstream (camdl)
- [ ] Investigate RNG seeding difference between `--obs-only` and
      `[synthetic]` (filed in agent-channel 2026-04-18)
- [ ] Ensure both paths produce identical data at the same seed
- [ ] Document the seeding strategy

### Book (camdl-book)
- [ ] Rewrite fitting chapter using `[synthetic]` feature throughout
- [ ] Re-run all SBC experiments with consistent data generation
- [ ] Re-evaluate the bias-sink narrative with clean data
- [ ] Consider adding this incident as a methodological note:
      "we initially observed +59% bias from 20 datasets; with 50
      datasets the bias was <1%; here's what happened"
- [ ] POMP cross-validation on 5-10 synthetic datasets for
      independent confirmation

### Methodology
- [ ] Determine minimum number of SBC datasets needed for reliable
      bias estimation at K=14 (N=20 was insufficient)
- [ ] Document the scout vs refine heating bias as a practical
      guideline: always run refine, scout alone overestimates k
