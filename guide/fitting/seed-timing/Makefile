# Rebuild pipeline for the seed-timing chapter. Run from this directory.
# Results land in results/ (git-ignored). See README.md for data provenance.
SEED ?= 42

.PHONY: help all sim survey if2 pmmh compare clean
help:
	@echo "make sim      regenerate synthetic seed_cases.tsv (camdl simulate at TRUTH)"
	@echo "make survey   likelihood-landscape survey + interactive HTML"
	@echo "make if2      IF2 (MLE) fits — these FAIL the gate (Method 1 section)"
	@echo "make pmmh     PMMH posteriors: informative, vague, const-rho, od"
	@echo "make compare  camdl compare across the 3 PMMH models (ELPD/CRPS/PIT)"
	@echo "make all      sim survey if2 pmmh compare"
	@echo "make clean    rm -rf results/*  (data fetched separately; see README)"

all: sim survey if2 pmmh compare

# Synthetic recovery data: TRUTH = beta .6 gamma .2 lambda 2 w 3 N0 5000 tau 30 rho .5 k 20
sim:
	camdl simulate models/seed_timing.camdl --backend chain_binomial --dt 1 --seed 11 \
	  --param beta=0.6 --param gamma=0.2 --param lambda=2.0 --param w=3.0 \
	  --param N0=5000 --param tau=30 --param rho=0.5 --param k=20 \
	  --obs-only data/seed_cases.tsv

survey:
	camdl survey models/seed_timing_report.camdl --fit fits/covid_wa_survey.toml \
	  --eval pfilter --eval-particles 500 --eval-replicates 5 --render \
	  --label "WA seed-timing landscape" --output results/surveys

if2:
	camdl fit run fits/covid_wa_if2.toml          --seed $(SEED) --stage scout --progress plain
	camdl fit run fits/covid_wa_if2_gammafix.toml --seed $(SEED) --stage scout --progress plain

pmmh:
	camdl fit run fits/covid_wa_pmmh.toml         --seed $(SEED) --stage posterior --progress plain
	camdl fit run fits/covid_wa_pmmh_vague.toml   --seed $(SEED) --stage posterior --progress plain
	camdl fit run fits/covid_wa_pmmh_constrho.toml --seed $(SEED) --stage posterior --progress plain
	camdl fit run fits/covid_wa_pmmh_od.toml      --seed $(SEED) --stage posterior --progress plain

# Generate prequential traces at each posterior mode, stage canonical artifacts
# to stable results/ paths (consumed by the chapter), and run `camdl compare`.
compare:
	python scripts/stage_and_compare.py

clean:
	rm -rf results/*
