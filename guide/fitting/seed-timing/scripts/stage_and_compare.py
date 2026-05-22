#!/usr/bin/env python3
"""Post-PMMH: generate prequential traces at each model's posterior mode, copy
canonical artifacts to stable results/ paths (consumed by the chapter), build
fits/compare.toml, and run `camdl compare`. Run from guide/fitting/seed-timing/
after `make pmmh`. Idempotent."""
import json, shutil, subprocess
from pathlib import Path

FITS = Path("results/fits")
MODELS = {
    "constrho": ("covid_wa_pmmh_constrho", "models/seed_timing_const_rho_prior.camdl",
                 {"N0": 1000000, "w": 5.0}, "constant rho"),
    "post":     ("covid_wa_pmmh", "models/seed_timing_report_prior.camdl",
                 {"N0": 1000000, "w": 5.0, "w_rep": 7.0, "t_rep": 37.0}, "rho(t) ramp"),
    "od":       ("covid_wa_pmmh_od", "models/seed_timing_report_od_prior.camdl",
                 {"N0": 1000000, "w": 5.0, "w_rep": 7.0, "t_rep": 37.0}, "rho(t) + overdispersion"),
    "vague":    ("covid_wa_pmmh_vague", "models/seed_timing_report_prior.camdl",
                 {"N0": 1000000, "w": 5.0, "w_rep": 7.0, "t_rep": 37.0}, "vague priors"),
}

def post(slug): return sorted(FITS.glob(f"{slug}-*"))[-1] / "real/fit_42/posterior"

Path("results/survey").mkdir(parents=True, exist_ok=True)
for tag, (slug, model, fixed, _) in MODELS.items():
    p = post(slug)
    shutil.copy(p / "draws.tsv", f"results/wa_{tag}_draws.tsv")
    shutil.copy(p / "pmmh_summary.json", f"results/wa_{tag}_summary.json")
    # prequential at the posterior mode (no --replicates: that overrides --save-prequential)
    mp = json.load(open(p / "pmmh_summary.json"))["map_params"]
    args = ["camdl", "pfilter", model, "--data", "data/covid_wa_growth.tsv",
            "--particles", "4000", "--dt", "1", "--seed", "1",
            "--save-prequential", str(p / "prequential")]
    for k, v in {**mp, **fixed}.items():
        args += ["--param", f"{k}={v}"]
    subprocess.run(args, capture_output=True, check=True)

sv = sorted(Path("results/surveys").rglob("landscape.html"))
if sv:
    shutil.copy(sv[-1], "results/survey/landscape.html")

# compare.toml over the three model-comparison arms (not vague)
def cfg_path(tag): return str(post(MODELS[tag][0]))
toml = 'baseline = "rho(t) ramp"\nmetrics = ["elpd", "crps", "pit_cov90"]\nformat = "md"\n'
for tag in ("constrho", "post", "od"):
    toml += f'\n[[model]]\nname = "{MODELS[tag][3]}"\npath = "{cfg_path(tag)}"\n'
Path("fits/compare.toml").write_text(toml)
out = subprocess.run(["camdl", "compare", "--config", "fits/compare.toml"],
                     capture_output=True, text=True)
Path("results/wa_compare.md").write_text(out.stdout)
print(out.stdout)
print("staged wa_{post,vague,constrho,od}_{draws.tsv,summary.json}, survey/landscape.html, wa_compare.md")
