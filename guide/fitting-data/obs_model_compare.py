"""Compare observation models — Poisson vs NegBin — at their respective MLEs.

Builds a side-by-side PPC figure: left panel Poisson, right panel NegBin.
Both fits use the same dynamics (chain-binomial SIR) and the same real data,
differ only in the observation layer.
"""
import subprocess, shlex
from pathlib import Path
import polars as pl, matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"

GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

POIS_MLE = next((ROOT/"results/fits").glob("fit-*/real/fit_42/refine/mle_params.toml"))
NB_MLE   = next((ROOT/"results/fits").glob("fit_negbin-*/real/fit_42/refine/mle_params.toml"))

def read_ll(mle_path: Path) -> float:
    for ln in mle_path.read_text().splitlines():
        if ln.strip().startswith("log_likelihood"):
            return float(ln.split("=")[1].strip())
    return float("nan")

def read_params(mle_path: Path) -> dict:
    d = {}
    for ln in mle_path.read_text().splitlines():
        s = ln.split("#",1)[0].strip()
        if s.startswith("[") or "=" not in s: continue
        k, v = s.split("=", 1)
        try: d[k.strip()] = float(v.strip())
        except: pass
    return d

# PPCs at each MLE — simulate forward under the model's own observation layer.
# Poisson uses the plain SIR model; NegBin uses the negbin model. Backend
# auto-match takes care of chain_binomial consistency.
PPC_POIS = OUT/"ppc_poisson.tsv"
PPC_NB   = OUT/"ppc_negbin.tsv"

subprocess.run(shlex.split(
    f"camdl simulate boarding_school_sir.ir.json --params {POIS_MLE} "
    f"--replicates 200 --seed 1 --scenario baseline --obs {PPC_POIS}"),
    check=True, capture_output=True, cwd=str(ROOT))
subprocess.run(shlex.split(
    f"camdl simulate boarding_school_sir_negbin.ir.json --params {NB_MLE} "
    f"--replicates 200 --seed 1 --scenario baseline --obs {PPC_NB}"),
    check=True, capture_output=True, cwd=str(ROOT))

obs = pl.read_csv(ROOT/"data/in_bed.tsv", separator="\t")
pois = pl.read_csv(PPC_POIS, separator="\t")
nb   = pl.read_csv(PPC_NB,   separator="\t")

pp = read_params(POIS_MLE); pp_ll = read_ll(POIS_MLE)
np_ = read_params(NB_MLE);  np_ll = read_ll(NB_MLE)

def ribbon(df, ax, color, label):
    q = df.group_by("time").agg(
        pl.col("in_bed").quantile(0.05).alias("lo"),
        pl.col("in_bed").quantile(0.5).alias("med"),
        pl.col("in_bed").quantile(0.95).alias("hi"),
    ).sort("time")
    ax.fill_between(q["time"], q["lo"], q["hi"], color=color, alpha=0.22, linewidth=0, edgecolor="none")
    ax.plot(q["time"], q["med"], color=color, linewidth=1.6, label=label)

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), sharey=True)

ribbon(pois, axes[0], "#3366cc",
       f"PPC envelope\nβ̂={pp['beta']:.3f}, γ̂={pp['gamma']:.3f}\nll = {pp_ll:.2f}")
axes[0].scatter(obs["time"], obs["in_bed"], color="black", s=30, zorder=5, label="observed")
axes[0].set_title("Poisson obs\n(likelihood = poisson(rate = prevalence(I)))",
                  fontsize=10, loc="left")
axes[0].set_xlabel("Day"); axes[0].set_ylabel("Boys in bed")
axes[0].legend(frameon=False, fontsize=9, loc="upper right")
axes[0].spines[["top","right"]].set_visible(False)

ribbon(nb, axes[1], "#c0392b",
       f"PPC envelope\nβ̂={np_['beta']:.3f}, γ̂={np_['gamma']:.3f}, k̂={np_['k']:.1f}\nll = {np_ll:.2f}")
axes[1].scatter(obs["time"], obs["in_bed"], color="black", s=30, zorder=5, label="observed")
axes[1].set_title("NegBin obs\n(likelihood = neg_binomial(mean = prevalence(I), r = k))",
                  fontsize=10, loc="left")
axes[1].set_xlabel("Day")
axes[1].legend(frameon=False, fontsize=9, loc="upper right")
axes[1].spines[["top","right"]].set_visible(False)

fig.suptitle("Observation model comparison — same SIR dynamics, same real data",
             fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(FIG/"obs_model_compare.png", bbox_inches="tight"); plt.close(fig)
print(f"wrote {FIG/'obs_model_compare.png'}")

# summary row
summary = pl.DataFrame([
    {"obs_model": "Poisson", "beta": pp["beta"], "gamma": pp["gamma"],
     "k": float("nan"), "ll": pp_ll, "n_params": 2, "aic": -2*pp_ll + 2*2,
     "bic": -2*pp_ll + 2*((14)**0.0) + 2 * __import__('math').log(14)},
    {"obs_model": "NegBin",  "beta": np_["beta"], "gamma": np_["gamma"],
     "k":     np_["k"],      "ll": np_ll,  "n_params": 3, "aic": -2*np_ll + 2*3,
     "bic": -2*np_ll + 3 * __import__('math').log(14)},
])
summary.write_csv(OUT/"obs_model_compare.tsv", separator="\t")
print(summary)
