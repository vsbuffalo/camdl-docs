"""Scout diagnostics for the v2 Poisson free-IC fit.

Two figures:
  - Trace plot (chains coloured by final clean log-lik, bound lines shown)
  - Rank-overlay plot (Vehtari et al. 2021) — easier to read for many
    chains; well-mixed chains show flat rank histograms, stuck chains
    show concentrated ones.
"""
from pathlib import Path
from plot_traces import load_chains, best_chain_by_final_ll, \
                        plot_traces, plot_rank_overlay

ROOT = Path(__file__).parent
fit_dir = next((ROOT/"results/fits").glob("fit_free_init_v2-*/real/fit_42"))

scout = load_chains(fit_dir/"scout")
best = best_chain_by_final_ll(scout)
print(f"best scout chain: {best}  (128 chains × {int(scout['iteration'].max())+1} iter)")

# v2 fit's estimate-block bounds — show as dashed reference lines
BOUNDS = {
    "beta":   (1.0, 4.5),
    "gamma":  (0.2, 0.9),
    "I0":     (1.0, 10.0),
    "R_init": (0.0, 250.0),
}

plot_traces(
    scout,
    params=["beta", "gamma", "I0", "R_init", "if2_perturbed_loglik"],
    labels=["β", "γ", "I(0)", "R(0)", "log-lik (if2-perturbed)"],
    # v2 R(0) bounds are [0, 250] — linear reads fine; symlog added
    # spurious negative ticks because it's a symmetric scale.
    log_params=set(),
    best_chain=best,
    color_mode="ll",
    bounds=BOUNDS,
    outpath=ROOT/"figures/traces_poisson_v2_scout.png",
    title=(f"Poisson free-IC v2 scout — 128 chains, "
           f"coloured by final clean log-lik (dotted = bounds)"),
)
print("wrote figures/traces_poisson_v2_scout.png")

plot_rank_overlay(
    scout,
    params=["beta", "gamma", "I0", "R_init"],
    labels=["β", "γ", "I(0)", "R(0)"],
    best_chain=best,
    outpath=ROOT/"figures/traces_poisson_v2_rank.png",
    title="Rank-overlay diagnostic — Poisson free-IC v2 scout",
)
print("wrote figures/traces_poisson_v2_rank.png")
