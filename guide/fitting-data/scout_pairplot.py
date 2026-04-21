"""Scout pair plot (lower-triangle chain trajectories + log-lik column)
for the v2 Poisson free-IC fit, using the shared plot_pair_grid helper."""
from pathlib import Path
from plot_traces import load_chains, best_chain_by_final_ll, plot_pair_grid

ROOT = Path(__file__).parent
fit_dir = next((ROOT/"results/fits").glob("fit_free_init_v2-*/real/fit_42"))

scout = load_chains(fit_dir/"scout")
best = best_chain_by_final_ll(scout)
print(f"best scout chain: {best}")

plot_pair_grid(
    scout,
    params=["beta", "gamma", "I0", "R_init", "if2_perturbed_loglik"],
    labels=["β", "γ", "I(0)", "R(0)", "log-lik"],
    best_chain=best,
    dll_filter=25.0,
    n_bins=14,
    ll_axis_dll=15,   # log-lik axes clipped to [best − 15, best + 2]
    outpath=ROOT/"figures/scout_pairplot.png",
    title=(f"Scout pair plot — v2 Poisson free-IC fit "
           f"(128 chains × 1000 iter, Δll < 25 from best)"),
)
print("wrote figures/scout_pairplot.png")
