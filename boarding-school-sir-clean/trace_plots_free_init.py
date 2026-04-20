"""Scout traces for the v1 Poisson free-IC fits.

v1 = initial attempts with wider bounds and shorter scout budget. Both
the standard and ic_free variants had their scout fail the tail-Rhat
gate. Showing scout here (not refine) because scout is the stage that
didn't converge — refine only "ran" because those early attempts used
--allow-nonconverged-scout, a flag we no longer use per CLAUDE.md.
"""
from pathlib import Path
from plot_traces import load_chains, best_chain_by_final_ll, plot_traces

ROOT = Path(__file__).parent

# v1 estimate-block bounds — wider than v2
BOUNDS_V1 = {
    "beta":   (0.5, 5.0),
    "gamma":  (0.1, 1.0),
    "I0":     (1.0, 20.0),
    "R_init": (0.0, 400.0),
}

for tag, title in [("fit_free_init", "free-IC standard, v1"),
                   ("fit_ic_free",   "IC-free (y₁ conditional), v1")]:
    try:
        fit_dir = next((ROOT/"results/fits").glob(f"{tag}-*/real/fit_42"))
    except StopIteration:
        print(f"skip {tag}: no fit dir found")
        continue
    scout = load_chains(fit_dir/"scout")
    best  = best_chain_by_final_ll(scout)
    n_chains = scout["chain"].n_unique()
    n_iter   = int(scout["iteration"].max()) + 1
    print(f"{tag}: {n_chains} chains × {n_iter} iter, best chain {best}")
    plot_traces(
        scout,
        params=["beta", "gamma", "I0", "R_init", "if2_perturbed_loglik"],
        labels=["β", "γ", "I(0)", "R(0)", "log-lik (if2-perturbed)"],
        best_chain=best,
        color_mode="ll",
        bounds=BOUNDS_V1,
        outpath=ROOT/f"figures/traces_{tag}.png",
        title=(f"Scout traces — {title} "
               f"({n_chains} chains × {n_iter} iter, chain {best} = best final ll)"),
    )
    print(f"  wrote figures/traces_{tag}.png")
