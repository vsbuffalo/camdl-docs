"""2D profile-likelihood grids for the NegBin fit. PF log-lik evaluated
on a regular grid; shows the ridge structure visually.

Grids:
  β × γ  — classic SIR identifiability ridge (R₀ = β/γ stays roughly
           constant along the ridge, outbreak shape weakly sensitive to
           β and γ individually when their ratio is preserved).
  β × k  — does NegBin dispersion trade off with transmission?

All other parameters fixed at the NegBin MLE.
"""
import subprocess, tempfile, time
from pathlib import Path
import numpy as np, polars as pl
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"
MODEL = "boarding_school_sir_negbin.ir.json"
DATA  = ROOT/"data/in_bed.tsv"

GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

# MLE to center on
MLE_PATH = next((ROOT/"results/fits").glob("fit_negbin-*/real/fit_42/refine/mle_params.toml"))
mle = {}
for ln in MLE_PATH.read_text().splitlines():
    s = ln.split("#",1)[0].strip()
    if s.startswith("[") or "=" not in s: continue
    k, v = s.split("=",1)
    try: mle[k.strip()] = float(v.strip())
    except: pass
print(f"NegBin MLE: β={mle['beta']:.3f}, γ={mle['gamma']:.3f}, k={mle['k']:.1f}")

N_GRID     = 20
N_PARTICLES = 1000
SEED        = 42

def eval_ll(params: dict) -> float:
    """Run camdl pfilter at params, return log-lik (stdout last line)."""
    with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as f:
        for k, v in params.items():
            f.write(f"{k} = {v}\n")
        path = f.name
    try:
        res = subprocess.run(
            ["camdl", "pfilter", MODEL, "--params", path,
             "--data", str(DATA), "--particles", str(N_PARTICLES),
             "--seed", str(SEED), "--dt", "1.0"],
            capture_output=True, text=True, check=False, cwd=str(ROOT))
        lines = [ln for ln in res.stdout.strip().splitlines() if ln.strip()]
        return float(lines[-1])
    except (ValueError, IndexError):
        return float("nan")

def run_grid(x_name, x_range, y_name, y_range, fixed, tag):
    xs = np.linspace(*x_range, N_GRID)
    ys = np.linspace(*y_range, N_GRID)
    lls = np.full((len(xs), len(ys)), np.nan)
    t0 = time.time()
    for i, xv in enumerate(xs):
        for j, yv in enumerate(ys):
            params = {**fixed, x_name: float(xv), y_name: float(yv)}
            lls[i, j] = eval_ll(params)
        done = (i + 1) * len(ys)
        total = len(xs) * len(ys)
        print(f"  [{tag}] {done}/{total}  elapsed {time.time()-t0:.0f}s")
    np.savez(OUT/f"profile_{tag}.npz", xs=xs, ys=ys, lls=lls,
             x_name=x_name, y_name=y_name, mle=mle)
    return xs, ys, lls

# Grid 1: β × γ (classic SIR ridge) — zoomed to where the basin lives
print("\n=== grid 1: β × γ (k fixed at MLE) ===")
xs1, ys1, lls1 = run_grid(
    "beta", (1.4, 3.0), "gamma", (0.4, 1.0),
    {"I0": 5, "N0": 763, "k": mle["k"]}, "beta_gamma")

# Grid 2: β × k (does dispersion trade with transmission?)
print("\n=== grid 2: β × k (γ fixed at MLE) ===")
xs2, ys2, lls2 = run_grid(
    "beta", (1.4, 3.0), "k", (5.0, 1500.0),
    {"I0": 5, "N0": 763, "gamma": mle["gamma"]}, "beta_k")

# Grid 3: γ × k (does dispersion trade with recovery?)
print("\n=== grid 3: γ × k (β fixed at MLE) ===")
xs3, ys3, lls3 = run_grid(
    "gamma", (0.4, 1.0), "k", (5.0, 1500.0),
    {"I0": 5, "N0": 763, "beta": mle["beta"]}, "gamma_k")

# Also load Poisson MLE for overlay (for the β×γ panel — k dimension absent)
POIS_MLE_PATH = next((ROOT/"results/fits").glob("fit-*/real/fit_42/refine/mle_params.toml"))
pois = {}
for ln in POIS_MLE_PATH.read_text().splitlines():
    s = ln.split("#",1)[0].strip()
    if s.startswith("[") or "=" not in s: continue
    k, v = s.split("=",1)
    try: pois[k.strip()] = float(v.strip())
    except: pass

# Plot all three — shared colorbar, all markers with legend
fig, axes = plt.subplots(1, 3, figsize=(17, 5.3),
                         gridspec_kw={"width_ratios": [1, 1, 1], "wspace": 0.25})

# Find global Δll range for shared colorbar
all_deltas = []
for lls in (lls1, lls2, lls3):
    ll_max = np.nanmax(lls)
    all_deltas.append(np.clip(ll_max - lls, 0, 30))

panels = [
    # xs, ys, lls, x_lab, y_lab, title,
    # refine_xy, grid_xy, pois_xy
    (xs1, ys1, lls1, "β", "γ",
     "β × γ  (k = {:.1f} fixed)".format(mle["k"]),
     (mle["beta"], mle["gamma"]), None, (pois["beta"], pois["gamma"])),
    (xs2, ys2, lls2, "β", "k",
     "β × k  (γ = {:.3f} fixed)".format(mle["gamma"]),
     (mle["beta"], mle["k"]), None, None),
    (xs3, ys3, lls3, "γ", "k",
     "γ × k  (β = {:.3f} fixed)".format(mle["beta"]),
     (mle["gamma"], mle["k"]), None, None),
]

# Fill in grid-max for each panel (replace the None slot)
for i in range(len(panels)):
    xs, ys, lls, x_lab, y_lab, title, refine_xy, _, pois_xy = panels[i]
    gi, gj = np.unravel_index(np.nanargmax(lls), lls.shape)
    grid_xy = (float(xs[gi]), float(ys[gj]))
    panels[i] = (xs, ys, lls, x_lab, y_lab, title, refine_xy, grid_xy, pois_xy)

im = None
for ax, (xs, ys, lls, x_lab, y_lab, title,
        refine_xy, grid_xy, pois_xy) in zip(axes, panels):
    ll_max = np.nanmax(lls)
    delta = np.clip(ll_max - lls, 0, 30)
    im = ax.pcolormesh(xs, ys, delta.T, cmap="viridis_r",
                       shading="auto", vmin=0, vmax=30)
    cs = ax.contour(xs, ys, (ll_max - lls).T,
                    levels=[1.92, 4.61, 9.21],
                    colors=["cornflowerblue", "#e67e22", "#c0392b"], linewidths=1.0)
    ax.clabel(cs, fmt={1.92: "95%", 4.61: "99%", 9.21: "99.9%"}, fontsize=8)
    ax.scatter([refine_xy[0]], [refine_xy[1]], marker="*", s=260,
               color="red", edgecolor="white", linewidth=1.0, zorder=6,
               label="refine MLE")
    ax.scatter([grid_xy[0]], [grid_xy[1]], marker="D", s=90,
               color="#2ecc71", edgecolor="white", linewidth=1.6, zorder=5,
               label="grid max")
    if pois_xy is not None:
        ax.scatter([pois_xy[0]], [pois_xy[1]], marker="s", s=90,
                   color="#3498db", edgecolor="white", linewidth=1.6, zorder=5,
                   label="Poisson MLE")
    ax.set_xlabel(x_lab); ax.set_ylabel(y_lab)
    ax.set_title(title, fontsize=10)
    ax.spines[["top","right"]].set_visible(False)

# Single shared colorbar in a dedicated axis so it doesn't overlap the right-most panel
cbar_ax = fig.add_axes([0.905, 0.18, 0.018, 0.68])
cbar = fig.colorbar(im, cax=cbar_ax,
                    label="Δ log-lik from panel max (clipped at 30)")
# Single legend — put it in the first panel only; markers are the same across
handles, labels = axes[0].get_legend_handles_labels()
# Replicate a clean legend below suptitle
fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
           fontsize=10, bbox_to_anchor=(0.5, 0.04))

fig.suptitle("2D profile log-likelihood grids — NegBin SIR, "
             f"{N_GRID}×{N_GRID} grid, {N_PARTICLES}-particle PF",
             fontsize=11, y=0.98)
fig.subplots_adjust(top=0.88, bottom=0.14, left=0.05, right=0.88, wspace=0.28)
fig.savefig(FIG/"profile_2d.png", bbox_inches="tight")
plt.close(fig)
print(f"\nwrote {FIG/'profile_2d.png'}")
