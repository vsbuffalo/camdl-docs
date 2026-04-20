"""True profile log-lik on (I(0), R(0)).

At each (I(0), R(0)) grid cell, Nelder-Mead over (β, γ) to maximise the
PF log-lik (fixed seed so the PF objective is deterministic). This is a
*true* profile — unlike the earlier slice that fixed β, γ.

~15 × 15 grid × ~50 NM calls/cell at 500 particles each ≈ 7 min wall.
"""
import subprocess, tempfile, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"
MODEL = "boarding_school_sir_init.ir.json"
DATA  = ROOT/"data/in_bed.tsv"
GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

N_PARTICLES = 500
SEED = 42
I0_GRID = np.linspace(1, 15, 15)
R0_GRID = np.linspace(0, 300, 15)
# (β, γ) init for Nelder-Mead — fixed-IC MLE; a reasonable central guess
BETA_INIT, GAMMA_INIT = 1.9, 0.66
BETA_BOUNDS  = (0.5, 5.0)
GAMMA_BOUNDS = (0.1, 1.0)

def eval_ll(beta, gamma, I0, R_init):
    with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as f:
        f.write(f"beta = {beta}\ngamma = {gamma}\n")
        f.write(f"I0 = {I0}\nR_init = {R_init}\nN0 = 763\n")
        path = f.name
    try:
        res = subprocess.run(
            ["camdl", "pfilter", MODEL, "--params", path,
             "--data", str(DATA), "--particles", str(N_PARTICLES),
             "--seed", str(SEED), "--dt", "1.0"],
            capture_output=True, text=True, check=False, cwd=str(ROOT))
        return float([ln for ln in res.stdout.strip().splitlines() if ln.strip()][-1])
    except (ValueError, IndexError):
        return float("nan")

def profile_at(I0, R0, beta0=BETA_INIT, gamma0=GAMMA_INIT):
    # Clip initial point to bounds
    b0 = float(np.clip(beta0,  *BETA_BOUNDS))
    g0 = float(np.clip(gamma0, *GAMMA_BOUNDS))
    def neg_ll(x):
        b, g = x
        if not (BETA_BOUNDS[0] <= b <= BETA_BOUNDS[1]): return 1e6
        if not (GAMMA_BOUNDS[0] <= g <= GAMMA_BOUNDS[1]): return 1e6
        ll = eval_ll(b, g, I0, R0)
        return -ll if np.isfinite(ll) else 1e6
    res = minimize(neg_ll, x0=[b0, g0], method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 0.1, "maxiter": 60})
    return -res.fun, float(res.x[0]), float(res.x[1]), res.nfev

cache = OUT/"profile_i0_r0_true.npz"
if cache.exists():
    d = np.load(cache)
    lls = d["lls"]; betas = d["betas"]; gammas = d["gammas"]
    print("loaded cached profile")
else:
    lls    = np.full((len(I0_GRID), len(R0_GRID)), np.nan)
    betas  = np.full_like(lls, np.nan)
    gammas = np.full_like(lls, np.nan)
    t0 = time.time()
    # Warm-start each cell with the best (β, γ) found at its I(0)-left
    # neighbour. Cuts the NM call budget roughly in half.
    prev_b, prev_g = BETA_INIT, GAMMA_INIT
    for j, R0 in enumerate(R0_GRID):
        col_b, col_g = BETA_INIT, GAMMA_INIT  # restart init per column
        for i, I0 in enumerate(I0_GRID):
            ll, bmax, gmax, nfev = profile_at(I0, R0, col_b, col_g)
            lls[i, j] = ll; betas[i, j] = bmax; gammas[i, j] = gmax
            col_b, col_g = bmax, gmax
        done = (j+1) * len(I0_GRID)
        total = len(I0_GRID) * len(R0_GRID)
        print(f"  col j={j} (R0={R0:.0f}): {done}/{total}  elapsed {time.time()-t0:.0f}s")
    np.savez(cache, lls=lls, betas=betas, gammas=gammas,
             I0_grid=I0_GRID, R0_grid=R0_GRID)
    print(f"wrote {cache}")

ll_max = np.nanmax(lls)
i_star, j_star = np.unravel_index(np.nanargmax(lls), lls.shape)
print(f"profile max ll = {ll_max:.2f} at I(0) = {I0_GRID[i_star]:.2f}, "
      f"R(0) = {R0_GRID[j_star]:.1f} with β = {betas[i_star, j_star]:.3f}, "
      f"γ = {gammas[i_star, j_star]:.3f}")

fig, axes = plt.subplots(1, 2, figsize=(14, 5.2),
                         gridspec_kw={"wspace": 0.25})

# Left panel: true profile log-lik
ax = axes[0]
delta = np.clip(ll_max - lls, 0, 30)
im = ax.pcolormesh(I0_GRID, R0_GRID, delta.T, cmap="viridis_r",
                   shading="auto", vmin=0, vmax=30)
cs = ax.contour(I0_GRID, R0_GRID, (ll_max - lls).T,
                levels=[1.92, 4.61, 9.21],
                colors=["cornflowerblue","#e67e22","#c0392b"], linewidths=1.0)
ax.clabel(cs, fmt={1.92:"95%", 4.61:"99%", 9.21:"99.9%"}, fontsize=8)
ax.scatter([I0_GRID[i_star]], [R0_GRID[j_star]], marker="*", s=260,
           color="red", edgecolor="white", linewidth=1.0, zorder=5,
           label="profile max")
ax.scatter([5], [0], marker="s", s=90, color="#3498db", edgecolor="white",
           linewidth=1.2, zorder=5, label="fixed-IC baseline")
ax.axhline(258, color="#f1c40f", linestyle="--", linewidth=1.2,
           label="Avilov sero R(0) = 258")
ax.set_xlabel("I(0)"); ax.set_ylabel("R(0)")
ax.set_title(f"True profile log-lik on (I(0), R(0))  "
             f"[β, γ optimised at each cell]", fontsize=10)
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.spines[["top","right"]].set_visible(False)
fig.colorbar(im, ax=ax, label="Δ log-lik from profile max (clipped at 30)")

# Right panel: how the profile-optimal (β, γ) moves with (I(0), R(0))
ax = axes[1]
# Scatter with β as color, γ as size — show how (β, γ) compensates
sc = ax.scatter(I0_GRID[:, None].repeat(len(R0_GRID), axis=1).ravel(),
                np.tile(R0_GRID, len(I0_GRID)),
                c=betas.ravel(), cmap="plasma", s=40,
                edgecolor="white", linewidth=0.5)
ax.scatter([I0_GRID[i_star]], [R0_GRID[j_star]], marker="*", s=260,
           color="red", edgecolor="white", linewidth=1.0, zorder=5)
ax.set_xlabel("I(0)"); ax.set_ylabel("R(0)")
ax.set_title("Optimal β̂(I(0), R(0)) across the grid", fontsize=10)
ax.spines[["top","right"]].set_visible(False)
fig.colorbar(sc, ax=ax, label="β̂ at profile max")

fig.suptitle("True profile likelihood — 15×15 grid, β and γ optimised via "
             "Nelder-Mead at each cell (500-particle PF)", fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(FIG/"profile_i0_r0_true.png", bbox_inches="tight")
plt.close(fig)
print(f"wrote {FIG/'profile_i0_r0_true.png'}")
