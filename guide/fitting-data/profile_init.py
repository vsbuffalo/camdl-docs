"""2D profile log-likelihood grid on (I(0), R(0)) at the free-init MLE's
(β, γ). Shows the initial-condition identifiability ridge."""
import subprocess, tempfile, time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

ROOT = Path(__file__).parent
FIG  = ROOT/"figures"
OUT  = ROOT/"output"
MODEL = "boarding_school_sir_init.ir.json"
DATA  = ROOT/"data/in_bed.tsv"
GREY = "#555555"
mpl.rcParams.update({"figure.dpi":150,"font.size":10,"axes.edgecolor":GREY,
    "axes.labelcolor":GREY,"xtick.color":GREY,"ytick.color":GREY,"text.color":GREY})

# Anchor β and γ at the free-init MLE
MLE_PATH = next((ROOT/"results/fits").glob("fit_free_init-*/real/fit_42/refine/mle_params.toml"))
mle = {}
for ln in MLE_PATH.read_text().splitlines():
    s = ln.split("#",1)[0].strip()
    if s.startswith("[") or "=" not in s: continue
    k, v = s.split("=", 1)
    try: mle[k.strip()] = float(v.strip())
    except: pass
print(f"free-init MLE: β={mle['beta']:.3f}, γ={mle['gamma']:.3f}, "
      f"I(0)={mle['I0']:.2f}, R(0)={mle['R_init']:.1f}")

# External reference: Avilov et al. 2024 sero-based R(0) ≈ 258
AVILOV_R0 = 258

N_GRID, N_PARTICLES, SEED = 20, 1000, 42
I_RANGE = (1.0, 20.0)
R_RANGE = (0.0, 400.0)

def eval_ll(params):
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
        return float([ln for ln in res.stdout.strip().splitlines() if ln.strip()][-1])
    except (ValueError, IndexError):
        return float("nan")

cache = OUT/"profile_init.npz"
if cache.exists():
    d = np.load(cache)
    xs, ys, lls = d["xs"], d["ys"], d["lls"]
    print("using cached grid")
else:
    xs = np.linspace(*I_RANGE, N_GRID)
    ys = np.linspace(*R_RANGE, N_GRID)
    lls = np.full((len(xs), len(ys)), np.nan)
    t0 = time.time()
    for i, I0 in enumerate(xs):
        for j, R0 in enumerate(ys):
            lls[i, j] = eval_ll({"beta": mle["beta"], "gamma": mle["gamma"],
                                 "I0": float(I0), "R_init": float(R0),
                                 "N0": 763})
        if (i+1) % 5 == 0:
            print(f"  {(i+1)*len(ys)}/{len(xs)*len(ys)}  {time.time()-t0:.0f}s")
    np.savez(cache, xs=xs, ys=ys, lls=lls)

ll_max = np.nanmax(lls)
i, j = np.unravel_index(np.nanargmax(lls), lls.shape)
print(f"grid max: ll = {ll_max:.2f} at I(0) = {xs[i]:.2f}, R(0) = {ys[j]:.1f}")

fig, ax = plt.subplots(figsize=(8.5, 5.5))
delta = np.clip(ll_max - lls, 0, 30)
im = ax.pcolormesh(xs, ys, delta.T, cmap="viridis_r", shading="auto",
                   vmin=0, vmax=30)
cs = ax.contour(xs, ys, (ll_max - lls).T, levels=[1.92, 4.61, 9.21],
                colors=["cornflowerblue","#e67e22","#c0392b"], linewidths=1.0)
ax.clabel(cs, fmt={1.92:"95%", 4.61:"99%", 9.21:"99.9%"}, fontsize=8)

# markers
ax.scatter([mle["I0"]], [mle["R_init"]], marker="*", s=260, color="red",
           edgecolor="white", linewidth=1.0, zorder=7, label="free-init MLE")
ax.scatter([xs[i]], [ys[j]], marker="D", s=110, color="#2ecc71",
           edgecolor="white", linewidth=1.4, zorder=6, label="grid max")
ax.axhline(AVILOV_R0, color="#f1c40f", linestyle="--", linewidth=1.3,
           zorder=5, label=f"Avilov et al. 2024 sero R(0) = {AVILOV_R0}")
# fixed-IC choice (I(0)=5, R(0)=0) for reference
ax.scatter([5], [0], marker="s", s=90, color="#3498db",
           edgecolor="white", linewidth=1.2, zorder=6,
           label="fixed-IC scenario (I(0)=5, R(0)=0)")

ax.set_xlabel("I(0)   initial infectious")
ax.set_ylabel("R(0)   initial recovered / immune")
ax.set_title(f"2D profile log-lik on (I(0), R(0))  "
             f"(β = {mle['beta']:.3f}, γ = {mle['gamma']:.3f} fixed at free-init MLE)",
             fontsize=10)
ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.spines[["top","right"]].set_visible(False)
fig.colorbar(im, ax=ax, label="Δ log-lik from grid max (clipped at 30)")
fig.tight_layout()
fig.savefig(FIG/"profile_init.png", bbox_inches="tight")
plt.close(fig)
print(f"wrote {FIG/'profile_init.png'}")
