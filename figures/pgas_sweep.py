"""Generate PGAS sweep figure: alternating parameter and trajectory updates."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from scipy.ndimage import gaussian_filter

np.random.seed(42)

fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(11, 4.2))
fig.subplots_adjust(wspace=0.35)

# --- Shared trajectory data ---
T = 80
t = np.arange(T)
# Base trajectory: SIR-ish infected curve
base = 120 * np.exp(-0.5 * ((t - 30) / 12) ** 2)
noise = np.cumsum(np.random.randn(T) * 2.5)
traj = base + noise
traj = np.clip(traj, 0, None)

# Renewed trajectory: same early, diverges mid, reconverges late
renewed = traj.copy()
# Renewal region: substeps 25-55 get redrawn
renew_start, renew_end = 25, 55
renewed[renew_start:renew_end] = traj[renew_start:renew_end] + (
    np.random.randn(renew_end - renew_start) * 8
    + 12 * np.sin(np.linspace(0, np.pi, renew_end - renew_start))
)
renewed = np.clip(renewed, 0, None)
# Smooth the join points
for i in range(3):
    renewed[renew_start + i] = (
        traj[renew_start + i] * (3 - i) / 3
        + renewed[renew_start + i] * i / 3
    )
    renewed[renew_end - 1 - i] = (
        traj[renew_end - 1 - i] * (3 - i) / 3
        + renewed[renew_end - 1 - i] * i / 3
    )

# Colors
C_TRAJ = "#1a1a1a"
C_TRAJ_LIGHT = "#999999"
C_RENEWED = "#e74c3c"
C_KEPT = "#1a1a1a"
C_PARAM_OLD = "#3366cc"
C_PARAM_NEW = "#e74c3c"
C_CONTOUR = "#3366cc"
C_BG = "#f8f8f8"

# ============================================================
# LEFT PANEL: Fix trajectory, update parameters (NUTS step)
# ============================================================
ax = ax_left

# -- Trajectory subplot (top portion) --
# Draw the fixed trajectory prominently
ax.plot(t, traj, color=C_TRAJ, linewidth=2.0, zorder=5)

# Light "ghost" trajectories in background to suggest this is one of many possible
for i in range(4):
    ghost = traj + np.cumsum(np.random.randn(T) * 1.5) + np.random.randn() * 10
    ghost = np.clip(ghost, 0, None)
    ax.plot(t, ghost, color="#cccccc", linewidth=0.6, alpha=0.5, zorder=1)

# Pin icon: small markers to show trajectory is "locked"
ax.plot(
    [5, T - 5],
    [traj[5], traj[T - 5]],
    "s",
    color=C_TRAJ,
    markersize=5,
    zorder=6,
)

ax.set_xlim(-2, T + 2)
ax.set_ylim(-15, 170)
ax.set_xlabel("time ($t$)", fontsize=10)
ax.set_ylabel("infected ($I_t$)", fontsize=10)
ax.set_title(
    "Step 1: Fix trajectory, update $\\theta$",
    fontsize=11,
    fontweight="bold",
    pad=10,
)

# -- Parameter inset (bottom-right corner) --
inset = ax.inset_axes([0.55, 0.55, 0.42, 0.42])

# Draw log-posterior contours
xx, yy = np.meshgrid(np.linspace(-2, 2, 100), np.linspace(-2, 2, 100))
# Tilted elliptical density
rho = 0.6
zz = (xx ** 2 - 2 * rho * xx * yy + yy ** 2) / (2 * (1 - rho ** 2))
zz = np.exp(-zz)
zz = gaussian_filter(zz, sigma=3)
inset.contour(
    xx, yy, zz, levels=5, colors=C_CONTOUR, linewidths=0.8, alpha=0.5
)

# NUTS trajectory: a winding path through parameter space
nuts_t = np.linspace(0, 4 * np.pi, 60)
nuts_x = 0.3 + 0.8 * np.cos(nuts_t) * np.exp(-nuts_t / 15)
nuts_y = -0.2 + 0.8 * np.sin(nuts_t) * np.exp(-nuts_t / 15)
inset.plot(nuts_x, nuts_y, color=C_PARAM_OLD, linewidth=1.0, alpha=0.6)
inset.plot(nuts_x[0], nuts_y[0], "o", color=C_PARAM_OLD, markersize=6, zorder=10)
inset.plot(nuts_x[-1], nuts_y[-1], "*", color=C_PARAM_NEW, markersize=10, zorder=10)

inset.set_xlim(-2, 2)
inset.set_ylim(-2, 2)
inset.set_xticks([])
inset.set_yticks([])
inset.set_xlabel("$\\beta$", fontsize=9)
inset.set_ylabel("$\\gamma$", fontsize=9)
inset.set_facecolor("white")
inset.patch.set_alpha(0.9)
for spine in inset.spines.values():
    spine.set_edgecolor("#cccccc")

# Arrow from trajectory to inset
ax.annotate(
    "",
    xy=(0.55, 0.55),
    xytext=(0.45, 0.35),
    xycoords="axes fraction",
    textcoords="axes fraction",
    arrowprops=dict(
        arrowstyle="->",
        color="#888888",
        lw=1.2,
        connectionstyle="arc3,rad=0.2",
    ),
)
ax.text(
    0.35,
    0.30,
    "NUTS explores\n$p(\\theta \\mid x_{0:T}, y)$",
    transform=ax.transAxes,
    fontsize=8,
    color="#666666",
    ha="center",
)

# ============================================================
# RIGHT PANEL: Fix parameters, update trajectory (CSMC-AS)
# ============================================================
ax = ax_right

# Draw the old trajectory in gray
ax.plot(t, traj, color=C_TRAJ_LIGHT, linewidth=1.5, linestyle="--", alpha=0.6, zorder=2, label="previous $x_{0:T}$")

# Draw the renewed trajectory: kept segments in black, renewed in red
# Before renewal region
ax.plot(t[:renew_start + 1], renewed[:renew_start + 1], color=C_KEPT, linewidth=2.0, zorder=5)
# Renewal region
ax.plot(
    t[renew_start:renew_end],
    renewed[renew_start:renew_end],
    color=C_RENEWED,
    linewidth=2.5,
    zorder=5,
)
# After renewal region
ax.plot(t[renew_end - 1:], renewed[renew_end - 1:], color=C_KEPT, linewidth=2.0, zorder=5)

# Shade the renewal region
ax.axvspan(
    renew_start,
    renew_end,
    alpha=0.08,
    color=C_RENEWED,
    zorder=0,
)

# Annotations for kept/renewed
ax.annotate(
    "kept",
    xy=(12, renewed[12]),
    xytext=(8, 145),
    fontsize=9,
    color=C_KEPT,
    fontweight="bold",
    arrowprops=dict(arrowstyle="->", color=C_KEPT, lw=0.8),
)
ax.annotate(
    "renewed",
    xy=(40, renewed[40]),
    xytext=(42, 150),
    fontsize=9,
    color=C_RENEWED,
    fontweight="bold",
    arrowprops=dict(arrowstyle="->", color=C_RENEWED, lw=0.8),
)
ax.annotate(
    "kept",
    xy=(68, renewed[68]),
    xytext=(65, 80),
    fontsize=9,
    color=C_KEPT,
    fontweight="bold",
    arrowprops=dict(arrowstyle="->", color=C_KEPT, lw=0.8),
)

ax.set_xlim(-2, T + 2)
ax.set_ylim(-15, 170)
ax.set_xlabel("time ($t$)", fontsize=10)
ax.set_ylabel("infected ($I_t$)", fontsize=10)
ax.set_title(
    "Step 2: Fix $\\theta$, update trajectory",
    fontsize=11,
    fontweight="bold",
    pad=10,
)

# Parameter pin in top-left
ax.text(
    0.05,
    0.92,
    "$\\theta^* = (\\hat\\beta, \\hat\\gamma, \\ldots)$ fixed",
    transform=ax.transAxes,
    fontsize=8.5,
    color=C_PARAM_OLD,
    fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_PARAM_OLD, alpha=0.9),
)

ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

# ============================================================
# Circular arrow between panels
# ============================================================
# Top arrow: left → right
fig.patches.append(
    FancyArrowPatch(
        (0.48, 0.88),
        (0.535, 0.88),
        transform=fig.transFigure,
        arrowstyle="->,head_width=4,head_length=4",
        color="#444444",
        linewidth=1.5,
        connectionstyle="arc3,rad=-0.15",
    )
)
# Bottom arrow: right → left
fig.patches.append(
    FancyArrowPatch(
        (0.535, 0.18),
        (0.48, 0.18),
        transform=fig.transFigure,
        arrowstyle="->,head_width=4,head_length=4",
        color="#444444",
        linewidth=1.5,
        connectionstyle="arc3,rad=-0.15",
    )
)
fig.text(0.507, 0.92, "CSMC-AS", fontsize=7.5, ha="center", color="#444444")
fig.text(0.507, 0.12, "NUTS", fontsize=7.5, ha="center", color="#444444")

plt.savefig(
    "inference/figures/pgas_sweep.svg",
    bbox_inches="tight",
    dpi=150,
    transparent=False,
    facecolor="white",
)
plt.savefig(
    "inference/figures/pgas_sweep.png",
    bbox_inches="tight",
    dpi=200,
    transparent=False,
    facecolor="white",
)
print("Saved pgas_sweep.svg and pgas_sweep.png")
