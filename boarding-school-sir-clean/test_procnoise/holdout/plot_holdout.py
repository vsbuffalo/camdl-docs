"""Holdout comparison: fit on days 0-7, evaluate on days 8-13."""
import sys
from pathlib import Path
import numpy as np, matplotlib.pyplot as plt

HERE = Path(__file__).parent
BOOK = HERE.parent.parent
sys.path.insert(0, str(BOOK))
from styles import apply_rc

apply_rc()

# Numbers from camdl pfilter (5000 particles, seed 1 throughout)
results = [
    # (label, k, train_ll, holdout_ll)
    ("Poisson\nfixed IC",                 2, -35.83, -22.66),
    ("NegBin\nfixed IC",                  3, -32.95, -24.80),
    ("Poisson + β proc-noise\nσ_se profile MLE", 3, -33.35, -22.41),
]
names = [r[0] for r in results]
ks    = np.array([r[1] for r in results])
train = np.array([r[2] for r in results])
hold  = np.array([r[3] for r in results])
total = train + hold
colors = ["#bdc3c7", "#3498db", "#e67e22"]

fig, axes = plt.subplots(1, 3, figsize=(13, 4.0),
                         gridspec_kw={"wspace": 0.32})

for ax, vals, ttl, ylbl in [
    (axes[0], train, "Train log-lik\n(days 0–7, 8 obs)", "log-likelihood"),
    (axes[1], hold,  "Holdout log-lik\n(days 8–13, 6 obs) — out-of-sample", "log-likelihood"),
    (axes[2], total, "Total log-lik (train + holdout)", "log-likelihood"),
]:
    bars = ax.bar(range(3), vals, color=colors)
    ax.set_xticks(range(3)); ax.set_xticklabels(names, fontsize=8.5)
    ax.set_ylabel(ylbl)
    ax.set_title(ttl, fontsize=10, pad=4)
    ax.set_ylim(min(vals)-0.8, max(vals)+0.5)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.1, f"{v:.2f}", ha="center", fontsize=9)
    # Mark best with a red star
    best = int(np.argmax(vals))
    ax.scatter([best], [vals[best]+0.35], marker="*", s=200,
               color="red", edgecolor="white", linewidth=0.8, zorder=5)
    ax.spines[["top","right"]].set_visible(False)

fig.suptitle("Holdout evaluation — fit on days 0–7, score on days 8–13 (post-peak decline)",
             fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(HERE.parent/"figures/holdout_comparison.png", bbox_inches="tight")
plt.close(fig)
print("wrote figures/holdout_comparison.png")

# Print table
print()
print(f"{'model':<40} {'k':>3} {'train ll':>10} {'holdout ll':>12} {'total':>10}")
print("-" * 80)
for r, v in zip(results, total):
    print(f"{r[0].replace(chr(10),' '):<40} {r[1]:>3} {r[2]:>10.2f} {r[3]:>12.2f} {v:>10.2f}")
