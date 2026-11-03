"""Vẽ đường hội tụ kiểu Hình 12: best makespan (% so với nghiệm ban đầu) theo vòng lặp.

    python scripts/plot_history.py results/history_random.csv results/fig_history_random.pdf [en]
"""
import sys
import csv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EN = len(sys.argv) > 3 and sys.argv[3] == "en"
SERIES = {"TSR": ("#2a78d6", "TSR (with redistribute)" if EN else "TSR (có redistribute)"),
          "TSO": ("#eb6834", "TSO (original)" if EN else "TSO (gốc)")}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

src, out = sys.argv[1], sys.argv[2]
with open(src, newline="") as f:
    rows = list(csv.DictReader(f))
it = [int(r["iteration"]) for r in rows]

fig, ax = plt.subplots(figsize=(6, 3.6), dpi=200)
for key, (color, label) in SERIES.items():
    y = [float(r[key]) for r in rows]
    ax.plot(it, y, color=color, lw=2, label=label)
    ax.annotate(f"{key} {y[-1]:.2f}%", xy=(it[-1], y[-1]), xytext=(4, 0),
                textcoords="offset points", va="center", fontsize=8, color=INK)
ax.set_xlabel("Iteration" if EN else "Vòng lặp", color=MUTED)
ax.set_ylabel("Best makespan (% of initial)" if EN else "Best makespan (% nghiệm ban đầu)", color=MUTED)
ax.grid(axis="y", color=GRID, lw=0.8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=MUTED, labelsize=8)
ax.legend(frameon=False, fontsize=8, labelcolor=INK)
ax.set_xlim(0, it[-1] * 1.12)
fig.tight_layout()
fig.savefig(out)
print("saved", out)
