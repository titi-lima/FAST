"""Draw verified complete-suite process times and measured FAST driver work."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import HERE, DEFAULT_OUTPUT
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=HERE / "data" / "verified-summary.json")
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT / "figures")
args = parser.parse_args()
root = args.output_dir.resolve()
root.mkdir(parents=True, exist_ok=True)
projects = ["requests", "rich", "click", "attrs"]
raw=json.loads(args.input.read_text())
indexed={p["project"]:p for p in raw["projects"]}
data=[indexed[name] for name in projects]
assert all(p["valid"] for p in data)
baseline = [p["baseline_mean_s"] for p in data]
baseline_std = [p["baseline_std_s"] for p in data]
fast = [p["fast_mean_s"] for p in data]
fast_std = [p["fast_std_s"] for p in data]
driver = [p["warm_fast_driver_mean_s"] for p in data]
plt.rcParams.update({"font.size": 7, "font.family": "DejaVu Sans", "pdf.fonttype": 42})
fig, (left, right) = plt.subplots(1, 2, figsize=(3.5, 2.05))
fig.subplots_adjust(left=.11, right=.995, bottom=.24, top=.72, wspace=.62)
x = np.arange(4)
base_bars = left.bar(x - .18, baseline, .36, yerr=baseline_std, capsize=1.5,
         label="Baseline", color="#537dae", error_kw={"linewidth": .6})
fast_bars = left.bar(x + .18, fast, .36, yerr=fast_std, capsize=1.5,
         label="+ fast-tcp", color="#c77853", error_kw={"linewidth": .6})
left.set(ylabel="Wall-clock time (s)", ylim=(0,max(baseline + fast)*1.5))
left.set_title("(a) Mean process time", fontsize=7, pad=25)
left.set_xticks(x, projects, rotation=30, ha="right", fontsize=6.4)
left.legend(fontsize=5.8, frameon=False, loc="lower center",
            bbox_to_anchor=(.5,1.01), ncol=2, handlelength=1.1,
            columnspacing=.8, handletextpad=.3, borderaxespad=0)
# Direct values keep the short rich/click bars interpretable on the same linear
# zero-based scale. Labels are staggered above each pair, clear of error bars.
for bars, values, deviations, offset in [
        (base_bars, baseline, baseline_std, 12),
        (fast_bars, fast, fast_std, 3)]:
    for bar, value, deviation in zip(bars, values, deviations):
        left.annotate(f"{value:.2f}",
                      (bar.get_x()+bar.get_width()/2, value+deviation),
                      xytext=(0,offset), textcoords="offset points",
                      ha="center", fontsize=5.2)
bars = right.bar(x, driver, width=.65, color="#65856b")
right.set(ylabel="Driver time (s)", ylim=(0, max(driver)*1.3))
right.set_title("(b) FAST driver time", fontsize=7, pad=25)
right.set_xticks(x, projects, rotation=30, ha="right", fontsize=6.4)
for bar, value in zip(bars, driver):
    right.annotate(f"{value:.3f}", (bar.get_x()+bar.get_width()/2, value),
                   xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6.5)
for ax in [left, right]:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=.15)
    ax.set_axisbelow(True)
fig.savefig(root / "case_study_bar.pdf", bbox_inches="tight", pad_inches=.02)
fig.savefig(root / "case_study_bar.png", dpi=300, bbox_inches="tight", pad_inches=.02)
plt.close(fig)
