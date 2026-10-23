"""Plot the captured passing executions, with no simulated timing inputs."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import HERE, DEFAULT_OUTPUT
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=HERE / "data" / "corrected-synthetic" / "measurements.json")
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT / "figures")
args = parser.parse_args()
root = args.output_dir.resolve()
root.mkdir(parents=True, exist_ok=True)
data=json.loads(args.input.read_text())
assert len(data)==20
x=[r["iteration"] for r in data]
plt.rcParams.update({"font.size":7,"font.family":"DejaVu Sans","pdf.fonttype":42})
fig,axes=plt.subplots(2,2,figsize=(3.5,2.7),layout="constrained")
axes[0,0].plot([r["tests"] for r in data],[1000*r["driver_s"] for r in data],"o-",ms=2,lw=.8,color="#537dae")
axes[0,0].set(title="(a) Measured driver time",xlabel="Collected tests",ylabel="Time (ms)")
axes[0,1].stackplot(x,[1000*r["preparation_s"] for r in data],
                    [1000*r["prioritization_s"] for r in data],
                    labels=["Preparation","Ranking"],colors=["#c77853","#65856b"],alpha=.85)
axes[0,1].set(title="(b) FAST components",xlabel="Mutation step",ylabel="Time (ms)")
axes[0,1].legend(fontsize=5.5,frameon=False,loc="upper left")
axes[1,0].plot(x,[1000*r["driver_s"]/r["tests"] for r in data],"o-",ms=2,lw=.8,color="#8a75a5")
axes[1,0].set(title="(c) Driver per test",xlabel="Mutation step",ylabel="Time (ms/test)")
axes[1,1].plot(x,[1000*r["prioritization_s"] for r in data],"o-",ms=2,lw=.8,color="#65856b")
axes[1,1].set(title="(d) Ranking only",xlabel="Mutation step",ylabel="Time (ms)")
for ax in axes.flat:
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(alpha=.15)
    ax.tick_params(labelsize=6)
    ax.set_ylim(bottom=0)
fig.savefig(root/"benchmark_grid.pdf",bbox_inches="tight",pad_inches=.02)
fig.savefig(root/"benchmark_grid.png",dpi=300,bbox_inches="tight",pad_inches=.02)
plt.close(fig)
