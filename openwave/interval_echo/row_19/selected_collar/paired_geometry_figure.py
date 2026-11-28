"""Render two independent observer comparisons from the archived readings."""

from fractions import Fraction
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def draw():
    mass = json.loads((HERE / "curved_shells.json").read_bytes())
    balls = json.loads((HERE / "paired_balls.json").read_bytes())
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.0, 3.35), layout="constrained")
    teal, orange = "#12675a", "#a34420"
    calls = mass["regimes"]["complete"]
    cut = np.sqrt(1 - np.mean(mass["cosine_shell_locations"], axis=1) ** 2)
    surface = [v["population_mass"][0] for v in calls]
    left.step(
        [*cut, 1],
        [*surface, surface[-1]],
        where="post",
        color="0.38",
        label="Address census",
    )
    for row in calls:
        r = np.sqrt(np.mean(row["squared_areal_radius"]))
        lo, hi = row["local_detector_normalized_mass"]
        left.errorbar(
            r,
            (lo + hi) / 2,
            yerr=(hi - lo) / 2,
            fmt="o",
            markersize=4,
            capsize=3,
            color=teal,
        )
    left.plot([], [], "o", color=teal, markersize=4, label="Local clock and energy")
    left.set(
        xlabel="Areal radius / screen radius",
        ylabel="Enclosed mass / mass unit",
        ylim=(-0.025, 0.78),
        xlim=(-0.025, 1.025),
        title="Mass inferred from aperture records",
    )
    left.legend(fontsize=8, loc="upper left", frameon=False)
    rows = balls["metric_results"]["fixed"]
    index = np.arange(len(rows))
    for key, name, offset, marker, color in (
        ("native", "Original word", -0.1, "o", teal),
        ("reversed_second", "Reordered word", 0.1, "s", orange),
    ):
        values = [
            1000
            * float(Fraction(row["readings"][key]["finite_dose"]["exact_changes"]["A"]))
            for row in rows
        ]
        right.scatter(
            index + offset, values, label=name, color=color, marker=marker, s=20
        )
    right.axhline(0, color=".35", linewidth=0.8)
    right.set_xticks(index, [str(row["center"]) for row in rows])
    right.set(
        xlabel="Center of timed region (marker ID)",
        ylabel="Boundary change (thousandths of a link)",
        title="Boundary at the same mean volume",
        ylim=(-1.07, 1.08),
    )
    right.legend(fontsize=8, loc="lower left", frameon=False, ncol=2, columnspacing=0.7)
    for ax in (left, right):
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(axis="y", color=".9", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.title.set_fontsize(11)
        ax.tick_params(labelsize=8)
    fig.savefig(
        HERE / "paired_geometry_comparison.pdf",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


if __name__ == "__main__":
    draw()
