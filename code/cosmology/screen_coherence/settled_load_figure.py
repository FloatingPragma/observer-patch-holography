"""Plot measured correlations and regional transfers from the settled-load census."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(receipt, output):
    data = json.loads(receipt.read_text())
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.4, 4), layout="constrained")
    rows = ("initial_load", "settled_load", "change")
    x = np.arange(3)
    for group, label, offset, marker, color in (
        ("contact_pairs", "Direct contact", -0.11, "o", "#21618d"),
        ("noncontact_control", "Paired without contact", 0.11, "s", "#b14b18"),
    ):
        values = [data[group][row]["correlation"] for row in rows]
        left.scatter(
            x + offset, values, label=label, color=color, marker=marker, zorder=3
        )
    left.axhline(0, color="0.65", linewidth=0.8)
    left.set(
        xticks=x,
        xticklabels=["Initial load", "Settled load", "Load change"],
        ylabel="Correlation between two observers",
        title="Agreement across observers",
        xlim=(-0.45, 2.45),
        ylim=(-0.39, 0.10),
    )
    left.legend(loc="lower left", fontsize=9, frameon=False)
    regions = data["regions"]
    size = np.array([r["observers_per_region"] for r in regions])
    measured = [r["variance_change_total_per_member"] for r in regions]
    right.loglog(
        size,
        [r["mean_boundary_links"] / (2 * n) for r, n in zip(regions, size)],
        color="0.45",
        label="Independent fair preparation",
    )
    right.scatter(size, measured, color="#21618d", zorder=3, label="Archived regions")
    right.set(
        xlabel="Observers per region",
        ylabel="Variance of total change / observers",
        title="Exchange across a region's boundary",
    )
    right.legend(loc="upper right", fontsize=9, frameon=False)
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.9", linewidth=0.65)
        ax.set_axisbelow(True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    plot(args.receipt, args.output)
