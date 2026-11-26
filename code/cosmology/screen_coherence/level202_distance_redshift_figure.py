"""Plot observer-zero clock and radar measurements from the published JSON receipt."""

import argparse
from collections import Counter
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def diagram(receipt, output):
    data = json.loads(receipt.read_text())
    ink = "#153548"
    teal = "#04746d"
    orange = "#b35022"
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), layout="constrained")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.23)
        ax.set_axisbelow(True)
        ax.set(xlim=(0.5, 20.5), xticks=[2, 6, 10, 14, 18])
        ax.set_xlabel("Initial radar distance (light-ticks)")
    for ax, epoch, color in zip(axes[:2], ("0", "96"), (ink, teal), strict=True):
        rows = data["panels"][epoch]
        counts = Counter(
            (r["radar_distance"], r["z_pulse"])
            for r in rows
            if r["radar_distance"] is not None and r["z_pulse"] is not None
        )
        for (dist, z), n in counts.items():
            ax.scatter(
                dist, z, s=13 * n, color=color, edgecolors="white", linewidths=0.45
            )
        ax.set_title(
            f"Query at tick {epoch}\n{sum(counts.values())} paired records",
            loc="left",
            fontsize=11,
            color=color,
            weight="bold",
            pad=12,
        )
        ax.set(yticks=[0, 0.125, 0.25], ylim=(-0.052, 0.31))
    axes[0].set_ylabel("Received period / source period − 1")
    axes[1].set_yticklabels([])
    ax = axes[2]
    continuations = data["vacuum_comparison"]["matched_long_baseline"]
    counts = Counter(
        (
            r["initial_radar_distance"],
            r["later_roundtrip_ticks"] / 2,
            r["predicted"]["next_radar_distance"],
        )
        for r in continuations
    )
    for (d, d2, model), n in counts.items():
        ax.scatter(d, d2, s=13 * n, color=teal, edgecolors="white", linewidths=0.5)
        ax.scatter(d, model, s=17, marker="x", color=orange, linewidths=1)
    ax.set(yticks=[0, 8, 16, 24, 32, 40], ylim=(-1, 41))
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("Later radar distance (light-ticks)")
    ax.set_title(
        "One doubling time later\n95 matched beacons",
        loc="left",
        fontsize=11,
        color=teal,
        weight="bold",
        pad=12,
    )
    ax.legend(
        handles=[
            Line2D([], [], marker="o", ls="", color=teal, label="Replay"),
            Line2D([], [], marker="x", ls="", color=orange, label="Comoving de Sitter"),
        ],
        loc="upper left",
        frameon=False,
        fontsize=8,
    )
    fig.suptitle(
        "Level 202  |  Observer 0’s distance–redshift measurements",
        fontsize=17,
        weight="bold",
        color=ink,
    )
    fig.savefig(
        output, dpi=150, bbox_inches="tight", metadata={"Software": "Matplotlib"}
    )
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    diagram(args.receipt, args.output)
