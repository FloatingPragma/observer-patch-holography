"""Plot the measured temperature spectrum and the declared cold-fluid completions."""

import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).parent


def main():
    data = json.loads((HERE / "results.json").read_text())
    sky = data["measured"]
    ell = sky["ell"]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), layout="constrained")
    for ax, key, title in zip(
        axes,
        ("raw_at_effective_ell", "aligned_at_effective_ell"),
        ("Fixed cosmological background", "First peak matched"),
        strict=True,
    ):
        for case, color, dash, label in (
            ("ordinary_cold_control", "#187778", "-", "Standard cold density"),
            ("maximal_one_layer", "#a44819", "--", "Maximum one-layer density"),
        ):
            ax.plot(
                ell,
                data["variants"][case][key],
                dash,
                color=color,
                linewidth=1.6,
                label=label,
                zorder=2,
            )
        ax.errorbar(
            ell,
            sky["measured"],
            sky["marginal_error"],
            fmt="o",
            markersize=2.8,
            color="#222222",
            elinewidth=0.65,
            label="Planck temperature",
            zorder=3,
        )
        ax.set(
            xlim=(115, 1535),
            ylim=(0, 15000 if key.startswith("raw") else 6800),
            xlabel="Angular multipole ℓ",
            title=title,
        )
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(axis="y", color="0.9", linewidth=0.6)
    axes[0].set_ylabel("Temperature power $D_\\ell$  (µK²)")
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)
    fig.savefig(
        HERE / "comparison.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
