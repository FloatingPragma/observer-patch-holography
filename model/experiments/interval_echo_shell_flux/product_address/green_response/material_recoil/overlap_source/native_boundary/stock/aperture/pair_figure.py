"""The contact-field distribution and the held mobile-probe comparison."""

from fractions import Fraction as F
import json
import matplotlib.pyplot as plt
import numpy as np

from .pair_work import HERE, SHOTS


def run():
    field = json.loads((HERE / "pair_field.json").read_bytes())
    held = json.loads((HERE / "pair_transfer.json").read_bytes())["ordinary_transfer"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(8, 3.2), layout="constrained")
    colors = ["#215775", "#aa492d"]
    for i, (name, row) in enumerate(field["spectrometers"]["1"].items()):
        left.bar(
            np.arange(5) + (i - 0.5) * 0.37,
            100 * np.array(row["counts"]) / SHOTS,
            width=0.36,
            color=colors[i],
            label=name.capitalize(),
        )
    left.set(
        xticks=np.arange(5),
        xticklabels=field["signed_levels"],
        ylabel="Preparations (%)",
        xlabel="Signed contact-field difference",
    )
    left.legend(frameon=False, fontsize=8)
    for j, t in enumerate(("1", "2")):
        v = field["independent_static_field"][t]
        measured = held[t]["rational_arm_arrivals"] / SHOTS
        error = measured - v["measured_tuned_arm_arrival_lower_99"]
        right.errorbar(
            j,
            100 * measured,
            yerr=100 * error,
            fmt="o",
            capsize=4,
            color=colors[1],
            label="Matched-energy preparation" if not j else None,
        )
        right.hlines(
            100 * float(F(v["any_independent_static_field_arrival_upper"])),
            j - 0.16,
            j + 0.16,
            color=colors[0],
            lw=3,
            label="Independent static field ceiling" if not j else None,
        )
    right.set(
        xlim=(-0.5, 1.5),
        xticks=[0, 1],
        xticklabels=["Single", "Double"],
        ylim=(29, 35.5),
        ylabel="Arrivals (%)",
        xlabel="Ordinary test load",
    )
    right.legend(frameon=False, fontsize=8, loc="upper center")
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(axis="y", alpha=0.18)
    fig.savefig(
        HERE / "pair_field.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)


if __name__ == "__main__":
    run()
