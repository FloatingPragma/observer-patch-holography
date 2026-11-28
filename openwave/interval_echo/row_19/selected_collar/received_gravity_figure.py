"""Plot the transported agreement signal and the independently calibrated probe responses."""

from pathlib import Path
from fractions import Fraction as F
import json

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def render():
    route = json.loads((HERE / "routed_contacts.json").read_bytes())
    row = next(
        r for r in route["profiles"] if (r["stations"], r["stop_bits"]) == (20, 5)
    )
    flight = json.loads((HERE / "routed_contact_flight.json").read_bytes())[
        "serial_agreement"
    ]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(8.1, 2.75), layout="constrained")
    shells = row["shells"]
    hops = [s["hops"] for s in shells]
    ax.plot(
        hops,
        [1000 * s["mean_odd"] / 144 for s in shells],
        color="#166254",
        lw=1.7,
        label="Finite route law",
    )
    ax.scatter(
        hops,
        [1000 * s["odd_counts"] / (s["sites"] * route["starts"] * 144) for s in shells],
        facecolors="white",
        edgecolors="#166254",
        zorder=3,
        s=32,
        label="Recorded visits",
    )
    ax.set(
        xticks=hops,
        xlabel="Hops from dispatch",
        ylabel="Excess writes per 1,000 trials",
        title="Where the agreement is recorded",
    )
    ax.legend(frameon=False, fontsize=8)
    data = flight["ordinary_doses"]["ordinary_one"]
    for j, (label, color, marker) in enumerate(
        [("Flight at 1/64", "#855126", "s"), ("Flight at 1/32", "#275882", "o")]
    ):
        ends = 100 * np.array(data["finite_dose_intervals"][j])
        mids = ends.mean(axis=1)
        bx.errorbar(
            mids,
            np.arange(3) + 0.15 * (j - 0.5),
            xerr=(ends[:, 1] - ends[:, 0]) / 2,
            fmt=marker,
            ms=4,
            c=color,
            capsize=2,
            label=label,
        )
    bx.scatter(
        [100 * float(F(k)) for k in data["required_fractions"]],
        np.arange(3),
        marker="|",
        c="black",
        s=60,
        label="At release",
    )
    bx.set(
        yticks=range(3),
        yticklabels=["Probe 1", "Probe 2", "Clock"],
        xlabel="Equivalent ordinary deposit (%)",
        title="What mass the instruments infer",
        ylim=(-0.35, 2.8),
        xlim=(0, 5.7),
    )
    bx.legend(
        frameon=False,
        fontsize=7,
        ncols=3,
        loc="upper center",
        columnspacing=0.7,
        handlelength=1.2,
    )
    for panel in (ax, bx):
        panel.spines[["top", "right"]].set_visible(False)
        panel.tick_params(labelsize=9)
        panel.title.set_fontsize(11)
    fig.savefig(
        HERE / "received_gravity_comparison.pdf",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


if __name__ == "__main__":
    render()
