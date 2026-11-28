"""Plot modular transfer at the sources and the held exterior clock."""

from fractions import Fraction as F
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def render():
    source = json.loads((HERE / "reference_arrival.json").read_bytes())
    flights = json.loads((HERE / "modular_arrival_flight.json").read_bytes())[
        "outcomes"
    ]
    arms = source["arms"]
    keys = ["opposed_serial", "same_serial", "shared"]
    q = np.array([np.mean(arms[k]["measured_modular_contrast"]) for k in keys]) * 1000
    serial = [
        144 * float(F(arms[k]["serial_acceptance"]) - F(25, 36)) * q[-1] for k in keys
    ]
    received = [
        np.mean(arms[k]["measured_arrival_change"])
        * np.mean(source["harmonic_normalization"])
        * 1000
        for k in keys
    ]
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.3, 3.05), layout="constrained")
    a.axhline(0, color="0.65", lw=0.65)
    a.plot([0, 3.65], [0, 3.65], color="0.55", lw=0.9)
    a.plot(q, serial, "x", color="#a95229", ms=6, label="Single contact pair")
    a.plot(
        q, received, "o", mfc="none", mec="#14745d", ms=7, label="Reference comparisons"
    )
    a.set(
        xlim=(0, 3.65),
        ylim=(-3.9, 4.1),
        xlabel="Measured source charge (millinats)",
        ylabel="Charge inferred at receipt (millinats)",
        title="Which source changes arrive?",
    )
    a.legend(loc="lower right", frameon=False, fontsize=8)
    b.axhline(0, color="0.65", lw=0.65)
    for direction, c in (("forward", "#14745d"), ("reversed", "#584b89")):
        pred = (
            np.array(
                [
                    flights[direction + "_" + k]["source_only_prediction"][-1]
                    for k in keys
                ]
            )
            * 1e6
        )
        held = (
            np.array(
                [
                    flights[direction + "_" + k]["observed_position_and_flight_change"][
                        -1
                    ]
                    for k in keys
                ]
            )
            * 1e6
        )
        xx = np.r_[0, q]
        b.fill_between(
            xx, np.r_[0, pred[:, 0]], np.r_[0, pred[:, 1]], color=c, alpha=0.22
        )
        b.errorbar(
            q,
            held.mean(axis=1),
            yerr=np.ptp(held, axis=1) / 2,
            fmt="o",
            ms=4.5,
            c=c,
            capsize=2,
            label=direction.capitalize() + " dispatch",
        )
    b.set(
        xlim=(0, 3.65),
        xlabel="Measured source charge (millinats)",
        ylabel="Clock change (millionths)",
        title="Clock after erasure",
    )
    b.legend(loc="upper left", frameon=False, fontsize=8)
    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.savefig(
        HERE / "reference_charge_transfer.pdf",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


if __name__ == "__main__":
    render()
