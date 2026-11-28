"""Compare information lost to thinning with material flight on the routed sources."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import measure

HERE = Path(__file__).parent


def draw():
    r = json.loads((HERE / "results.json").read_text())
    states = r["reference"][1]["states"]
    flight = r["probes"]["flight"]
    fig, (access, motion) = plt.subplots(
        1, 2, figsize=(7.1, 3.05), layout="constrained"
    )
    family = states["FG0"]["thinning"]
    for key, color, label in (
        ("unflagged_loss", "#176864", "Missing flag discarded"),
        ("flagged_loss", "#305b8a", "Missing flag retained"),
    ):
        access.plot(
            [x["p"] for x in family],
            [x[key] for x in family],
            color=color,
            marker="o",
            ms=3,
            label=label,
        )
    access.set(
        xlabel="Probability of erasing the memory",
        ylabel="Relative entropy lost (nats)",
        title="Same loss across routed sources",
        xlim=(0, 1),
    )
    access.legend(frameon=False, fontsize=8, loc="upper left")
    data = np.array([flight[f"FG{i}"]["mean"][2] for i in range(len(measure.ANGLES))])
    theta = np.linspace(0, np.pi / 2, 150)
    motion.plot(
        np.degrees(theta),
        data[0] * np.cos(theta) ** 2 + data[-1] * np.sin(theta) ** 2,
        color="#ae472e",
    )
    motion.scatter(np.degrees(measure.ANGLES), data, color="#ae472e", s=20)
    motion.axhline(0, color=".45", ls="--", lw=0.9)
    motion.text(87, -0.025, "inward", ha="right", fontsize=8, color=".35")
    motion.text(3, 0.024, "outward", fontsize=8, color=".35")
    motion.set(
        xlabel="Native exchange angle (degrees)",
        ylabel="Change in mean radius (hops)",
        title="Material response at time four",
        xticks=[0, 30, 60, 90],
    )
    for axis in (access, motion):
        axis.spines[["right", "top"]].set_visible(False)
        axis.tick_params(labelsize=8)
        axis.title.set_fontsize(10)
        axis.xaxis.label.set_size(9)
        axis.yaxis.label.set_size(9)
    fig.savefig(
        HERE / "recovery_flight.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)


if __name__ == "__main__":
    draw()
