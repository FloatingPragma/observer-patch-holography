"""Plot the stationary directional response and two measured radial mass profiles."""

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def receipt(name):
    return json.loads((HERE / (name + ".json")).read_bytes())


def render():
    full = receipt("full_board_polarization")["outcomes"]
    predicted = receipt("polarization_predictions")["position_predictions"]
    tracking = receipt("material_polarization")
    compact = receipt("local_trapping")["surveys"][0]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(8.2, 2.8), layout="constrained")
    for y, (tag, color) in enumerate([("forward", "#156552"), ("reversed", "#94532e")]):
        data = full[tag]
        intervals = np.array(data["held_flight_enclosed"][-1]) * 1e6
        prediction = np.array(predicted[tag][-1]) * 1e6
        for observed, bounds in [(False, prediction), (True, intervals)]:
            ax.errorbar(
                bounds.mean(axis=1),
                np.arange(3) + 0.26 * (y - 0.5) + 0.10 * (observed - 0.5),
                xerr=(bounds[:, 1] - bounds[:, 0]) / 2,
                fmt="o" if observed else "s",
                ms=3.8,
                mfc=color if observed else "white",
                capsize=2,
                c=color,
                label=("Dispatch " + ("A" if y == 0 else "B")) if observed else None,
            )
    ax.axvline(0, color=".45", lw=0.7)
    ax.set(
        yticks=range(3),
        yticklabels=["Probe 1", "Probe 2", "Clock"],
        xlabel="Change in signed reading (millionths)",
        title="After the record is erased",
        ylim=(-0.4, 2.7),
    )
    ax.legend(frameon=False, fontsize=8, ncols=2, loc="upper center")
    base = next(
        v for v in tracking["full_boards"] if v["cells"] == 20 and v["stop_bits"] == 5
    )
    mass = np.array([s["enclosed"] for s in base["shells"]] + [0])
    hidden = np.cumsum(compact["shell_t6"])
    for series, label, color, marker in [
        (mass, "All ordinary material", "#94532e", "o"),
        (hidden, "Trapped modes", "#156552", "s"),
    ]:
        bx.plot(
            range(len(series)),
            series / max(series),
            marker=marker,
            ms=3.5,
            lw=1.4,
            color=color,
            label=label,
        )
    bx.set(
        xticks=range(6),
        xlabel="Outermost included shell (hops)",
        ylabel="Enclosed change / own maximum",
        title="Early response around dispatch",
        ylim=(-0.04, 1.08),
    )
    bx.legend(frameon=False, fontsize=8, loc="center right")
    for axis in (ax, bx):
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=9)
        axis.title.set_fontsize(11)
    fig.savefig(
        HERE / "ordinary_polarization.pdf",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


if __name__ == "__main__":
    render()
