"""The direct source reading and its dependence on regional correlations."""

from fractions import Fraction as F
import json

import matplotlib.pyplot as plt
import numpy as np

from . import fluctuation_balls as lab


def run():
    local = json.loads((lab.HERE / "fluctuation_return.json").read_bytes())[
        "closed_source"
    ]
    reg = json.loads((lab.HERE / "fluctuation_balls.json").read_bytes())
    _, w, *_ = lab.energy.pair.setup()
    board, homes, seams, *_ = lab.energy.returned(w)
    _, dist, _ = lab.radar.chart(board, seams)
    r = np.arange(reg["diameter"] + 1)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(8, 3.3), layout="constrained")
    colors = ("#255a74", "#ad4a26")
    for idx, (name, x) in enumerate(local.items()):
        offset = (-0.09, 0.09)[idx]
        for k in (1, 2):
            c = x["original_middle_write"]["clocks"][str(k)]
            value = 100 * c["transits"] / c["trials"]
            low, high = 100 * np.array(c["interval"])
            ax.errorbar(
                k + offset,
                value,
                yerr=[[value - low], [high - value]],
                fmt="o",
                ms=5,
                color=colors[idx],
                capsize=3,
                label=name.capitalize() if k == 1 else None,
            )
            ax.plot(
                [k + offset - 0.055, k + offset + 0.055],
                [100 * c["pooled"]] * 2,
                color="black",
                lw=1,
            )
        m = ((dist[homes[1], homes, None] <= r) * (1 << np.arange(4))[:, None]).sum(0)
        var = np.array([float(F(v)) for v in reg["arms"][name]["variances"]])
        bx.step(
            r,
            np.sqrt(var[m]).cumsum(),
            where="post",
            color=colors[idx],
            lw=2,
            label=f"Middle · {name}",
        )
    bx.plot(
        r, np.zeros(len(r)), "--", color="#353535", lw=1.5, label="Equidistant · both"
    )
    ax.plot([], [], color="black", label="Finite prediction", lw=1)
    ax.set(
        title="Material receiver",
        xlabel="Attempts per appointment",
        ylabel="Arrivals before return (%)",
        xticks=[1, 2],
        xlim=(0.63, 2.37),
        ylim=(13, 30),
    )
    bx.set(
        title="Nested material apertures",
        xlabel="Radius (ordinary echo hops)",
        ylabel="Cumulative RMS (load units)",
        xlim=(0, r[-1]),
        ylim=(-0.3, 8.6),
    )
    for axis in (ax, bx):
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.18)
        axis.legend(frameon=False, fontsize=8, loc="upper left")
    fig.savefig(
        lab.HERE / "fluctuation_return.pdf",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


if __name__ == "__main__":
    run()
