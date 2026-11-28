"""Plot the independent material prediction and the galaxy's required exterior stock."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from . import exterior_marks as marks


def draw():
    base = marks.HERE
    lab = marks.budget.local.lab.HERE
    theory = json.loads((lab / "poisson_flight.json").read_bytes())["forecast"][
        "written"
    ]
    obs = json.loads((lab / "poisson_precision.json").read_bytes())["decoded"][
        "written"
    ]
    inputs, charts, _, masses = marks.budget.load_inputs()
    name = "NGC2841"
    data = marks.transport.mass_intervals(inputs[name], masses[name])
    r, lo, hi = np.array(data, dtype=float).T
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.15, 2.85), layout="constrained")
    x = np.array([row["radius"] for row in obs])
    y = np.array([row["from_position"] for row in obs]) * 1e8
    intervals = np.array([row["terminal_interval"] for row in obs]) * 1e8
    left.axhline(0, color="0.65", lw=0.8)
    left.errorbar(
        x - 0.05,
        y,
        yerr=np.vstack((y - intervals[:, 0], intervals[:, 1] - y)),
        fmt="o",
        ms=4,
        capsize=2.5,
        color="#24546b",
        label="Terminal counts, joint 99% bounds",
    )
    left.scatter(
        x + 0.09,
        np.array([row["linear_displacement"] for row in theory]) * 1e8,
        marker="x",
        s=32,
        color="#ac442b",
        label="Ordinary-calibrated prediction",
        zorder=3,
    )
    left.set(
        title="Held written contact",
        xlabel="Enclosing boundary (native hops)",
        ylabel=r"Extra enclosed probability ($10^{-8}$)",
        xticks=x,
        ylim=(-0.2, 7.0),
    )
    left.legend(loc="lower left", fontsize=6.7, frameon=False)
    edge = charts[name][0][-1]
    right.axvspan(0, edge, color="0.91", label="Mapped baryonic origins")
    right.vlines(
        r, lo, hi, color="#24546b", lw=1.0, label="SPARC two-error mass intervals"
    )
    for i in (27, 49):
        right.vlines(r[i], lo[i], hi[i], color="#ac442b", lw=2.2)
    right.set(
        title="NGC 2841",
        xlabel="Orbital radius (kpc)",
        ylabel=r"Required extra mass inside radius / $M_b$",
        xlim=(0, 68),
        ylim=(0, 1.06 * max(hi)),
    )
    right.legend(loc="upper left", fontsize=6.7, frameon=False)
    path = Path(base) / "poisson_comparison.pdf"
    fig.savefig(path, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    return path


if __name__ == "__main__":
    print(draw())
