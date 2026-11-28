"""Plot the measured-charge release and a withheld outer-galaxy rotation curve."""

import json

import matplotlib.pyplot as plt
import numpy as np

from . import pair_galaxies as galaxies


def draw():
    root = galaxies.HERE
    lab = galaxies.transfer.budget.local.lab.HERE
    release = json.loads((lab / "pair_probe.json").read_bytes())["decoded"]["64"]
    curves = json.loads((root / "pair_galaxies.json").read_bytes())["tests"]
    data = galaxies.transfer.budget.load_inputs()[0]["NGC2841"]
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.15, 2.85), layout="constrained")
    cuts = np.array([2, 3, 4, 6])
    left.axhline(0, color="0.65", lw=0.8)
    for k, offset, color, marker, label in (
        ("from_charge", -0.075, "#9f462e", "x", "From connected charge"),
        ("held", 0.075, "#285c70", "o", "From material readings"),
    ):
        bounds = np.array([release[f"{r}_undisturbed"][k] for r in cuts]) * 1e7
        y = bounds.mean(1)
        left.errorbar(
            cuts + offset,
            y,
            yerr=(bounds[:, 1] - bounds[:, 0]) / 2,
            fmt=marker,
            markersize=3.5,
            color=color,
            label=label,
            capsize=2,
            lw=1,
        )
    left.set(
        title="Separate probe at release",
        xlabel="Enclosing boundary (native hops)",
        ylabel=r"Change in inward acceleration ($10^{-7}$)",
        xticks=cuts,
    )
    left.legend(frameon=False, fontsize=7)
    r = np.array(data["radius_kpc"])
    right.errorbar(
        r,
        data["observed_kms"],
        yerr=2 * np.array(data["error_kms"]),
        fmt=".",
        ms=2,
        color="#202020",
        alpha=0.75,
        elinewidth=0.65,
        label="SPARC, twice quoted errors",
    )
    for tag, color, style, label in (
        ("atomic", "#285c70", "-", "Atomic gas pairs"),
        ("mapped_baryons", "#9f462e", "--", "Mapped baryon pairs"),
    ):
        right.plot(
            r,
            curves[tag]["galaxies"]["NGC2841"]["forecast_kms"],
            style,
            color=color,
            lw=1.2,
            label=label,
        )
    right.set(
        title="Withheld NGC 2841",
        xlabel="Orbital radius (kpc)",
        ylabel="Circular speed (km/s)",
    )
    right.legend(frameon=False, fontsize=6.7)
    path = root / "pair_comparison.pdf"
    fig.savefig(path, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    return path


if __name__ == "__main__":
    print(draw())
