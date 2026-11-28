"""Plot the measured shell tide and the pre-existing DDO 154 rotation data."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    here = Path(__file__).resolve().parent
    data = json.loads((here / "tidal_rotation.json").read_bytes())["distributed_disks"]
    tab = json.loads((here / "source_profile.json").read_bytes())["resolved_maps"][
        "inputs"
    ]["DDO154"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.1, 3.1), constrained_layout=True)
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
        ax.axhline(0, color="#8a8a8a", lw=0.7)
        ax.tick_params(labelsize=8)
    for name, color, label, marker in [
        ("source", "#155776", "Canonical repair", "o"),
        ("positive", "#984920", "Matched branching excess", "s"),
    ]:
        row = data[name]
        left.plot(
            row["shells"],
            1000 * np.array(row["native_inward_tide"]),
            color=color,
            label=label,
            marker=marker,
            ms=3.2,
            lw=1.4,
        )
    left.set(
        xlabel="Radius from station 219 (hops)",
        ylabel=r"Inward tide ($10^{-3}$ per clock unit$^2$)",
        xticks=[1, 4, 7, 10, 13, 15],
    )
    left.legend(frameon=False, fontsize=8)
    r, v, e, b = [
        np.asarray(tab[k])
        for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
    ]
    p = np.array(
        data["source"]["forecasts"]["6/mapped_baryons/1"]["galaxies"]["DDO154"][
            "predicted_signed_v2_kms2"
        ]
    )
    y = v**2 - b**2
    right.errorbar(
        r,
        y,
        yerr=[v**2 - np.maximum(0, v - 2 * e) ** 2, (v + 2 * e) ** 2 - v**2],
        color="#333333",
        fmt="o",
        ms=3.4,
        lw=0.9,
        capsize=1.8,
        label="SPARC excess (two errors)",
    )
    right.plot(r, p - b**2, color="#155776", lw=1.5, label="Distributed repair profile")
    right.set(
        xlabel="DDO 154 radius (kpc)",
        ylabel=r"$v^2-v_{\rm baryon}^2$ (km$^2$ s$^{-2}$)",
    )
    right.legend(frameon=False, fontsize=8, loc="center right")
    fig.savefig(
        here / "tidal_profile.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
