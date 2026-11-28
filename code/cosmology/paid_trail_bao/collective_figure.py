"""Plot the collective source threshold and the unretuned ordinary-equivalent disk."""

from fractions import Fraction as Q
import json
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from . import tidal_scales as scales


def main():
    here = scales.HERE
    lab = json.loads((scales.budget.LAB / "collective_relay_lab.json").read_bytes())
    mass = json.loads((here / "pooled_ordinary.json").read_bytes())[
        "mapped_baryons_first_capacity"
    ]
    sky = json.loads((here / "source_profile.json").read_bytes())["resolved_maps"][
        "inputs"
    ]["DDO154"]
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(8.2, 3.05), constrained_layout=True
    )
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    x = np.arange(5)
    y = [-1000 * float(Q(lab["dose"][f"{k}/0"]["change"])) for k in x]
    left.axhline(0, color="#888", lw=0.8)
    left.plot(x, y, "o-", color="#155776", ms=4, label="Four route sources")
    lo, hi = -1000 * np.array(lab["held_change_interval"])[::-1]
    left.errorbar(
        [4.13],
        [(lo + hi) / 2],
        yerr=[(hi - lo) / 2],
        fmt="D",
        ms=3.5,
        capsize=3,
        color="#984920",
        label="All eight, held route read",
    )
    left.set(
        xlabel="Number of repaired sources on the route",
        ylabel=r"Inward change ($10^{-3}$ per clock unit$^2$)",
        xticks=x,
        ylim=(-0.10, 2.65),
    )
    left.legend(frameon=False, fontsize=8, loc="upper left")
    r, v, e, b = [
        np.array(sky[k])
        for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
    ]
    right.errorbar(
        r,
        v,
        yerr=2 * e,
        color="#333",
        fmt="o",
        ms=3,
        lw=0.8,
        capsize=2,
        label="SPARC (two errors)",
    )
    for pool, color, ls in (("disk", "#155776", "-"), ("tile", "#984920", "--")):
        y = np.array(mass[f"{pool}/spherical"]["DDO154"]["predicted_excess_v2"])
        right.plot(
            r,
            np.sqrt(b * b + y),
            ls,
            color=color,
            lw=1.3,
            label=f"{pool.capitalize()} pool",
        )
    right.set(xlabel="DDO 154 radius (kpc)", ylabel="Circular speed (km/s)")
    right.legend(frameon=False, fontsize=8)
    fig.savefig(
        here / "collective_source.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
