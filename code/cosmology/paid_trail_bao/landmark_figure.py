"""Plot two-probe endpoint predictions and the disk predicted from the same repair."""

from fractions import Fraction as Q
import json
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import tidal_scales as scales


def main():
    here = scales.HERE
    native = json.loads((scales.budget.LAB / "landmark_arrival.json").read_bytes())[
        "panels"
    ]["25/36"]
    models = json.loads((here / "landmark_abundance.json").read_bytes())
    data = json.loads((here / "source_profile.json").read_bytes())["resolved_maps"][
        "inputs"
    ]["DDO154"]
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(8.35, 3.4), constrained_layout=True
    )
    colors = ["#176880", "#af552a"]
    old = native["fixed_stations"]["104"]
    common = 1e4 * float(Q(native["matched_partial_fixed_stations"]["104"]))
    for j in (0, 1):
        left.scatter(
            [0 + (j - 0.5) * 0.12, 1 + (j - 0.5) * 0.12],
            [1e4 * float(Q(old[str(1 + j)])), common],
            marker=["o", "D"][j],
            s=34,
            c=colors[j],
            label=["Single", "Double"][j],
            zorder=4,
        )
    lo, hi = 1e4 * np.array(native["held_transfer"]["matched_common_change"])
    left.vlines(1, lo, hi, color="#62676d", linewidth=2, label="From measured source")
    left.hlines([lo, hi], 0.94, 1.06, color="#62676d", linewidth=1)
    left.axhline(0, color="#888", linewidth=0.7)
    left.set(
        xlim=(-0.42, 1.48),
        ylim=(lo - 1.5, 2.2),
        xticks=[0, 1],
        xticklabels=["Canonical\nrepair", "Energy-matched\nentropy repair"],
        ylabel="Change in arrivals at station 36\nper 10,000 launches",
        title="Two freely traveling probes",
    )
    left.legend(loc="lower left", fontsize=7, frameon=False)
    radius = np.array(data["radius_kpc"])
    baryons = np.array(data["baryons_kms"])
    right.errorbar(
        radius,
        data["observed_kms"],
        yerr=2 * np.asarray(data["error_kms"]),
        fmt=".",
        ms=3.4,
        lw=0.6,
        color="#33383b",
        label="SPARC, twice quoted errors",
    )
    right.plot(radius, baryons, ":", color="#65686b", lw=1.3, label="Baryons")
    for key, label, color in zip(
        ("one_proposal", "renewed_until_accepted"),
        ("One proposal", "Renewed proposals"),
        colors,
        strict=True,
    ):
        v = models[key]["mapped_baryons_first_capacity"]["disk/spherical"]["DDO154"]
        right.plot(
            radius,
            np.sqrt(baryons**2 + v["predicted_excess_v2"]),
            color=color,
            lw=1.6,
            label=label,
        )
    right.set(
        xlabel="Radius (kpc)",
        ylabel="Circular speed (km/s)",
        title="DDO 154, imaged source population",
    )
    right.legend(loc="center right", bbox_to_anchor=(0.99, 0.51), fontsize=7, frameon=False)
    for ax in [left, right]:
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
        ax.title.set_size(10)
        ax.yaxis.label.set_size(9)
    fig.savefig(
        here / "landmark_source.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
