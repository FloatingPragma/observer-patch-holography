"""Plot the measured source and its differential clock prediction in DDO 154."""

import json

import matplotlib.pyplot as plt
import numpy as np

from .closed_inventory import HERE, LAB


def run():
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.2, 3.35), layout="constrained")
    source = json.loads((LAB / "sector_return_clocks.json").read_text())[
        "same_home_coincidences"
    ]
    x = np.arange(4)
    values = np.array(source["fresh_estimate_and_standard_error"]) * 1e6
    a.errorbar(
        x - 0.06,
        values[:, 0],
        yerr=2 * values[:, 1],
        fmt="o",
        color="#19715c",
        label="Model records ±2 SE",
        capsize=3,
    )
    a.scatter(
        x + 0.06,
        np.array(source["prediction"]) * 1e6,
        color="#1f344c",
        marker="x",
        label="Finite prediction",
        zorder=4,
    )
    a.axhline(0, color="#50555a", lw=0.7)
    a.set(
        xticks=x,
        xticklabels=source["homes"],
        xlabel="Home address",
        ylabel="Change in joint occupancy (ppm)",
        title="Two reports from one register",
    )
    a.legend(frameon=False, fontsize=8, loc="lower right")
    s = json.loads((HERE / "sector_return_lapse.json").read_text())["profiles"][
        "DDO154"
    ]
    sky = json.loads((HERE / "sector_return_redshift.json").read_text())["ten_imaged"][
        "galaxies"
    ]["DDO154"]
    use = [i for i, z in enumerate(s["clocks"]) if z is not None]
    r = np.array(s["radius_kpc"])[use]
    raw = np.array([s["clocks"][i]["extra_log_frequency"] for i in use])
    q = raw - raw[0, ::-1]
    lo = np.array(sky["minimum_cumulative_log_rate_difference"])[use]
    lo -= lo[0]
    assert np.all(q[1:, 0] > 0) and np.all(lo[1:] > 0)
    b.plot(
        r[1:], lo[1:], color="#1f344c", marker="o", ms=3, label="Minimum from rotation"
    )
    b.fill_between(r[1:], q[1:, 0], q[1:, 1], color="#19715c", alpha=0.25)
    b.plot(
        r[1:],
        q[1:].mean(1),
        color="#19715c",
        marker="s",
        ms=3,
        label="Departure-clock prediction",
    )
    b.set(
        yscale="log",
        xlabel="Radius (kpc)",
        ylabel="Extra log frequency ratio",
        title="DDO 154 relative to 1.48 kpc",
    )
    b.legend(frameon=False, fontsize=8)
    for ax in (a, b):
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(axis="y", alpha=0.18)
    fig.savefig(
        HERE / "sector_return_clocks.pdf",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


if __name__ == "__main__":
    run()
