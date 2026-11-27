"""Plot the central and anchor-shifted forecasts together with the three galaxy clocks."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(receipt, output, windows=None):
    d = json.loads(receipt.read_text())
    measured = d["inputs"]
    h = d["vacuum_H"]
    f = d["central_forecast_using_lowest_z"]
    low, high = d["forecasts_at_minus_plus_one_anchor_sigma"]
    if windows:
        window_result = json.loads(windows.read_text())
        f = window_result["background"]
        low, high = window_result["clock_anchor_variation"]
    blue, black, gold = "#136c70", "#273a53", "#aa5c17"
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.1, 3.7), layout="constrained")
    z = np.linspace(0, 1, 240)

    def expansion(z, branch):
        return np.sqrt(h * h + branch["matter_coefficient"] * (1 + z) ** 3)

    left.fill_between(
        z, expansion(z, low), expansion(z, high), color=blue, alpha=0.14, lw=0
    )
    left.plot(z, expansion(z, f), color=blue, label="Flat dust + assigned screen Λ")
    if windows:
        left.plot(
            measured["redshifts"][1:],
            window_result["reported_H_for_three_bins"][1:],
            ls="none",
            marker="D",
            markerfacecolor="none",
            markeredgecolor=blue,
            label="Background through fitted bins",
        )
    for i, ((zz, value), sig) in enumerate(
        zip(
            zip(measured["redshifts"], measured["H_observed"]),
            np.hypot(measured["stat_sigma"], measured["syst_sigma"]),
        )
    ):
        left.errorbar(
            zz,
            value,
            yerr=sig,
            marker="s" if i == 0 else "o",
            color=gold if i == 0 else black,
            linestyle="none",
            capsize=3,
            label=["Clock used for calibration", "Comparison clocks", None][i],
        )
    left.set(
        ylabel="Expansion rate (km/s/Mpc)",
        xlabel="Redshift",
        xlim=(0, 1),
        ylim=(45, 153),
    )
    left.set_title(
        "Expansion from galaxy spectra", loc="left", fontweight="bold", fontsize=12
    )
    left.legend(frameon=False, fontsize=8, loc="upper left")
    right.axhline(0, color="#888888", lw=0.8)
    z = np.linspace(0, 4.1, 400)

    def drift(branch):
        return (
            29979245800
            * 315576000
            / 3.0856775814913673e19
            * (branch["H0"] - expansion(z, branch) / (1 + z))
        )

    endpoints = np.stack([drift(low), drift(high)])
    right.fill_between(
        z, endpoints.min(axis=0), endpoints.max(axis=0), color=blue, alpha=0.14, lw=0
    )
    right.plot(z, drift(f), color=blue)
    for label, value in f["drift_cm_s_per_ten_julian_years"].items():
        right.plot(int(label), value, "o", color=blue, ms=4)
        right.annotate(
            f"{value:+.2f}",
            (int(label), value),
            xytext=(-4, 9),
            textcoords="offset points",
            ha="right",
            fontsize=9,
        )
    right.set(
        xlabel="Absorber redshift",
        ylabel="Equivalent velocity change in 10 years (cm/s)",
        xlim=(0, 4.15),
    )
    right.set_title(
        "Conditional spectral-line drift", loc="left", fontweight="bold", fontsize=12
    )
    for ax in (left, right):
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
    fig.savefig(
        output, dpi=180, bbox_inches="tight", metadata={"Software": "Matplotlib"}
    )
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--windows", type=Path, help="Use the retained-bin reconstruction"
    )
    args = parser.parse_args()
    plot(args.receipt, args.output, args.windows)
