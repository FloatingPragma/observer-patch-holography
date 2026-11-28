"""Plot the 13 DESI distances against the four low-redshift calibrations."""

import json

import matplotlib.pyplot as plt
import numpy as np

from .measure import HERE


def plot(result):
    z = np.asarray(result["z"])
    obs = np.asarray(result["observed"])
    sigma = np.asarray(result["sigma"])
    kind = np.asarray(result["kind"])
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(7.8, 5.8), layout="constrained")
    for ax, mask, label in (
        (axs[0], kind != "DH_over_rs", r"Transverse $D_M/r_d$ and low-$z$ $D_V/r_d$"),
        (axs[1], kind == "DH_over_rs", r"Radial $D_H/r_d$"),
    ):
        ax.axvspan(0.2, 1.0, color="#efefef")
        ax.axvline(1.0, color="#888888", lw=0.8, ls=":")
        ax.axhline(0, color="#434343", lw=0.7)
        selected = np.flatnonzero(mask)
        x = z[selected]
        ax.fill_between(
            x, -np.ones_like(x), np.ones_like(x), color="#e3eff6", alpha=0.9
        )
        for key, label_model, color, mark, dx in (
            ("20", "20 owners", "#aa582f", "x", -0.014),
            ("320", "320 owners", "#835493", "+", 0.014),
            ("standard_control", "Standard control", "#292c36", "o", -0.005),
            ("80", "80 owners", "#087f75", "s", 0.005),
        ):
            prediction = np.asarray(
                result["stocks"][key]["runs"]["0.02237"]["low_z_H0"]["predicted"]
            )
            ax.plot(
                x + dx,
                ((prediction - obs) / sigma)[selected],
                linestyle="none",
                marker=mark,
                markersize=6,
                mfc="none",
                mew=1.2,
                color=color,
                label=label_model,
            )
        ax.set_title(label, fontsize=10.5, loc="left")
        ax.set_ylabel("Residual / measured error")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
    axs[0].legend(fontsize=8.5, ncol=4, loc="upper left")
    axs[-1].set_xlabel("Effective redshift")
    axs[-1].set_xlim(0.19, 2.44)
    fig.suptitle("DESI DR2 distances after calibrating H₀ only at z ≤ 1", fontsize=13)
    fig.savefig(HERE / "distances.pdf")
    plt.close(fig)


if __name__ == "__main__":
    plot(json.loads((HERE / "results.json").read_text()))
