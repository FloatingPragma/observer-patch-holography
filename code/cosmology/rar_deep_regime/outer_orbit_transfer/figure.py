"""Plot the outer-orbit shape constraint and lens-calibrated acceleration intervals."""

import json

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from .. import lensing_profile
from .linear_profile import orbit_interval
from .measure import HERE, orbit


def plot(path):
    rows = {
        row["name"]: row
        for row in json.loads((HERE / "inputs.json").read_text())["galaxies"]
    }
    panels = json.loads((HERE / "results.json").read_text())["panels"]
    lens = json.loads((HERE / "linear.json").read_text())[
        "lens_two_statistical_error_intervals_a0_m_s2"
    ]
    fig, ((ax, bx), shear_axes) = plt.subplots(
        2,
        2,
        figsize=(9, 6.5),
        layout="constrained",
        gridspec_kw={"width_ratios": [1.22, 1]},
    )
    colors = {"NGC3626": "#ae3a28", "NGC4262": "#1e658a"}
    for name in panels["primary"]["names"]:
        g, f = orbit(rows[name])
        x = np.sqrt(1e-11 / g)
        lo, hi = (orbit(rows[name], k)[1] for k in (-1, 1))
        ax.errorbar(
            x,
            f,
            yerr=[[f - lo], [hi - f]],
            fmt="o",
            ms=4,
            capsize=2.5,
            color=colors.get(name, "#808486"),
        )
        if name in colors:
            ax.annotate(
                name.replace("NGC", "NGC "),
                (x, f),
                xytext=(6, 5) if name == "NGC3626" else (6, -1),
                textcoords="offset points",
                fontsize=8,
                color=colors[name],
            )
    gx, fu = orbit(rows["NGC3626"], 1)
    xs = np.linspace(0.3, 1.12, 100)
    ax.plot(
        xs, fu * np.maximum(1, xs / np.sqrt(1e-11 / gx)), color=colors["NGC3626"], lw=1
    )
    ax.set(
        xlabel=r"Scaled radius $\sqrt{10^{-11}\,\mathrm{m}\,\mathrm{s}^{-2}/g_b}$",
        ylabel=r"Enclosed halo / total baryonic mass",
        title="Nine exterior orbits",
    )
    items = ["KiDS-bright", "GAMA", *panels["secondary_deep"]["names"], "NGC3626"]
    for i, name in enumerate(reversed(items)):
        if name in lens:
            lower, upper = lens[name]
            color = "#604474"
        else:
            lower, upper = orbit_interval(rows[name], 1.0)
            lower, upper = lower**2, max(0, upper) ** 2
            color = colors.get(name, "#28715e")
        bx.plot(
            np.array([lower, upper]) / 1e-10,
            [i, i],
            lw=3,
            solid_capstyle="round",
            color=color,
        )
    bx.axvspan(
        max(b[0] for b in lens.values()) / 1e-10,
        min(b[1] for b in lens.values()) / 1e-10,
        color="#604474",
        alpha=0.13,
    )
    bx.set(
        yticks=list(range(len(items))),
        yticklabels=[n.replace("NGC", "NGC ") for n in reversed(items)],
        xlabel=r"Acceleration scale $a_0$ ($10^{-10}$ m s$^{-2}$)",
        title="Linear exterior profile",
    )
    surveys = json.loads(lensing_profile.RECEIPT.read_text())["profiles"]
    scale = json.loads((HERE / "linear.json").read_text())["panels"]["secondary_deep"][
        "distance_and_speed"
    ]["a0_m_s2"]
    for axis, name in zip(shear_axes, ("KiDS-bright", "GAMA")):
        profile = surveys[name]
        low = np.asarray(profile["lower_bin_edge_m_s2"])
        mask = low >= 1e-13
        center = np.asarray(profile["g_bar_m_s2"])[mask]
        low = low[mask]
        high = 2 * center - low
        observed = np.asarray(profile["tangential_esd"])[mask]
        sigma = np.sqrt(np.diag(profile["covariance"]))[mask]

        def shear(g):
            return (4 / np.pi * g + np.sqrt(g * scale)) / lensing_profile.FOUR_G

        axis.fill_between(
            center / 1e-12,
            shear(low),
            shear(high),
            color="#604474",
            alpha=0.15,
            label="Within-bin range",
        )
        axis.plot(
            center / 1e-12,
            shear(center),
            color="#604474",
            label="Profile at bin centre",
        )
        axis.errorbar(
            center / 1e-12,
            observed,
            yerr=sigma,
            color="#282c32",
            fmt="o",
            ms=4,
            capsize=2,
            label="Measured",
        )
        axis.set(
            xscale="log",
            title=name,
            xlabel=r"Baryonic acceleration ($10^{-12}$ m s$^{-2}$)",
        )
        axis.set_xticks([0.2, 0.5, 1, 2, 4], labels=["0.2", "0.5", "1", "2", "4"])
        axis.xaxis.set_minor_formatter(NullFormatter())
    shear_axes[0].set_ylabel(r"Excess surface density ($M_\odot$ pc$^{-2}$)")
    shear_axes[0].legend(fontsize=8, frameon=False)
    for axis in (ax, bx, *shear_axes):
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="x", alpha=0.16)
        axis.set_axisbelow(True)
    fig.savefig(path, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    plot(HERE / "comparison.pdf")
