"""Show galaxy exceedances and the outer-orbit thresholds in NGC 3741."""

import json

import matplotlib.pyplot as plt
import numpy as np

from .measure import HERE


def plot(r):
    fig, (a, b) = plt.subplots(2, 1, figsize=(7.1, 6.5), layout="constrained")
    styles = [
        ("nominal", "Catalogue center", "#232f55", "-"),
        ("one_error", "One-error favorable box", "#267b8a", "--"),
        ("two_errors", "Two-error favorable box", "#b15a28", "-."),
        ("extreme", "Extreme box", "#866992", ":"),
    ]
    x = np.linspace(0, 17, 1801)
    for box, label, color, ls in styles:
        q = np.array([row["bounds"][box]["q_required_min"] for row in r["galaxies"]])
        # Reverse cumulative counts of galaxies above an arbitrary mass budget.
        a.step(
            x,
            np.sum(q[:, None] > x, axis=0),
            where="post",
            label=label,
            color=color,
            ls=ls,
            lw=1.6,
        )
    for name, d in r["source_caps"].items():
        a.axvline(d["cap"], color=".55", alpha=0.6, lw=0.8)
        a.text(
            d["cap"] + 0.12,
            98,
            f"{name} seams",
            fontsize=8,
            rotation=90,
            color=".25",
            va="center",
        )
    a.legend(fontsize=8.5)
    a.set(
        xlim=(0, 17),
        ylim=(-2, 120),
        xlabel="Available spherical stock / galaxy baryon mass",
        ylabel="Galaxies above budget",
        title="118 individual outer rotation measurements",
    )
    curve = r["worst_example_rotation_curve"]
    for box, label, color, ls in styles:
        b.plot(
            [p["r_kpc"] for p in curve],
            [p["q_required"][box] for p in curve],
            label=label,
            ls=ls,
            color=color,
            lw=1.6,
        )
    for name, cap in r["source_caps"].items():
        b.axhline(cap["cap"], color=".55", lw=0.8)
        b.text(0.12, cap["cap"] + 0.2, f"{name} seams", fontsize=8, color=".25")
    extra = r["NGC3741_controls"]["removed_pressure_correction"]
    b.scatter(
        [p["catalogue_radius_kpc"] for p in extra],
        [p["q_required"]["two_errors"] for p in extra],
        label="Two-error box, full drift correction removed",
        marker="x",
        color="#8c271d",
        s=44,
        zorder=4,
    )
    b.legend(fontsize=8, loc="upper left", ncol=2)
    b.set(
        xlim=(0, 7.2),
        ylim=(-1, 18),
        xlabel="Published radius (kpc)",
        ylabel="Required stock / baryon mass",
        title="NGC 3741 at the measured radii",
    )
    for ax in (a, b):
        ax.grid(axis="y", alpha=0.14)
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(HERE / "orbits.pdf")
    plt.close(fig)


if __name__ == "__main__":
    plot(json.loads((HERE / "results.json").read_text()))
