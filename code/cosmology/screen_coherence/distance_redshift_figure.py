"""Plot observer zero's paired clock and distance readings, and a separate measured reference."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np

PANTHEON_COMMIT = "7fc680548d1ea9fe6ee983a89d8b5635bb5784a8"
PANTHEON_SHA = "1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8"
PANTHEON_URL = f"https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/{PANTHEON_COMMIT}/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
INK = "#172e3b"
MUTED = "#435968"
BLUE = "#12577d"
TEAL = "#09766c"


def style_axes(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#d8e1e3", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    for edge in ("left", "bottom"):
        ax.spines[edge].set_color("#9fb0b5")
    ax.tick_params(colors=MUTED, labelsize=9)


def diagram(receipt, out):
    data = json.loads(receipt.read_text())
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.25), sharex=True, sharey="row")
    fig.subplots_adjust(
        left=0.12, right=0.97, bottom=0.18, top=0.8, wspace=0.13, hspace=0.20
    )
    fig.text(
        0.075,
        0.943,
        "Observer 0’s distance–redshift diagram",
        fontsize=23,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.077,
        0.897,
        "Level 13  ·  Material-signal replay  ·  274 enrolled beacons on 28 remote observers",
        fontsize=11.3,
        color=MUTED,
    )
    for col, (epoch, label, color) in enumerate(
        (("0", "Query at tick 0", BLUE), ("96", "Settled · query at tick 96", TEAL))
    ):
        rows = [
            r
            for r in data["panels"][epoch]
            if r["z_pulse"] is not None and r["radar_distance"] is not None
        ]
        ax, flash = axes[:, col]
        ax.set_title(
            f"{label}  ·  {len(rows)} paired",
            fontsize=12,
            weight="bold",
            color=INK,
            loc="left",
            pad=11,
        )
        grouped = Counter((r["z_pulse"], r["radar_distance"]) for r in rows)
        for (z, r), n in grouped.items():
            ax.scatter(
                z,
                r,
                s=14 * n**0.78,
                c=color,
                alpha=0.82,
                lw=0.6,
                edgecolors="white",
                zorder=3,
            )
        for r in rows:
            d = r["flash"]["equivalent_distances"][-1]
            if d is not None:
                flash.scatter(r["z_pulse"], d, s=13, color=color, alpha=0.47, lw=0)
        for item in (ax, flash):
            style_axes(item)
            item.axvline(0, color="#a2b1b5", lw=0.8, zorder=0)
            item.set_xlim(-0.12, 2.40)
            item.set_xticks(
                [0, 0.5, 1, 1.5, 2, 2.25], ["0", "0.5", "1", "1.5", "2", "2.25"]
            )
        ax.set(ylim=(0.5, 16.6), yticks=[2, 5, 8, 11, 14])
        flash.set_yscale("log")
        flash.set_ylim(1.15, 180)
        flash.yaxis.set_major_locator(FixedLocator([2, 5, 10, 30, 100]))
        flash.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        flash.set_xlabel(
            "One-way pulse redshift  $z_{\\rm pulse}$", fontsize=11, color=INK
        )
    axes[0, 0].set_ylabel("Radar distance\n(relay units)", color=INK, fontsize=11)
    axes[1, 0].set_ylabel(
        "Flash-count distance\n(local reference = 1)", color=INK, fontsize=11
    )
    for n in (1, 5, 15):
        axes[0, 1].scatter([], [], s=14 * n**0.78, c=TEAL, label=str(n), alpha=0.82)
    axes[0, 1].legend(
        title="Beacons per radar point",
        loc="upper right",
        ncol=3,
        fontsize=8,
        title_fontsize=9,
        frameon=False,
        columnspacing=0.9,
        handletextpad=0.2,
    )
    fig.text(
        0.078,
        0.082,
        "Redshift comes from two incoming clock marks. Radar distance comes from a separate round trip.",
        color=MUTED,
        fontsize=9.8,
    )
    fig.text(
        0.078,
        0.05,
        "Lower panels use a 32-tick flash-count gate. The two distance scales have different calibrations.",
        color=MUTED,
        fontsize=9.8,
    )
    fig.savefig(out, dpi=185, facecolor="white")
    plt.close(fig)


def observed(path, out):
    raw = (
        path.read_bytes()
        if path
        else urllib.request.urlopen(PANTHEON_URL, timeout=30).read()
    )
    if hashlib.sha256(raw).hexdigest() != PANTHEON_SHA:
        raise ValueError("Unexpected Pantheon+SH0ES public release")
    from io import BytesIO

    rows = np.genfromtxt(BytesIO(raw), names=True, usecols=("zHD", "MU_SH0ES"))
    z, dist = rows["zHD"], 10 ** ((rows["MU_SH0ES"] - 25) / 5)
    if not np.all(np.isfinite(dist) & (dist > 0) & np.isfinite(z) & (z > 0)):
        raise ValueError("Invalid observational inputs")
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    fig.subplots_adjust(left=0.145, right=0.96, bottom=0.21, top=0.81)
    fig.text(
        0.075,
        0.92,
        "Observed supernova distance–redshift diagram",
        fontsize=17,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.076,
        0.858,
        f"Pantheon+ and SH0ES  ·  {len(z):,} light-curve entries  ·  All released rows",
        fontsize=10.8,
        color=MUTED,
    )
    ax.scatter(z, dist, s=10, color="#953e60", alpha=0.53, lw=0)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="Hubble-diagram redshift  $z_{\\rm HD}$",
        ylabel="Luminosity distance (Mpc)",
    )
    style_axes(ax)
    fig.text(
        0.075,
        0.063,
        "Distances use the released SH0ES-calibrated distance moduli. Multiple observations of a supernova are retained.",
        fontsize=8.9,
        color=MUTED,
    )
    fig.savefig(out, dpi=180, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--observation-table", type=Path)
    args = parser.parse_args()
    diagram(args.receipt, args.receipt.with_suffix(".png"))
    observed(
        args.observation_table, args.receipt.with_name("pantheon_distance_redshift.png")
    )
