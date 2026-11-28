"""Plot recorded enclosed fractions and the closest certified speed-box example."""

from fractions import Fraction as F
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cosmology.rar_deep_regime.brightness_pairs import midpoint_halo as sky
from . import midpoint_bulk as bulk


def run():
    source = json.loads((bulk.HERE / "midpoint_bulk.json").read_bytes())["levels"][-1]
    panel = [
        r
        for r in source["observations"]
        if r["center"] is None and r["address_law"] == "solid_angle"
    ]
    data = json.loads((sky.HERE / "midpoint_halo.json").read_bytes())["galaxies"]
    gal = "UGC09133"
    row = data[gal]
    relaxed = row["additional_speed_halfwidth_kms"]
    params = list(map(lambda v: float(F(v)), relaxed["high_certificate"]["parameters"]))
    radii = np.array([float(F(s)) for s in row["radii_kpc"]])
    raw = (sky.sample.measure.DATA / "table2.dat").read_text().splitlines()
    r, v, e, pred = [], [], [], []
    for line in row["lines"]:
        t = raw[line - 1]
        R = float(t[19:25])
        r.append(R)
        v.append(float(t[26:32]))
        e.append(2 * float(t[33:38]))
        bary = np.array([float(t[a:b]) for a, b in ((39, 45), (46, 52), (53, 59))])
        squared = np.dot(params[:3], bary * abs(bary)) + np.dot(
            params[3:], np.minimum(R, radii**2 / R)
        )
        pred.append(np.sqrt(squared))
    with plt.rc_context(
        {
            "font.size": 9,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "pdf.fonttype": 42,
        }
    ):
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.05), layout="constrained")
        left, right = axes
        grid = np.linspace(0, 1, 120)
        left.plot(grid, grid**2, color="#185460", label="Round-screen prediction")
        x = [s["radius"] for s in panel]
        bands = np.array([s["separation_read"]["enclosed_fraction"] for s in panel])
        left.errorbar(
            x,
            bands.mean(axis=1),
            yerr=(bands[:, 1] - bands[:, 0]) / 2,
            fmt="o",
            ms=5,
            capsize=3,
            color="#b75229",
            label="Endpoint-separation counts",
        )
        left.set(
            xlim=(0, 1),
            ylim=(0, 1.05),
            xlabel="Radius / screen radius",
            ylabel="Fraction of pair midpoints enclosed",
            title="5,120-station angular survey",
        )
        left.legend(frameon=False, fontsize=8, loc="upper left")
        right.errorbar(
            r,
            v,
            yerr=e,
            fmt=".",
            color="#313d4d",
            capsize=2,
            label="SPARC, twice the speed error",
        )
        right.plot(r, pred, color="#b75229", label="Mixture with expanded speed errors")
        i = [row["lines"].index(n) for n in row["mixture_witness"]["lines"]]
        right.scatter(
            np.array(r)[i],
            np.array(v)[i],
            s=48,
            facecolors="none",
            edgecolors="#185460",
            label="Three exclusion rings",
            zorder=5,
        )
        right.set(
            xlabel="Catalogue radius (kpc)",
            ylabel="Circular speed (km/s)",
            title="UGC 09133 outer rotation curve",
        )
        right.legend(frameon=False, fontsize=7)
        out = Path(__file__).with_name("midpoint_bulk_comparison.pdf")
        fig.savefig(out, metadata={"CreationDate": None, "ModDate": None})
        plt.close(fig)
    return out


if __name__ == "__main__":
    print(run())
