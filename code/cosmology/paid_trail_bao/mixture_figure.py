"""Distance likelihood and the galaxy setting the resolved two-error limit."""

import json

import matplotlib.pyplot as plt
import numpy as np

from . import rings


def plot():
    root = rings.mixture.bao.HERE
    cosmic = json.loads((root / "mixture.json").read_bytes())
    radial = json.loads((root / "rings.json").read_bytes())
    name = "F571-8"
    detail = radial["two_error_followup"][name]
    threshold = detail["distance_fit_at_feasible_cap"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.5, 3.45), layout="constrained")
    best = cosmic["best"]
    qlo, qhi = [fit["q"] for fit in cosmic["profile_95"]]
    fits = sorted(
        [*cosmic["profile_grid"], *cosmic["profile_95"], threshold],
        key=lambda f: f["q"],
    )
    a.axvspan(qlo, qhi, color="#cee6df", label="95% cosmic interval")
    a.plot(
        [f["q"] for f in fits],
        [f["chi2"] - best["chi2"] for f in fits],
        ".-",
        color="#164f54",
        lw=1.6,
        label="DESI + BBN",
    )
    a.scatter(
        [threshold["q"]],
        [detail["delta_chi2"]],
        color="#b84c2b",
        marker="D",
        zorder=4,
        label="Minimum for all rings",
    )
    a.set(
        title="Universal trail mixture",
        xlabel="Paid cold-to-baryon ratio",
        ylabel=r"Profile $\Delta\chi^2$",
        xlim=(4.15, 6.9),
        ylim=(-0.8, 33),
    )
    a.legend(fontsize=8, loc="upper left", frameon=False)
    meta, groups = rings.o.catalogue()
    kept, _ = rings.o.select(meta, groups)
    rows = kept[name]
    r = np.array([v["rad_kpc"] for v in rows])
    data = np.array([v["vobs_kms"] for v in rows])
    errors = np.array([v["e_vobs_kms"] for v in rows])
    b.errorbar(
        r,
        data,
        yerr=2 * errors,
        fmt="o",
        color="#222222",
        markersize=3.5,
        capsize=2.5,
        label="SPARC, ±2 errors",
        zorder=5,
    )
    for cap, color, style, label in (
        (radial["caps"]["cosmic_95"], "#164f54", "--", "Cosmic cap 5.975"),
        (detail["minimum_cap_bracket"][1], "#b84c2b", "-", "Feasible cap 6.785"),
    ):
        f = rings.shape(rows, meta[name], rings.o.BOXES["two_errors"], cap)
        vel2 = [
            f["scaled_stellar"] * rings.o.signed_square(row["vdisk_kms"])
            + f["scaled_gas"] * rings.o.signed_square(row["vgas_kms"])
            + 1e9 * rings.o.G * h / radius
            for row, h, radius in zip(rows, f["scaled_enclosed_mass"], r, strict=True)
        ]
        b.plot(r, np.sqrt(vel2), ls=style, lw=1.9, color=color, label=label)
    b.set(
        title="F571-8", xlabel="Catalogue radius (kpc)", ylabel="Circular speed (km/s)"
    )
    b.legend(fontsize=8, loc="lower right", frameon=False)
    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(root / "mixture.pdf", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    plot()
