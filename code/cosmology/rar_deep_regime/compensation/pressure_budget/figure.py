"""Show conserved prepared stock alongside independent observational abundance scales."""

import json
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = Path(__file__).parent


def draw():
    observed = json.loads((HERE / "results.json").read_text())
    chains = json.loads(
        (
            HERE.parents[2] / "cold_acoustic_polarization" / "chain_results.json"
        ).read_text()
    )
    fig, axes = plt.subplots(
        2, 1, figsize=(7.2, 4.9), sharex=True, gridspec_kw={"height_ratios": [1, 1.42]}
    )
    fig.subplots_adjust(left=0.23, right=0.97, top=0.84, bottom=0.14, hspace=0.40)
    cosmic, lens = axes
    for row in observed["commissioned_trails"]:
        q = float(Fraction(row["conserved_ratio"]))
        for axis in axes:
            axis.axvline(
                q, color="#978171", linestyle=(0, (3, 3)), linewidth=1.1, zorder=1
            )
        cosmic.text(
            q,
            1.18,
            f"{row['live_seams']} records\n{q:.2f}",
            transform=cosmic.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=8,
        )
    for i, k in enumerate(("EE", "TE")):
        low, lo, mid, hi, high = chains[k]["ratio_quantiles_025_16_50_84_975"]
        cosmic.plot([low, high], [i, i], color="#245c80", lw=2)
        cosmic.plot([lo, hi], [i, i], color="#245c80", lw=7, solid_capstyle="butt")
        cosmic.plot(mid, i, "o", color="white", mec="#245c80", ms=5)
    cosmic.set_yticks([0, 1], ["EE + lowE", "TE + lowE"])
    cosmic.set_ylim(-0.48, 1.48)
    cosmic.set_title(
        "Planck cold-to-baryon ratio · 68% and 95% credible intervals",
        loc="left",
        fontsize=9.5,
    )
    thresholds = observed["cases"]["two_systematic_double_baryons"][
        "additional_mass_over_baryons"
    ]
    for i, q in enumerate(thresholds):
        lens.annotate(
            "",
            xy=(11.15, i),
            xytext=(q, i),
            arrowprops=dict(arrowstyle="->", color="#17675e", lw=1.45),
        )
        lens.plot(q, i, marker="|", ms=10, mew=1.6, color="#17675e")
        lens.text(
            q - 0.12,
            i,
            f"{q:.2f}",
            ha="right",
            va="center",
            fontsize=8.5,
            color="#14584e",
        )
    lens.set_yticks(range(4), ["9.0–10.5", "10.5–10.8", "10.8–11.1", "11.1–11.5"])
    lens.set_ylim(3.55, -0.55)
    lens.set_ylabel("Catalogue log₁₀(Mᵦ/M☉)", fontsize=9)
    lens.set_title("Galaxy lensing · lower bounds at 291 kpc", loc="left", fontsize=9.5)
    lens.set_xlabel(
        "Additional mass-equivalent energy per unit of ordinary matter", fontsize=10
    )
    for axis in axes:
        axis.set_xlim(0, 11.35)
        axis.xaxis.set_major_locator(MultipleLocator(1))
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=8.5)
    fig.savefig(
        HERE / "abundance.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)


if __name__ == "__main__":
    draw()
