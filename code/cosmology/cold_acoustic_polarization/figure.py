"""Plot the independent held polarization spectra and the source density grid."""

import argparse
import json

import matplotlib.pyplot as plt

from .measure import HERE


def plot(report, trails=None):
    bands = report["selected_bands_uK2"]
    predictions = bands["predictions"].copy()
    styles = [
        ("federation_matching", "External matching cap 3/7", "#ad5828", ":"),
        ("max_service", "Single-service cap 1", "#186b76", "--"),
        ("standard_control", "Standard cold ratio 5.36", "#282c79", "-"),
    ]
    title = "Planck polarization predicted from temperature calibration"
    output = "polarization.pdf"
    if trails is not None:
        for board in trails["boards"]:
            predictions[str(board["owners"])] = board["predicted_bands_uK2"]
        styles = [
            ("20", "20 owners, ratio 2.10", "#b25729", ":"),
            ("320", "320 owners, ratio 10.52", "#805391", "-."),
            ("standard_control", "Standard ratio 5.36", "#202331", "-"),
            ("80", "80 owners, ratio 5.13", "#137e77", "--"),
        ]
        title = "Planck polarization for three commissioned trail lengths"
        output = "trail_polarization.pdf"
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, layout="constrained")
    for i, (ax, ch) in enumerate(zip(axes, ("TE", "EE"), strict=True), start=1):
        ax.errorbar(
            bands["ell"],
            bands["measured"][i],
            yerr=bands["sd"][i],
            fmt=".",
            color="#3f454b",
            alpha=0.65,
            lw=0.55,
            ms=2.3,
            zorder=1,
            label="Planck PR3",
        )
        for key, label, color, style in styles:
            ax.plot(
                bands["ell"],
                predictions[key][i],
                style,
                color=color,
                lw=1.55,
                label=label,
            )
        ax.set_ylabel(rf"$D_\ell^{{{ch}}}$ ($\mu\mathrm{{K}}^2$)")
        ax.grid(axis="y", alpha=0.17)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(100, 1500)
    axes[0].axhline(0, color=".78", lw=0.7, zorder=0)
    axes[-1].set_xlabel(r"Multipole $\ell$")
    axes[0].legend(ncol=2, fontsize=8.5, loc="upper right", frameon=True)
    fig.suptitle(title, fontsize=13)
    fig.savefig(HERE / output)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trails", action="store_true")
    args = parser.parse_args()
    trails = (
        json.loads((HERE / "trail_results.json").read_text()) if args.trails else None
    )
    plot(json.loads((HERE / "results.json").read_text()), trails)
