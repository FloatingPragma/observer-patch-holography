"""Measured flux-clock delay and matched line widths in ten galaxies."""

import json
import matplotlib.pyplot as plt

from .closed_inventory import HERE, LAB


def run():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(8.2, 3.3), layout="constrained")
    clock = json.loads((LAB / "renewal_clock.json").read_bytes())["readers"]
    for t, row in clock.items():
        a, b = (row["arms"][k] for k in ("original", "repaired"))
        for dx, field, ci, color, label in [
            (
                -0.1,
                "predicted_slots_per_tick",
                "slots_interval_99",
                "#a34122",
                "Source held until tick",
            ),
            (
                0.1,
                "renew_each_attempt_slots",
                "fresh_attempt_slots_interval_99",
                "#215a7a",
                "Fresh source each attempt",
            ),
        ]:
            pred = 100 * (b[field] / a[field] - 1)
            lo = 100 * (b[ci][0] / a[ci][1] - 1)
            hi = 100 * (b[ci][1] / a[ci][0] - 1)
            ax.errorbar(
                int(t) + dx,
                pred,
                yerr=[[pred - lo], [hi - pred]],
                fmt="o",
                ms=4,
                capsize=3,
                color=color,
                label=label if t == "1" else None,
            )
    ax.axhline(0, lw=0.8, color="#888888")
    ax.set(
        xticks=(1, 2),
        xticklabels=("Single", "Double"),
        xlim=(0.6, 2.5),
        ylabel="Change in mean time per tick (%)",
    )
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    rows = json.loads((HERE / "clock_width.json").read_bytes())["galaxies"]
    letters = {"DDO154": "DDO 154", "NGC2841": "NGC 2841", "NGC3521": "NGC 3521"}
    for n, tab in rows.items():
        r = tab["outer_rings"][0]["renewal"]
        x, y = r["observed_rms_kms"], r["required_rms_kms"]
        bx.scatter(x, y, s=26, color="#215a7a")
        if n in letters:
            bx.annotate(
                letters[n],
                (x, y),
                xytext=(4, 1),
                textcoords="offset points",
                fontsize=7.5,
            )
    bx.plot([0, 50], [0, 50], color="#777777", ls="--", lw=1, label="Equal widths")
    bx.set(
        xlim=(0, 48),
        ylim=(0, 225),
        xlabel="Detected H I width (km/s)",
        ylabel="Minimum renewal width (km/s)",
    )
    bx.legend(frameon=False, fontsize=8, loc="upper left")
    for ax_ in (ax, bx):
        ax_.spines[["top", "right"]].set_visible(False)
        ax_.set_axisbelow(True)
        ax_.grid(axis="y", alpha=0.18)
    fig.savefig(
        HERE / "clock_width.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)


if __name__ == "__main__":
    run()
