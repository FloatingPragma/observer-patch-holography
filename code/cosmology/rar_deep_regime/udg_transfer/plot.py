"""Plot the external transfer of the SPARC deep law to gas-rich diffuse galaxies."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
result = json.loads((HERE / "results.json").read_text())
source = json.loads((HERE / "observations.json").read_text())["table"]
dark = "#243853"
orange = "#b64720"
lower, upper = result["training_interval"]
scale = result["acceleration_scale"]
fig, (sample, profile) = plt.subplots(1, 2, figsize=(9.4, 3.5), layout="constrained")

index = np.arange(len(source))[::-1]
predicted = np.array([r["predicted_kms"] for r in result["six_galaxies"]])
sample.errorbar(
    [r[5] for r in source],
    index,
    xerr=np.array([[r[6], r[7]] for r in source]).T,
    fmt="o",
    color=dark,
    capsize=3,
    label="Measured gas speed",
)
sample.plot(predicted, index, "D", color=orange, label="Limiting-speed prediction")
sample.set_yticks(index, [f"AGC {r[0]}" for r in source])
sample.set(
    xlabel="Circular speed  (km/s)",
    title="Six external gas-rich galaxies",
    xlim=(8, 77),
)
sample.legend(
    frameon=False, fontsize=8, loc="center", bbox_to_anchor=(0.52, -0.26), ncol=2
)

rings = result["resolved_rings"]["figure8"]
r = np.array([x["radius_kpc"] for x in rings])
vb = np.array([x["baryonic_kms"] for x in rings])
v = np.array([x["predicted_kms"] for x in rings])
bands = [np.sqrt(vb**2 + (v**2 - vb**2) * np.sqrt(a / scale)) for a in (lower, upper)]
profile.fill_between(r, *bands, color=orange, alpha=0.14)
profile.plot(r, v, color=orange, label="Transferred deep law")
profile.plot(r, vb, ":", color=dark, label="Baryonic force")
profile.errorbar(
    r,
    [x["observed_kms"] for x in rings],
    yerr=[x["plotted_error_kms"] for x in rings],
    fmt="o",
    ms=4,
    color=dark,
    capsize=3,
    label="Measured gas speed",
)
profile.set(
    xlabel="Radius  (kpc)",
    ylabel="Circular speed  (km/s)",
    title="AGC 114905 at the optical inclination",
    ylim=(0, 72),
    xlim=(0, 10.8),
)
profile.legend(frameon=False, fontsize=8)
for ax in (sample, profile):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#aaaaaa", alpha=0.19)
fig.savefig(HERE / "comparison.pdf", metadata={"CreationDate": None, "ModDate": None})
