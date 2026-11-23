"""Evaluate two native proxies for the per-cut inputs of the OPH acceleration scale."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import measure


def moments(law, score):
    probabilities = law / law.sum()
    mean = float(probabilities @ score)
    squared = (score - mean) ** 2
    variance = float(probabilities @ squared)
    return variance, float(probabilities @ ((squared - variance) ** 2))


def evaluate():
    reader = json.loads((measure.HERE / "results.json").read_text())
    old = json.loads((measure.HERE.parent / "shared_variance.json").read_text())[
        "connected_modular_charge"
    ]
    # Independent full-word cohorts. The microscopic reference in results.json
    # was acquired with a separate seed and does not use any of these outcomes.
    cohorts = [old["rows"][0]["reference_counts"], *old["held_counts"]]
    energy = np.add.outer([0, 32, 33], [0, 32, 33]).ravel()
    laws, _, _ = measure.source_reader()
    ref_counts = np.array(reader["reference"]["counts"])
    empirical = measure.modular_score(ref_counts / ref_counts.sum())
    ideal = measure.modular_score(laws[0])
    rows = []
    for a, b in ((32, 32), (33, 33), (65, 65), (32, 65), (64, 65)):
        mask = (energy == a)[:, None] & (energy == b)[None, :]
        arms = []
        for k, hist in enumerate(cohorts):
            counts = np.asarray(hist)[mask]
            exact, _ = moments(laws[k][mask], ideal[mask])
            fixed, _ = moments(laws[k][mask], empirical[mask])
            sampled, v4 = moments(counts, empirical[mask])
            arms.append(
                dict(
                    arm=old["reference_order"][k],
                    selected_shots=int(counts.sum()),
                    total_shots=int(np.asarray(hist).sum()),
                    population=exact,
                    population_conditional_on_reference=fixed,
                    measured=sampled,
                    standard_error=float(np.sqrt(v4 / counts.sum())),
                )
            )
        rows.append(dict(word_energy_sums=[a, b], results=arms))
    densities = []
    for panel in reader["panels"]:
        measures = []
        for s, scope in zip(
            panel["settings"], panel["geometry"]["annulus_counts"], strict=True
        ):
            arms = s["arms"]["separated"]["reads"][1]
            width = s["relative_radius"] / 2
            n, _, _ = scope
            N = panel["geometry"]["cells"]
            sample = sum(
                (int(ij[0] == "0") + int(ij[1] == "0")) * sum(c)
                for ij, c in arms.items()
            )
            events = sum(sum(v) for v in arms.values())
            up = 2 * n / N
            # Distinct occupied radial bins differ from the marked-word count
            # only if both land in the same annulus and on the same radius.
            pair = n * (n - 1) / (N * (N - 1))
            measures.append(
                dict(
                    relative_radius=s["relative_radius"],
                    width=width,
                    conditional_annulus_station_count=n,
                    expected_words_per_routed_length=up / width,
                    sampled_words_per_routed_length=sample / events / width,
                    distinct_occupied_shells_per_routed_length_bounds=[
                        (up - pair) / width,
                        up / width,
                    ],
                )
            )
        densities.append(dict(cells=panel["geometry"]["cells"], measurements=measures))
    return dict(
        variance_units="nats squared of the connected modular score",
        density_units="occupied native word radii per home-range unit",
        source_scale_units="No independently measured baryonic mass, proper metre, mass-valued anomaly, or universal Newton response enters these records.",
        energy_conditioned=rows,
        mark_density=densities,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = evaluate()
    if args.output:
        with args.output.open("x") as output:
            json.dump(result, output, indent=2, allow_nan=False)
            output.write("\n")
    elif result != json.loads((measure.HERE / "source_scale.json").read_text()):
        raise ValueError("A native scale proxy differs from its archived result")
    print(
        "Conditioned connected-score variances:",
        [r["measured"] for r in result["energy_conditioned"][0]["results"]],
    )
    print(
        "Finest radial density bounds:",
        [
            p["distinct_occupied_shells_per_routed_length_bounds"]
            for p in result["mark_density"][-1]["measurements"]
        ],
    )
