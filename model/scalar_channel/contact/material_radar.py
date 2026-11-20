"""Replay the coupled pulse, material clock, and six-station experiments offline."""

import argparse
import json
from pathlib import Path

import numpy as np

from ..carrier import require
from . import coherent_radar, contact_response, intervention, joint_clock, radar, station_clock
from .__main__ import HERE, check_sources, digest, load_json, protocol as check_contact

SOURCE_FILES = {"model/scalar_channel/contact/" + name for name in (
    "sources.json", "radar_protocol.json", "joint_source.py", "joint_clock.py", "radar.py",
    "coherent_radar.py", "contact_response.py", "station_clock.py", "intervention.py", "material_radar.py")}


def protocol():
    check_contact()
    check_sources(HERE / "radar_sources.json", SOURCE_FILES)
    return load_json(HERE / "radar_protocol.json")


def experiment(spec, measured=None):
    radar_report, radar_data = radar.run(spec)
    if measured is None:
        clock_report, clock_data = joint_clock.run(spec)
    else:
        with np.load(HERE / "clock_weights.npz", allow_pickle=False) as weights:
            clock_report, clock_data = joint_clock.run(spec, acquisition=(measured, weights))
        require(clock_report["sampling_within_bound"], "archived counts fail the simultaneous source-law test")
    coherent_report, amplitudes = coherent_radar.run(spec)
    # On tracing the coherent mark, every individual sector must have the same
    # outbound/return law as the separately executed stationary-mark program.
    labels = amplitudes["coherent_radar_labels"]
    sectors = amplitudes["coherent_radar_sector_weights"]
    for i, (origin, target) in enumerate(amplitudes["coherent_radar_pairs"]):
        a, b = (spec["endpoints"].index(int(x)) for x in (origin, target))
        expected = radar_data["radar_joint_probabilities"][:, a, b]
        np.testing.assert_allclose(sectors[i], expected[:, labels[:, 0]-1, labels[:, 1]-1], rtol=0, atol=5e-11)
    pulse_report, pulse_counts = intervention.run(spec)
    local_counts = station_clock.acquire(spec, radar_data["radar_leg_probabilities"])
    emission = station_clock.emission_counts(spec)
    counts = {k: v for k, v in clock_data.items() if k.endswith("_counts") and k != "coarse_counts"}
    counts.update(radar_joint_counts=radar_data["radar_joint_counts"], local_clock_counts=local_counts,
                  emission_clock_counts=emission, intervention_counts=pulse_counts)
    # These are sufficient statistics for every coherent-radar history and
    # quadrature, including coherent residual vectors' overlap on timeout.
    for key in ("coherent_radar_sector_weights", "coherent_radar_overlap",
                "coherent_radar_labels", "coherent_radar_pairs"):
        counts[key] = amplitudes[key]
    return {"material_and_propagation": pulse_report, "exact_response": contact_response.run(spec),
            "coherent_roundtrip": coherent_report, "radar": radar_report, "joint_clock": clock_report,
            "receiving_clock": station_clock.score(spec, local_counts, emission)}, counts


def verify_counts(path, expected):
    with np.load(path, allow_pickle=False) as saved:
        require(len(saved.files) == len(expected) and set(saved.files) == set(expected),
                "missing, duplicate or unexpected count family")
        for key, sampled in expected.items():
            counts = saved[key]
            require(counts.dtype == sampled.dtype and counts.shape == sampled.shape,
                    f"recorded array type or shape differs: {key}")
            if key in ("coherent_radar_sector_weights", "coherent_radar_overlap"):
                require(np.allclose(counts, sampled, rtol=0, atol=5e-11), f"coherent record differs: {key}")
            else:
                require(counts.dtype == np.int64 and np.array_equal(counts, sampled),
                        f"recorded counts or labels differ from source replay: {key}")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--record", type=Path, help="acquire all counts at a new path")
    mode.add_argument("--read-archive", action="store_true",
                      help="verify the archived multinomial acquisition with its source laws and read it directly")
    cli.add_argument("--details", action="store_true", help="print all fold and control scores")
    args = cli.parse_args()
    spec = protocol()
    path = HERE / "radar_counts.npz"
    if args.record is None:
        require(digest(path) == path.with_suffix(".sha256").read_text(encoding="ascii").strip(),
                "archived count digest differs")
    if args.read_archive:
        with np.load(path, allow_pickle=False) as measured:
            report, counts = experiment(spec, measured)
    else:
        report, counts = experiment(spec)
    if args.record is not None:
        with args.record.open("xb") as f:
            np.savez_compressed(f, **counts)
        verify_counts(args.record, counts)
    else:
        verify_counts(path, counts)
    if not args.details:
        clock = report["joint_clock"]
        report = {"local_intervention": report["material_and_propagation"],
                  "mark_visible_dimensions": report["exact_response"]["mark_cyclic_dimension"],
                  "selected_on_five_stations": clock["held_out_model_selection"]["selected_models"],
                  "coherent_return_vs_reset_history_TV":
                      report["coherent_roundtrip"]["maximum_history_TV_from_independent_repreparation"],
                  "successful_3d_folds": [s["geometry_from_born_laws"]["accepted_pairs"]
                                          for s in report["radar"]["settings"]]}
        if args.read_archive:
            acquisition = clock["original_acquisition"]
            families = acquisition["families"].values()
            report["archived_clock_laws"] = {
                "weight_digest": acquisition["weight_digest"],
                "locally_requantized_bins_different": sum(f["locally_requantized_bins_different"] for f in families),
                "maximum_source_total_variation": max(f["maximum_total_variation"] for f in families)}
    print(json.dumps(report, allow_nan=False, indent=2))


if __name__ == "__main__":
    main()
