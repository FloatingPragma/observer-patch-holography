"""Certify the radar distances and replay an active and a silent time-bin reflector."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from ..carrier import require
from ..cage import construction as cage
from . import certified_ranges, metric_certificate, radar, tagged_radar
from . import material_radar
from .__main__ import HERE, check_sources, digest
from .joint_source import radar_source

SOURCES = {"model/scalar_channel/contact/" + name for name in (
    "radar_checks.py", "tagged_radar.py", "certified_ranges.py", "metric_certificate.py", "radar_sources.json")}


def preserved_cages(spec):
    """Check exact reducing pairs in both mark sectors, against the full source."""
    car = cage.cell()
    sources = [radar_source(spec, mark) for mark in (0, 1)]
    masks = np.zeros((2, len(car["edges"])), dtype=int)
    masks[:, spec["endpoints"]] = 1
    masks[1, spec["mark"]] = 1
    rows = []
    for ring in cage.rings(car):
        profile = cage.ring_profile(car, ring)
        frame, _, _ = cage.modes(car, profile)
        support = np.flatnonzero(profile)
        counts = masks[:, support].sum(axis=1)
        retained = (counts == 0) | (counts == 4)
        for h, occupied in zip(sources, counts):
            send = h @ frame
            leakage = send - frame @ (frame.T @ send)
            require(abs(np.sum(leakage**2) - occupied*(4-occupied)/16) < 1e-12,
                    "full-source ring leakage differs from its exact potential variance")
        touches = spec["mark"] in support
        require(not (retained.all() and touches), "protected pair intersects an active mark")
        rows.append({"ring": list(ring), "contacted_seams": counts.tolist(),
                     "retained": retained.tolist(), "contains_mark": touches})
    return rows


def precision_acquisition(spec, measurements):
    """New fixed-budget cohorts. Multinomial batches sample sufficient statistics."""
    weights = measurements["radar_acquisition_weights"]
    counts = np.zeros_like(weights)
    confidence = []
    for mark in (0, 1):
        for a, b in np.ndindex(6, 6):
            if a == b:
                continue
            rng = np.random.Generator(np.random.PCG64([spec["seed"], 8639, mark, a, b]))
            p = weights[mark, a, b].ravel() / radar.LATTICE
            draws = rng.multinomial(2**28, p, size=64)
            counts[mark, a, b] = draws.sum(axis=0).reshape(13, 13)
        mean, hits = radar.distance_statistics(counts[mark], spec["dt"])
        bias = radar._lattice_conditional_bias(weights[mark], counts[mark], spec["dt"], spec["ticks"])
        _, lo, hi = radar.confidence_envelopes(mean, hits, spec["dt"], spec["ticks"], bias)
        confidence.append((mean, lo, hi))
    return counts, confidence


def confirm(spec):
    exact = certified_ranges.run(spec)
    radar_summary, measurements = radar.run(spec)
    with np.load(HERE / "radar_counts.npz", allow_pickle=False) as saved:
        require(np.array_equal(saved["radar_joint_counts"], measurements["radar_joint_counts"]),
                "geometry confidence intervals do not use the original radar acquisition")
    fits = []
    precision_counts, precision = precision_acquisition(spec, measurements)
    for index, record in enumerate(exact["settings"]):
        scale = Fraction(1, exact["scale"])
        low = np.array(record["distance_lower"], dtype=object) * scale
        high = np.array(record["distance_upper"], dtype=object) * scale
        mask = ~np.eye(len(low), dtype=bool)
        require(np.all(np.asarray(record["distance_defined"])[mask]), "undefined exact radar range")
        mean = np.array(radar_summary["settings"][index]["distances_from_born_laws"])
        require(np.all(mean >= low.astype(float)) and np.all(mean <= high.astype(float)),
                "floating-point range escapes its integer enclosure")
        rigorous = metric_certificate.certify((low+high)/2, low, high, numerical_guard=0)
        sampled = metric_certificate.certify(
            radar_summary["settings"][index]["distances_from_counts"],
            measurements["radar_distance_confidence_lower"][index],
            measurements["radar_distance_confidence_upper"][index], confidence_failure_probability=1e-6)
        resolved = metric_certificate.certify(*precision[index], confidence_failure_probability=1e-6)
        fits.append({"source_enclosures": rigorous["thresholds"], "simultaneous_counts": sampled["thresholds"],
                     "precision_counts": resolved["thresholds"]})
    tag_report, archive = tagged_radar.run(spec)
    scale = radar.common_multiplier_intersection(*precision[0][1:], *precision[1][1:])
    return {"certified_ranges": exact, "geometry": fits, "precision_common_scale": scale,
            "cages": preserved_cages(spec), "tagged_echo": tag_report}, {
                "tagged_time_counts": archive["tagged_time_counts"], "precision_counts": precision_counts}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--details", action="store_true")
    cli.add_argument("--record", type=Path, help="acquire tagged records at a new path")
    args = cli.parse_args()
    check_sources(HERE / "radar_checks_sources.json", SOURCES)
    require(digest(HERE / "radar_counts.npz") == (HERE / "radar_counts.sha256").read_text().strip(),
            "radar receipt differs")
    report, counts = confirm(material_radar.protocol())
    if args.record:
        with args.record.open("xb") as f:
            np.savez_compressed(f, **counts)
        material_radar.verify_counts(args.record, counts)
    else:
        path = HERE / "tagged_counts.npz"
        require(digest(path) == path.with_suffix(".sha256").read_text().strip(), "tagged receipt differs")
        material_radar.verify_counts(path, counts)
    if not args.details:
        report = {"excluded_relative_errors": [[x["relative_tolerance"] for x in row["source_enclosures"]
                                                if x["global_3d_fit_excluded"]] for row in report["geometry"]],
                  "common_protected_rings": sum(all(r["retained"]) for r in report["cages"]),
                  "tagged_time_counts": int(counts["tagged_time_counts"].sum())}
    print(json.dumps(report, allow_nan=False, indent=2))


if __name__ == "__main__":
    main()
