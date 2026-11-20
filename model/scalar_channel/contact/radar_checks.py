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


def precision_acquisition(spec, measurements, saved=None):
    """New fixed-budget cohorts. Multinomial batches sample sufficient statistics."""
    weights = measurements["radar_acquisition_weights"]
    counts = np.zeros_like(weights)
    if saved is not None:
        require(saved.dtype == np.int64 and saved.shape == weights.shape and
                np.all((0 <= saved) & (saved <= 2**34)) and
                np.all(saved[:, np.arange(6), np.arange(6)] == 0),
                "invalid archived precision counts")
    log_tail = np.log(2 * 60 * 169 / 1e-6)
    confidence = []
    for mark in (0, 1):
        for a, b in np.ndindex(6, 6):
            if a == b:
                continue
            rng = np.random.Generator(np.random.PCG64([spec["seed"], 8639, mark, a, b]))
            p = weights[mark, a, b].ravel() / radar.LATTICE
            if saved is None:
                draws = rng.multinomial(2**28, p, size=64)
                counts[mark, a, b] = draws.sum(axis=0).reshape(13, 13)
            else:
                counts[mark, a, b] = saved[mark, a, b]
            row = counts[mark, a, b].ravel()
            tolerance = np.sqrt(2 * 2**34 * p*(1-p)*log_tail) + 2*log_tail/3
            require(row.sum() == 2**34 and np.all(row[p == 0] == 0) and
                    np.all(abs(row - 2**34*p) <= tolerance),
                    "precision records fail the source-law sampling test")
        mean, hits = radar.distance_statistics(counts[mark], spec["dt"])
        bias = radar._lattice_conditional_bias(weights[mark], counts[mark], spec["dt"], spec["ticks"])
        _, lo, hi = radar.confidence_envelopes(mean, hits, spec["dt"], spec["ticks"], bias)
        confidence.append((mean, lo, hi))
    return counts, confidence


def confirm(spec, saved=None):
    exact = certified_ranges.run(spec)
    radar_summary, measurements = radar.run(spec)
    with np.load(HERE / "radar_counts.npz", allow_pickle=False) as original_radar:
        require(np.array_equal(original_radar["radar_joint_counts"], measurements["radar_joint_counts"]),
                "geometry confidence intervals do not use the original radar acquisition")
    fits = []
    precision_counts, precision = precision_acquisition(spec, measurements, saved)
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
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--record", type=Path, help="acquire tagged records at a new path")
    mode.add_argument("--read-archive", action="store_true", help="score original precision counts directly")
    args = cli.parse_args()
    check_sources(HERE / "radar_checks_sources.json", SOURCES)
    require(digest(HERE / "radar_counts.npz") == (HERE / "radar_counts.sha256").read_text().strip(),
            "radar receipt differs")
    spec = material_radar.protocol()
    path = HERE / "tagged_counts.npz"
    if args.record is None:
        require(digest(path) == path.with_suffix(".sha256").read_text().strip(), "tagged receipt differs")
    if args.read_archive:
        with np.load(path, allow_pickle=False) as receipt:
            report, counts = confirm(spec, receipt["precision_counts"])
    else:
        report, counts = confirm(spec)
    if args.record:
        with args.record.open("xb") as f:
            np.savez_compressed(f, **counts)
        material_radar.verify_counts(args.record, counts)
    else:
        material_radar.verify_counts(path, counts)
    if not args.details:
        report = {"excluded_relative_errors": [[x["relative_tolerance"] for x in row["source_enclosures"]
                                                if x["global_3d_fit_excluded"]] for row in report["geometry"]],
                  "common_protected_rings": sum(all(r["retained"]) for r in report["cages"]),
                  "tagged_time_counts": int(counts["tagged_time_counts"].sum())}
    print(json.dumps(report, allow_nan=False, indent=2))


if __name__ == "__main__":
    main()
