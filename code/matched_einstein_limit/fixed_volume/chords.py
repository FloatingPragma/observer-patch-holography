"""Recover volumes from cosine-weighted boundary travel on the ultrastatic comparison."""
import argparse
import json
import math

import numpy as np
from numpy.polynomial.legendre import leggauss

from . import measure


def clock_tapes(curvature, radius, tick, cosines):
    """First-exit round trips in a totally geodesic radial plane, at fixed incidence."""
    if curvature == 0:
        half_chords = radius * cosines
    elif curvature > 0:
        k = math.sqrt(curvature)
        if k * radius >= math.pi / 2:
            raise ValueError("Positive-curvature aperture must be strictly convex")
        half_chords = np.arctan(np.tan(k * radius) * cosines) / k
    else:
        k = math.sqrt(-curvature)
        half_chords = np.arctanh(np.tanh(k * radius) * cosines) / k
    # A reflected light ray in these ultrastatic units clocks twice its spatial length.
    return np.rint(4 * half_chords / tick).astype(np.int64)


def boundary_volume(area, tick, tapes, cosines, weights):
    lengths = tick * tapes / 2
    mean = float((weights * 2 * cosines) @ lengths)
    return area * mean / 4, mean


def run():
    data = json.loads((measure.HERE / "measurements.json").read_text())
    nodes, weights = leggauss(16)
    cosines, weights = (1 + nodes) / 2, weights / 2
    rows = []
    for panel in data["panels"]:
        ell = panel["flat_radius"]
        tick = panel["proper_clock_tick"]
        budget = panel["count_budget"]
        target, flat_area = measure.spatial_measures(0, ell)
        errors = []
        area_errors = []
        decoded = []
        for arm in panel["cases"]:
            k = arm["prepared_sectional_curvature"]
            counts = arm["record"]["regions"]["volume_threshold"]
            radius = counts["radial_round_trip_ticks"] * tick / 2
            measured_area = counts["surface_count"] * flat_area / budget
            tapes = clock_tapes(k, radius, tick, cosines)
            volume, mean = boundary_volume(measured_area, tick, tapes, cosines, weights)
            check = target * counts["held_interior_count"] / budget
            errors.append(abs(volume - check) / target)
            exact_volume, exact_area = measure.spatial_measures(k, radius)
            area_errors.append(abs(volume / measured_area - exact_volume / exact_area) / radius)
            decoded.append(dict(sectional_curvature=k, chord_round_trip_ticks=tapes.tolist(),
                                mean_length=mean, boundary_volume=volume, interior_count_volume=check))
        rows.append(dict(reference_radius=ell, maximum_fractional_difference=max(errors),
                         maximum_clock_ratio_error_per_radius=max(area_errors), cases=decoded))
    return dict(angular_quadrature_order=len(nodes), panels=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = run()
    dest = measure.HERE / "chords.json"
    if parser.parse_args().write:
        dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from model.experiments.interval_echo_einstein_station.collar_pairing.measure import compare
        compare(result, json.loads(dest.read_text()))
    print([p["maximum_fractional_difference"] for p in result["panels"]])
