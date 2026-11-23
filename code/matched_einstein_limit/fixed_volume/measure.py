"""Observe geodesic-ball area contrast at fixed proper volume in ultrastatic geometries."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


HERE = Path(__file__).resolve().parent
RADII = (0.5, 0.25, 0.125, 0.0625)
CURVATURES = (-1.0, -0.5, 0.0, 0.5, 1.0)
SEED = 2026112315


def spatial_measures(curvature, radius):
    """Proper three-volume and boundary area on the forward source."""
    if curvature == 0:
        return 4 * math.pi * radius**3 / 3, 4 * math.pi * radius**2
    scale = math.sqrt(abs(curvature))
    angle = scale * radius
    if curvature > 0:
        area_radius = math.sin(angle) / scale
        volume = math.pi * (2 * angle - math.sin(2 * angle)) / scale**3
    else:
        area_radius = math.sinh(angle) / scale
        volume = math.pi * (math.sinh(2 * angle) - 2 * angle) / scale**3
    return volume, 4 * math.pi * area_radius**2


def clock_triangle(curvature, ell, tick):
    # Cosine law for equal perpendicular radial sides in a totally geodesic two-plane.
    if curvature == 0:
        link = math.sqrt(2) * ell
    elif curvature > 0:
        k = math.sqrt(curvature)
        link = 2 * math.asin(math.sin(k * ell) / math.sqrt(2)) / k
    else:
        k = math.sqrt(-curvature)
        link = 2 * math.asinh(math.sinh(k * ell) / math.sqrt(2)) / k
    return dict(radial_ticks=round(2 * ell / tick), mirror_link_ticks=round(2 * link / tick))


def acquire(curvature, ell, counts, tick, rng):
    target_volume, flat_area = spatial_measures(0, ell)
    interior_density = counts / target_volume
    surface_density = counts / flat_area
    # After intensity transformation, ordered radial hits have exponential spacings.
    arrival_volume = rng.gamma(shape=counts, scale=1.0) / interior_density
    stopped_radius = brentq(lambda r: spatial_measures(curvature, r)[0] - arrival_volume,
                            0, 2 * ell, xtol=1e-14)
    regions = {}
    for label, radius in (('volume_threshold', stopped_radius), ('fixed_radius', ell)):
        volume, area = spatial_measures(curvature, radius)
        regions[label] = dict(radial_round_trip_ticks=round(2 * radius / tick),
                              held_interior_count=int(rng.poisson(interior_density * volume)),
                              surface_count=int(rng.poisson(surface_density * area)))
    return dict(regions=regions, triangle=clock_triangle(curvature, ell, tick))


def interpret(ell, count_budget, tick, record, reference):
    flat_area = 4 * math.pi * ell**2
    area_unit = flat_area / count_budget
    coefficient = 4 * math.pi * ell**4 / 15
    triangle = record['triangle']
    side = tick * triangle['radial_ticks'] / 2
    diagonal = tick * triangle['mirror_link_ticks'] / 2
    # b² = 2r² - K r⁴/3 + O(r⁶) for a right-angle equal-sided triangle.
    sectional = 3 * (2 * side**2 - diagonal**2) / side**4
    inferred = dict(radar_einstein_contraction=3 * sectional)
    for label, region in record['regions'].items():
        difference = (region['surface_count'] - reference['regions'][label]['surface_count']) * area_unit
        inferred[label] = dict(area_contrast=difference, area_einstein_contraction=-difference / coefficient,
                               held_volume_count_standardized=(region['held_interior_count'] - count_budget)
                               / math.sqrt(2 * count_budget))
    return inferred


def run():
    panels = []
    for level, ell in enumerate(RADII):
        count_budget = 2**(28 + 8 * level)
        tick = 2**-18 * ell**5
        def stream(index):
            return np.random.default_rng(np.random.SeedSequence([SEED, level, index]))
        baseline = acquire(0, ell, count_budget, tick, stream(0))
        cases = []
        for i, curvature in enumerate(CURVATURES, start=1):
            record = acquire(curvature, ell, count_budget, tick, stream(i))
            inferred = interpret(ell, count_budget, tick, record, baseline)
            cases.append(dict(prepared_sectional_curvature=curvature, record=record, inferred=inferred))
        panels.append(dict(flat_radius=ell, count_budget=count_budget, proper_clock_tick=tick,
                           flat_comparison=baseline, cases=cases))
    return dict(protocol_sha256=hashlib.sha256((HERE / 'PROTOCOL.md').read_bytes()).hexdigest(),
                seed=SEED, panels=panels)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write a new receipt')
    options = parser.parse_args()
    result = run()
    if options.output:
        with options.output.open('x') as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    elif result != json.loads((HERE / 'measurements.json').read_text()):
        raise SystemExit('fixed-volume measurements differ from receipt')
    for panel in result['panels']:
        residual = max(abs(case['inferred']['volume_threshold']['area_einstein_contraction'] -
                           case['inferred']['radar_einstein_contraction']) for case in panel['cases'])
        print(f"radius {panel['flat_radius']:.4f}, largest area-radar contraction difference {residual:.7f}")
