"""Fixed-volume repair on shrinking balls with one sector ensemble across observers."""

import argparse
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from ..fixed_volume.measure import spatial_measures
from ..frame import FRAME


HERE = Path(__file__).resolve().parent
CURVATURES = (-1.0, -0.5, 0.0, 0.5, 1.0)
REFERENCE = 2
RADII = (0.25, 0.125, 0.0625, 0.03125)
BRANCHES = ('maximally_symmetric', 'curved_plane_product')
LAWS = ('reversible', 'clear_on_contraction')
SEED = 712952


@lru_cache(None)
def quadrature(stretch_squared, order):
    mu, mu_weights = np.polynomial.legendre.leggauss(order)
    phi = 2 * np.pi * np.arange(2 * order) / (2 * order)
    angular_factor = stretch_squared * np.cos(phi)**2 + np.sin(phi)**2
    projected_length = np.sqrt((1 - mu[:, None]**2) * angular_factor)
    transverse_factor = stretch_squared / angular_factor
    angular_weights = mu_weights[:, None] / (4 * order)
    radial_nodes, radial_weights = np.polynomial.legendre.leggauss(order)
    return (projected_length, transverse_factor, angular_weights,
            (1 + radial_nodes) / 2, radial_weights / 2)


def product_measures(curvature, radius, stretch_squared, order=12):
    """Volume and area on the spacelike exponential slice in R^(1,1) x M_K^2."""
    if radius == 0 or curvature == 0:
        return spatial_measures(0, radius)
    length, transverse, weights, nodes, radial_weights = quadrature(stretch_squared, order)
    radii = radius * np.append(nodes, 1.0)
    angle = math.sqrt(abs(curvature)) * radii[:, None, None] * length
    scale = (np.sin(angle) if curvature > 0 else np.sinh(angle)) / angle
    density = np.sqrt(1 + (scale**2 - 1) * transverse)
    area = 4 * math.pi * radii**2 * np.sum(weights * density, axis=(1, 2))
    return float(radius * (radial_weights @ area[:-1])), float(area[-1])


def geometry(branch, curvature, radius, u):
    if branch == 'maximally_symmetric':
        return spatial_measures(curvature, radius)
    if branch == 'curved_plane_product':
        return product_measures(curvature, radius, float(1 + u[1]**2 + u[2]**2))
    raise ValueError('unknown comparison geometry')


def equal_volume_surface(branch, curvature, ell, u, target):
    radius = brentq(lambda r: geometry(branch, curvature, r, u)[0] - target,
                    0, 1.4 * ell, xtol=1e-14)
    return geometry(branch, curvature, radius, u)[1]


def sector_kernel(full_tiles, alphabet, law):
    """Uniformly refresh every allowed symbol, propose a uniform sector, then refresh the bit."""
    counts = np.asarray(full_tiles, dtype=np.int64)
    loss = np.maximum(counts[:, None] - counts[None, :], 0)
    if law == 'reversible':
        acceptance = np.exp(-loss * math.log1p(1 / alphabet))
    elif law == 'clear_on_contraction':
        acceptance = np.ones(loss.shape)
    else:
        raise ValueError('unknown repair')
    transition = acceptance / len(counts)
    transition[np.diag_indices(len(counts))] += 1 - transition.sum(axis=1)
    return transition


def terminal_records(kernel, trials, rng):
    records = []
    for launch in (0, len(kernel) - 1):
        probabilities = np.linalg.matrix_power(kernel, 1024)[launch]
        probabilities /= probabilities.sum()
        outcome = np.repeat(probabilities / 2, 2)
        counts = np.zeros(len(outcome), dtype=np.int64)
        remaining = trials
        while remaining:
            size = min(remaining, 2**50)
            counts += rng.multinomial(size, outcome)
            remaining -= size
        records.append(counts.reshape(-1, 2).tolist())
    return records


def survey(branch, ell, u, budget, rng):
    target_volume, flat_area = spatial_measures(0, ell)
    records = []
    for curvature in CURVATURES:
        stopped_volume = rng.gamma(budget) * target_volume / budget
        area = equal_volume_surface(branch, curvature, ell, u, stopped_volume)
        records.append(dict(held_interior=int(rng.poisson(budget * stopped_volume / target_volume)),
                            surface=int(rng.poisson(budget * area / flat_area))))
    return records


def read(ell, budget, geometry_records, sector_records):
    flat_area = 4 * math.pi * ell**2
    surface = np.asarray([r['surface'] for r in geometry_records])
    area_slope = (surface - surface[REFERENCE]) * (flat_area / budget)
    counts = np.asarray(sector_records)
    bulk_slope = np.log(counts[REFERENCE, 0] / counts)
    total = (bulk_slope + area_slope[:, None]) / ell**4
    return dict(bulk_slopes=bulk_slope.tolist(), area_slopes=area_slope.tolist(),
                residuals_over_ell4=total.tolist())


def panel(branch, ell, level, branch_id):
    def stream(*keys):
        return np.random.default_rng(np.random.SeedSequence([SEED, branch_id, level, *keys]))

    alphabet = int(ell**-6)
    budget, trials = 2**(26 + 8 * level), 2**(26 + 10 * level)
    target = spatial_measures(0, ell)[0]
    source_areas = np.array([equal_volume_surface(branch, k, ell, FRAME[0], target)
                             for k in CURVATURES])
    tiles = np.floor(alphabet * source_areas).astype(np.int64)
    source_records = {}
    for law_id, law in enumerate(LAWS):
        kernel = sector_kernel(tiles, alphabet, law)
        source_records[law] = terminal_records(kernel, trials, stream(1, law_id))
    observers = []
    max_z = 0.
    maxima = {law: 0. for law in LAWS}
    max_area = 0.
    for frame_id, u in enumerate(FRAME):
        records = survey(branch, ell, u, budget, stream(2, frame_id))
        max_z = max(max_z, *(abs(r['held_interior'] - budget) / math.sqrt(2 * budget)
                             for r in records))
        interpretations = {}
        for law in LAWS:
            interpretations[law] = [read(ell, budget, records, launch)
                                    for launch in source_records[law]]
            for observed in interpretations[law]:
                maxima[law] = max(maxima[law], float(np.max(np.abs(observed['residuals_over_ell4']))))
                max_area = max(max_area, float(np.max(np.abs(observed['area_slopes']))) / ell**4)
        observers.append(dict(frame_id=frame_id, timelike=[str(x) for x in u],
                              records=records, interpretations=interpretations))
    return dict(branch=branch, ell=ell, ordinary_symbols=alphabet, geometry_budget=budget,
                sector_trials_per_launch=trials, full_tiles=tiles.tolist(),
                terminal_sector_bit_counts=source_records, observers=observers,
                max_held_interior_z=max_z, max_abs_area_over_ell4=max_area,
                max_abs_residual_over_ell4=maxima)


def run():
    panels = [panel(branch, ell, level, branch_id) for branch_id, branch in enumerate(BRANCHES)
              for level, ell in enumerate(RADII)]
    return dict(seed=SEED, protocol_sha256=hashlib.sha256(
                (HERE / 'REFINEMENT_PROTOCOL.md').read_bytes()).hexdigest(),
                curvature_sectors=CURVATURES, repair_attempts=1024, panels=panels)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='record a fresh measurement')
    args = parser.parse_args()
    result = run()
    for row in result['panels']:
        print(row['branch'], row['ell'], row['max_held_interior_z'], row['max_abs_residual_over_ell4'])
    if args.output:
        with args.output.open('x') as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write('\n')
    elif json.loads(json.dumps(result)) != json.loads((HERE / 'refinement.json').read_text()):
        raise SystemExit('measurement differs from recorded receipt')
