"""Check whether finite reversible repair forces fixed-proper-volume entropy stationarity."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.csgraph import connected_components

from ..fixed_volume.measure import acquire, spatial_measures

HERE = Path(__file__).resolve().parent
TILES = (3, 4, 5)
SHOTS = 2**20
COUNTER_BUDGET = 2**32


def geometry_sectors(rng):
    volume, flat_area = spatial_measures(0., 1.)
    tile_area = flat_area / 4

    def boundary_area(curvature):
        radius = brentq(lambda r: spatial_measures(curvature, r)[0]-volume, 0., 2.)
        return spatial_measures(curvature, radius)[1]

    sectors = []
    for tiles in TILES:
        if tiles == 4:
            curvature = 0.
        else:
            bounds = (0., 2.) if tiles == 3 else (-2., 0.)
            curvature = brentq(lambda k: boundary_area(k) - tiles*tile_area, *bounds)
        record = acquire(curvature, 1., COUNTER_BUDGET, 2.**-36, rng)
        measured = record["regions"]["volume_threshold"]
        volume_error = (measured["held_interior_count"]-COUNTER_BUDGET)/np.sqrt(2*COUNTER_BUDGET)
        assert abs(volume_error) < 6
        sectors.append(dict(tiles=tiles, curvature=float(curvature), record=record,
                            area_reading=flat_area*measured["surface_count"]/COUNTER_BUDGET,
                            volume_reading=volume*measured["held_interior_count"]/COUNTER_BUDGET))
    return sectors


def repair_kernel(states, upward_acceptance):
    index = {state: i for i, state in enumerate(states)}
    matrix = np.zeros((len(states), len(states)))
    for i, (n, bulk, bits) in enumerate(states):
        # Nine proposals: quiet, bulk flip, five possible tile flips, two sector changes.
        moves = [(n, bulk, bits), (n, 1-bulk, bits)]
        moves.extend((n, bulk, bits ^ (1 << bit) if bit < n else bits) for bit in range(5))
        for small in (3, 4):
            if n == small:
                changed = n+1
            elif n == small+1 and bits < 2**small:
                changed = n-1
            else:
                changed = n
            moves.append((changed, bulk, bits))
        for destination in moves:
            accept = upward_acceptance if destination[0] > n else 1.
            matrix[i,index[destination]] += accept/len(moves)
            matrix[i,i] += (1-accept)/len(moves)
    np.testing.assert_allclose(matrix.sum(axis=1), 1., atol=1e-14)
    assert connected_components(matrix, directed=True, connection="strong")[0] == 1
    assert np.all(matrix.diagonal() > 0)
    return matrix


def run():
    rng = np.random.default_rng(511423)
    sectors = geometry_sectors(rng)
    states = [(n, bulk, bits) for n in TILES for bulk in range(2) for bits in range(2**n)]
    sector_indices = np.array([2*(n-3)+bulk for n,bulk,_ in states])
    areas = np.repeat([sector["area_reading"] for sector in sectors], 2)
    eta = np.log(2)/np.pi
    trials = []
    for name, acceptance in (("reversible", 1.), ("biased_exposure", .5)):
        kernel = repair_kernel(states, acceptance)
        invariant = np.array([acceptance**n for n,_,_ in states])
        invariant /= invariant.sum()
        np.testing.assert_allclose(invariant[:,None]*kernel,
                                   (invariant[:,None]*kernel).T, atol=1e-15)
        endpoint = np.linalg.matrix_power(kernel, 512)
        for launch in (0, len(states)-1):
            sums = np.bincount(sector_indices, weights=endpoint[launch], minlength=6)
            exact = np.bincount(sector_indices, weights=invariant, minlength=6)
            np.testing.assert_allclose(sums, exact, atol=1e-7, rtol=0)
            observed = rng.multinomial(SHOTS, sums/sums.sum())
            # Fresh copies also read each boundary word, conditional on the reported sector and bit.
            joint_counts = rng.multinomial(SHOTS, endpoint[launch]/endpoint[launch].sum())
            edge_entropies = []
            for group in range(6):
                counts = joint_counts[sector_indices == group]
                conditional = counts[counts > 0]/counts.sum()
                edge_entropies.append(float(-conditional @ np.log(conditional)))
            np.testing.assert_allclose(edge_entropies, np.repeat(TILES, 2)*np.log(2), atol=.001)
            gradient = -np.log(observed/SHOTS) + eta*areas
            residuals = (gradient[1:]-gradient[0]).tolist()
            largest = max(abs(v) for v in residuals)
            if name == "reversible":
                assert largest < .03
            else:
                assert largest > .6
            trials.append(dict(repair=name, initial=states[launch], terminal_counts=observed.tolist(),
                               joint_counts_by_sector_material_bit_edge_word=joint_counts.tolist(),
                               measured_edge_entropies=edge_entropies,
                               five_directional_derivatives=residuals, maximum_absolute_derivative=largest))
    return dict(seed=511423, trials_per_launch=SHOTS, steps=512, surface_density=eta,
                target_volume=4*np.pi/3, geometry=sectors, observed=trials)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import check_replay
        check_replay(json.loads((HERE/"results.json").read_text()), json.loads(json.dumps(result)))
    for row in result["observed"]:
        print(row["repair"], row["initial"], row["maximum_absolute_derivative"])
