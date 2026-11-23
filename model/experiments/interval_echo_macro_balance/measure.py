"""Coarse star-charge forecasts on the free Interval Echo federation."""

import argparse
import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components, shortest_path
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.executable_model.carrier import Carrier, CELL_SEAMS
from model.executable_model.refinement import DEPOSIT_SHA256, load_refinement
from model.experiments.interval_echo_angular import modular
from model.experiments.interval_echo_einstein_station import cap
from model.experiments.interval_echo_space.shell_source import laplacian, prepare

HERE = Path(__file__).parent
SHOTS = 2**24


def instruments():
    parent, _, _ = modular.common_reference()
    stars, effects, _ = modular.regional_operators(parent)
    xyz = cap.chord_coordinates()
    x, y, z = xyz.T
    design = np.column_stack((np.ones(12), x, y, z, x*x-z*z, y*y-z*z, x*y, x*z, y*z))
    design_receipt = json.loads((modular.HERE / "results.json").read_text())
    trained = design_receipt["angular_design"]["best_nine_example_native_ports"]
    held = sorted(set(range(12)) - set(trained))
    forecast = np.linalg.solve(design[trained].T, design[held].T).T
    orthogonal_residual = np.eye(12) - design @ np.linalg.pinv(design)
    antipodes = np.argmin(xyz @ xyz.T, axis=1)
    np.testing.assert_allclose(orthogonal_residual[antipodes], -orthogonal_residual, atol=1e-13)
    residual = effects[held] - np.einsum("ht,tij->hij", forecast, effects[trained])
    basis = np.column_stack((np.ones(12)/np.sqrt(12), xyz/2))
    projector = basis @ basis.T
    fast = np.eye(12) - projector
    slow_error = max(np.linalg.norm(projector @ r @ projector, 2) for r in residual)
    cross_norm = max(np.linalg.norm(projector @ r @ fast, 2) for r in residual)
    fast_norm = max(np.linalg.norm(fast @ r @ fast, 2) for r in residual)
    local = cell().laplacian().toarray()/2
    np.testing.assert_allclose(basis.T @ basis, np.eye(4), atol=2e-14)
    if np.linalg.eigvalsh(local - 3*fast)[0] < -1e-12 or slow_error > 1e-10:
        raise ValueError("The local kinetic lower bound or angular identity failed")
    return SimpleNamespace(stars=stars, effects=effects, trained=trained, held=held,
                           forecast=forecast, orthogonal_residual=orthogonal_residual, slow_error=slow_error,
                           cross_norm=cross_norm, fast_norm=fast_norm)


def acquire(wave, meter, seed, shots):
    """Pool disjoint cell readings at each star setting. Quiet is a recorded outcome."""
    cells = wave.reshape(-1, 12)
    descendants = len(cells)//20
    measured = np.zeros((20, 12))
    variances = np.zeros_like(measured)
    counters = []
    for pole, (ports, effect) in enumerate(zip(meter.stars, meter.effects, strict=True)):
        rng = np.random.default_rng([1093108, *seed, pole])

        def setting(probabilities, weights):
            grouped = probabilities.reshape(20, descendants, -1).sum(axis=1)
            full = np.r_[grouped.ravel(), max(0., 1-grouped.sum())]
            counts = rng.multinomial(shots, full/full.sum())
            counters.append(counts)
            observed = counts[:-1].reshape(grouped.shape)
            measured[:, pole] += observed @ weights / shots
            variances[:, pole] += (grouped @ (weights**2) - (grouped @ weights)**2)/shots

        setting(abs(cells[:, list(ports)])**2, np.diag(effect)[list(ports)])
        for i, j in combinations(ports, 2):
            chances = np.column_stack((abs(cells[:, i] + cells[:, j])**2,
                                       abs(cells[:, i] - cells[:, j])**2))/2
            setting(chances, np.array([effect[i, j], -effect[i, j]]))
    return measured, variances, np.concatenate(counters)


def parent_graph(car, level):
    contacts = [(a//12, b//12) for a, b in car.external]
    descendants = 4**level
    boundary = Counter((a//descendants, b//descendants) for a, b in contacts
                       if a//descendants != b//descendants)
    expected = {(a//12, b//12): 2**level for a, b in federation(0).external}
    if boundary != expected:
        raise ValueError("The twenty grouped regions have different shared boundaries")
    dual = laplacian(car.cells, contacts)
    adjacency = (sparse.diags(dual.diagonal()) - dual).tocsr()
    for group in range(20):
        region = slice(group * descendants, (group + 1) * descendants)
        if connected_components(adjacency[region, region], directed=False, return_labels=False) != 1:
            raise ValueError("A parent address groups disconnected cells")
    return adjacency


def source(car, level, sharp):
    if sharp:
        target = np.zeros(car.n)
        target[0] = 1.
    else:
        adjacency = parent_graph(car, level)
        seeds = np.arange(4**level)
        distance = shortest_path(adjacency, indices=seeds, unweighted=True, directed=False).min(axis=0)
        target = np.repeat(np.maximum(0., 1-2.**(-level)*distance/2), 12)
        target /= np.linalg.norm(target)
    return prepare(SimpleNamespace(car=car, active=car.seams), target)


def analyze(car, wave, meter, energy, seed, shots):
    cells = wave.reshape(car.cells, 12)
    local_reads = np.einsum("ci,pij,cj->cp", cells.conj(), meter.effects, cells).real
    exact = local_reads.reshape(20, car.cells//20, 12).sum(axis=1)
    measured, variances, counts = acquire(wave, meter, seed, shots)
    expected = exact[:, meter.trained] @ meter.forecast.T
    delta = exact[:, meter.held] - expected
    difference = measured[:, meter.held] - measured[:, meter.trained] @ meter.forecast.T
    errors = np.sqrt(variances[:, meter.held] + variances[:, meter.trained] @ (meter.forecast.T**2))
    centered = exact - exact.mean(axis=1, keepdims=True)
    angular_rms = np.sqrt(np.mean(centered**2))
    held_rms = np.sqrt(np.mean(delta**2))
    best_tensor_rms = np.sqrt(np.mean((exact @ meter.orthogonal_residual)**2))
    isotropic = exact[:, meter.held] - exact[:, meter.trained].mean(axis=1, keepdims=True)
    informative = angular_rms > 1e-12
    # Each block bound applies to the most discrepant direction in that block.
    # Cauchy–Schwarz controls the sum over disjoint blocks for one normalized load.
    bound = meter.slow_error + 2*meter.cross_norm*np.sqrt(energy/3) + meter.fast_norm*energy/3
    block_residual = np.max(abs(delta), axis=1).sum()
    if block_residual > bound + 1e-10:
        raise ValueError("Pooled readings exceed the kinetic-energy bound")
    resolved = errors > 1e-14
    max_standardized = max(abs(difference[resolved])/errors[resolved], default=0.)
    return dict(angular_signal_rms=float(angular_rms), held_error_rms=float(held_rms),
                isotropic_held_error_rms=float(np.sqrt(np.mean(isotropic**2))),
                relative_held_error=float(held_rms/angular_rms) if informative else None,
                best_all_direction_relative_error=float(best_tensor_rms/angular_rms) if informative else None,
                relative_target_pass=bool(held_rms <= .05*angular_rms) if informative else None,
                max_sampled_held_difference_in_standard_errors=float(max_standardized),
                sum_worst_held_absolute_error=float(block_residual), angular_bound=float(bound),
                sampled_held_differences=difference.round(11).tolist(),
                sampled_held_standard_errors=errors.round(11).tolist()), counts


def deposited_carrier(path, level):
    deposited = load_refinement(path, level=level)
    internal = [(a + 12*i, b + 12*i) for i in range(deposited.cells) for a, b in CELL_SEAMS]
    return Carrier(f"deposited level {level}", 12*deposited.cells, internal + list(deposited.external))


def run(holdout_board=None, holdout_level=None, precision=False):
    shots = 2**30 if precision else SHOTS
    seed_family = [903717] if precision else []
    meter = instruments()
    observations, records = [], {}
    for level in ([holdout_level] if holdout_board else range(3)):
        car = deposited_carrier(holdout_board, level) if holdout_board else federation(level)
        kinetic = candidate.generator(car) - 6*candidate.KAPPA*sparse.eye(car.n)
        difference = kinetic - car.laplacian()/2
        if np.max(abs(difference.data), initial=0.) > 1e-12:
            raise ValueError("The source generator differs from the registered one-load law")
        for sharp in ([False] if holdout_board else (False, True)):
            source_wave, gates = source(car, level, sharp)
            a, b = np.array(car.seams).T
            odd_chances = abs(source_wave[a] - source_wave[b])**2/2
            energy = float(odd_chances.sum())
            rng = np.random.default_rng([175112, *seed_family, level, int(sharp)])
            energy_counts = rng.binomial(shots, np.minimum(1., odd_chances))
            source_name = "sharp" if sharp else "smooth"
            records[f"{level}_{source_name}_antisymmetric_launch_counts"] = energy_counts
            np.testing.assert_allclose(energy, np.vdot(source_wave, kinetic @ source_wave), atol=1e-12)
            schedules = [("native", 1., 4. * 4**level)]
            if not sharp:
                schedules += [("input", 1., 0.), ("half_exchange", .5, 4. * 4**level)]
            for setting, (law, rate, dwell) in enumerate(schedules):
                wave = expm_multiply((-1j*rate*dwell)*kinetic, source_wave) if dwell else source_wave
                result, counts = analyze(car, wave, meter, energy, [*seed_family, level, int(sharp), setting], shots)
                key = f"{level}_{source_name}_{law}"
                records[key] = counts
                observations.append(dict(cells=car.cells, cells_per_parent=4**level, source=source_name,
                                         law=law, exposure=dwell, source_gates=gates,
                                         native_kinetic_energy=energy,
                                         sampled_native_kinetic_energy=float(energy_counts.sum()/shots),
                                         kinetic_energy_standard_error=float(np.sqrt((odd_chances*(1-odd_chances)).sum()/shots)),
                                         count_key=key, **result))
    result = dict(shots_per_setting=shots, calibration_ports=meter.trained, held_ports=meter.held,
                slow_block_error=meter.slow_error, slow_fast_residual_norm=meter.cross_norm,
                fast_residual_norm=meter.fast_norm, observations=observations)
    if holdout_board:
        result["deposited_board_sha256"] = DEPOSIT_SHA256[holdout_level]
    return result, records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--holdout-board", type=Path)
    parser.add_argument("--holdout-level", type=int, choices=(3, 4))
    parser.add_argument("--precision", action="store_true", help="Independent level-four count confirmation")
    arguments = parser.parse_args()
    if (arguments.holdout_board is None) != (arguments.holdout_level is None):
        parser.error("Supply the deposited file and its holdout level together")
    if arguments.precision and arguments.holdout_level != 4:
        parser.error("The precision acquisition uses the level-four screened holdout")
    result, counts = run(arguments.holdout_board, arguments.holdout_level, arguments.precision)
    arguments.output.mkdir(parents=True, exist_ok=True)
    suffix = f"_holdout_{arguments.holdout_level}" if arguments.holdout_board else ""
    if arguments.precision:
        suffix += "_precision"
    (arguments.output/f"results{suffix}.json").write_text(json.dumps(result, indent=2)+"\n")
    np.savez_compressed(arguments.output/f"counts{suffix}.npz", **counts)
