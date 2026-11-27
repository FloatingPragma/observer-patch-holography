"""Recover the native one-cell slow band from local transition records.

Run python -m model.evaluation.patch_spatial_readout to replay the receipt.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as contrast
from model.candidate_laws.interval_echo import candidate as echo
from model.executable_model import cell


def transitions(generator, duration):
    return np.abs(expm(-1j * duration * generator)) ** 2


def sign_control(generator, exposure):
    """Compare sharp port records and a supplied coherent preparation under H and -H*."""
    alternative = -generator.conj()
    error = max(float(np.max(np.abs(transitions(generator, n * exposure) -
                                    transitions(alternative, n * exposure)))) for n in (1, 2, 3))
    probe = np.zeros(len(generator), dtype=complex)
    probe[[0, 1]] = np.array([1, 1j]) / np.sqrt(2)
    probabilities = [abs((expm(-2j * exposure * h) @ probe)[0]) ** 2 for h in (generator, alternative)]
    if error > 1e-12 or abs(probabilities[0] - probabilities[1]) <= 1e-12:
        raise ArithmeticError("The sign control disagrees with the sharp-port or coherent prediction")
    return {
        "alternative": "Negative complex conjugate of the registered generator in the port basis.",
        "largest_sharp_port_error_below": 1e-12,
        "coherent_probe": "Supplied (|0> + i|1>)/sqrt(2) in the original port labels, followed by a port read at twice the calibration exposure.",
        "registered_probability_at_port_0": round(float(probabilities[0]), 12),
        "alternative_probability_at_port_0": round(float(probabilities[1]), 12),
    }


def slow_projector(table):
    eigenvalues, modes = np.linalg.eigh(table)
    level = max(eigenvalues[eigenvalues < 1 - 1e-9])
    selected = modes[:, np.abs(eigenvalues - level) < 1e-9]
    return selected @ selected.T, selected.shape[1]


def infer(calibration):
    """The reader sees anonymous square arrival tables and a uniform-hopping model class."""
    n = len(calibration)
    off_diagonal = np.sort(calibration[np.triu_indices(n, k=1)])
    gap = np.argmax(np.diff(off_diagonal))
    threshold = np.mean(off_diagonal[gap:gap + 2])
    neighbors = calibration > threshold
    np.fill_diagonal(neighbors, False)
    laplacian = np.diag(neighbors.sum(axis=1)) - neighbors
    # The source supplies a weak-exposure branch. The hop angle is inferred.
    def residual(angle):
        return np.mean(transitions(laplacian, angle)[neighbors] - calibration[neighbors])

    angle = brentq(residual, 0, .2, xtol=1e-14)
    projector, rank = slow_projector(calibration)
    scale = np.sqrt(projector.diagonal())
    gram = projector / np.outer(scale, scale)
    return angle, laplacian, projector, rank, gram


def measure():
    generators = [candidate.generator(cell()).toarray() for candidate in (contrast, echo)]
    if echo.native_law is not contrast.law or not np.array_equal(*generators):
        raise ValueError("The candidate generators differ on the one-load source.")
    h = generators[0]
    labels = np.random.default_rng(38611).permutation(len(h))

    def record(t):
        table = transitions(h, t)[np.ix_(labels, labels)]
        if not np.allclose(table.sum(axis=0), 1, rtol=0, atol=1e-12):
            raise ArithmeticError("The detector omits an outcome")
        return table

    exposure = 1 / 8
    calibration = record(exposure)
    angle, laplacian, projector, rank, gram = infer(calibration)
    held = []
    for multiplier in (2, 3):
        forecast = transitions(laplacian, multiplier * angle)
        observed = record(multiplier * exposure)
        independent, held_rank = slow_projector(observed)
        probability_error = float(np.max(np.abs(forecast - observed)))
        projector_error = float(np.linalg.norm(projector - independent))
        if max(probability_error, projector_error) > 1e-10:
            raise ArithmeticError("Fresh held records disagree with calibration")
        held.append({"relative_exposure": multiplier, "independently_read_rank": held_rank,
                     "max_probability_error_below": 1e-10, "projector_frobenius_error_below": 1e-10,
                     "records": observed.round(12).tolist()})
    return {
        "candidates": [contrast.NAME, echo.NAME], "cell_count": 1, "source_loads": 1,
        "instrument": "A one-cell observer prepares one sharp load at each port on separate trials, then reads its final port. All twelve exclusive outcomes are retained. Quiet has zero weight by load conservation. Tables are complete Born probabilities, with columns indexing preparations and rows indexing detections.",
        "inference_inputs": "Relabelled calibration probabilities, equal-hop graph model, weak exposure branch [0,0.2], and held ratios of supplied dwell settings. No graph, coordinates, Hamiltonian, spatial dimension or chosen local triplet are passed to the inverse.",
        "scope": "One accessible cell with an external exposure dial. This identifies three independent slow response components. Propagation of three spatial coordinates between cells and an independently ranged spacetime were not measured.",
        "fitted_edges": int(-np.triu(laplacian, k=1).sum()),
        "fitted_hop_angle": round(angle, 12), "calibration_rank": rank,
        "normalized_gram": gram.round(10).tolist(),
        "best_rank_two_relative_gram_error": round(float(
            np.linalg.norm(np.linalg.svd(gram, compute_uv=False)[2:]) / np.linalg.norm(gram)), 12),
        "calibration_records": calibration.round(12).tolist(), "held": held,
        "generator_sign_control": sign_control(h, exposure),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="write measurement to a new file")
    args = parser.parse_args()
    result = measure()
    if args.record:
        with args.record.open("x") as f:
            f.write(json.dumps(result, indent=2) + "\n")
    elif result != json.loads(Path(__file__).with_suffix(".json").read_text()):
        raise SystemExit("Native cell record replay differs from stored observations")
    print(f"Measured rank {result['calibration_rank']}, held ranks "
          f"{[x['independently_read_rank'] for x in result['held']]}")
