"""Directional marginal modular readings on a serviced cell of the federation."""

import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation, rotations
from model.experiments.interval_echo_einstein_station import cap
from model.experiments.interval_echo_space.shell_source import prepare

HERE = Path(__file__).parent
SHOTS = 2**18


def common_reference():
    modes = cap.exact_stationary_source(1, cap.standard_count())["mode_probabilities"]
    odds = np.array(modes[1:3]) / modes[0]
    native = cap.NativeFilter()
    sweeps = 4
    filters = [
        (cap.seam_indices(edge), native.success(odds, 1 / (4 * sweeps)))
        for edge in cell().seams
    ]
    parent, _, acceptance, _, _ = cap.prepared_charge(filters, sweeps)
    mixed = np.zeros_like(parent)
    for rotation in rotations():
        indices = [
            cap.WORD_INDEX[tuple(sorted(rotation[p] for p in word))]
            for word in cap.WORDS
        ]
        mixed[np.ix_(indices, indices)] += parent / len(rotations())
    return mixed, parent, acceptance


def regional_operators(parent):
    stars = [tuple(sorted((p, *cell().neighbours[p]))) for p in range(12)]
    effects, floors = [], []
    for ports in stars:
        region = cap.CapRestriction(cap.WORDS, ports)
        reading, eigenvalues = region.charge(parent)
        effects.append(region.lift(reading)[1:13, 1:13])
        floors.append(min(eigenvalues))
    return stars, np.array(effects), min(floors)


def read_star(psi, ports, operator, rng):
    """Estimate a vacuum-centred real observable. All complementary clicks count zero."""
    ids = list(ports)
    diagonal = np.diag(operator)[ids]
    probs = abs(psi[ids]) ** 2
    p = np.r_[probs, max(0.0, 1 - probs.sum())]
    counts = rng.multinomial(SHOTS, p / p.sum())
    record = hashlib.sha256(np.asarray(counts, dtype="<i8").tobytes())
    total = float(diagonal @ counts[:-1] / SHOTS)
    variance = float((diagonal**2 @ probs) - (diagonal @ probs) ** 2)
    for i, j in combinations(ports, 2):
        plus, minus = (abs(psi[i] + sign * psi[j]) ** 2 / 2 for sign in (1, -1))
        chances = np.array([plus, minus, max(0.0, 1 - plus - minus)])
        outcomes = rng.multinomial(SHOTS, chances / chances.sum())
        record.update(np.asarray(outcomes, dtype="<i8").tobytes())
        weight = operator[i, j]
        total += weight * (outcomes[0] - outcomes[1]) / SHOTS
        variance += weight**2 * (plus + minus - (plus - minus) ** 2)
    return total, variance / SHOTS, record.hexdigest()


def run():
    parent, untwirled, acceptance = common_reference()
    stars, effects, floor = regional_operators(parent)
    n = cap.chord_coordinates()
    basis = np.column_stack((np.ones(12) / np.sqrt(12), n / 2))
    np.testing.assert_allclose(basis.T @ basis, np.eye(4), atol=3e-14)
    x, y, z = n.T
    design = np.column_stack(
        (np.ones(12), x, y, z, x * x - z * z, y * y - z * z, x * y, x * z, y * z)
    )
    trained = json.loads((HERE / "results.json").read_text())["angular_design"][
        "best_nine_example_native_ports"
    ]
    held = sorted(set(range(12)) - set(trained))
    forecast = design[held] @ np.linalg.inv(design[trained])
    compressed = basis.T @ effects @ basis
    # Stabilizer and transitivity give the same four axial coefficients at each pole.
    block = compressed[0]
    scalar = block[0, 0]
    mixed = block[0, 1:] @ n[0]
    along = n[0] @ block[1:, 1:] @ n[0]
    transverse = (np.trace(block[1:, 1:]) - along) / 2
    anisotropy = along - transverse
    predicted = np.zeros_like(compressed)
    predicted[:, 0, 0] = scalar
    predicted[:, 0, 1:] = predicted[:, 1:, 0] = mixed * n
    predicted[:, 1:, 1:] = transverse * np.eye(3) + anisotropy * np.einsum(
        "pi,pj->pij", n, n
    )
    # On trace-zero variations a state-independent scalar is unobservable.
    variations = (
        compressed
        - np.trace(compressed, axis1=1, axis2=2)[:, None, None] * np.eye(4) / 4
    )

    def rank(operators):
        return int(np.linalg.matrix_rank(operators.reshape(12, -1), tol=1e-9))

    def closure(operators):
        differences = operators[held] - np.einsum(
            "ij,jkl->ikl", forecast, operators[trained]
        )
        return max(float(np.linalg.norm(block, 2)) for block in differences)

    singular_values = np.linalg.svd(variations.reshape(12, -1), compute_uv=False)
    geometry = federation(0)
    board = SimpleNamespace(car=geometry, active=geometry.seams)
    pilots = [
        (basis[:, 0] + (n @ n[0]) / 2) / np.sqrt(2),
        (basis[:, 0] + 1j * (n @ n[5]) / 2) / np.sqrt(2),
        np.eye(12)[0],
    ]
    rows = []
    for pilot, target in enumerate(pilots):
        full = np.zeros(geometry.n, complex)
        full[:12] = target
        prepared, gates = prepare(board, full)
        for setting, time in enumerate((0.75, 2.0)):
            # External contacts and nonuniform degree pulses are parked in flight.
            # The internal seam sum and cell generator differ only by a scalar.
            psi = expm(-1j * time * candidate.response(cell())) @ prepared[:12]
            detections = [
                read_star(
                    psi, star, op, np.random.default_rng([39107, pilot, setting, pole])
                )
                for pole, (star, op) in enumerate(zip(stars, effects, strict=True))
            ]
            readings, variances = np.array([row[:2] for row in detections]).T
            exact = np.einsum("i,pij,j->p", psi.conj(), effects, psi).real
            expected = forecast @ readings[trained]
            errors = np.sqrt(variances[held] + forecast**2 @ variances[trained])
            rows.append(
                dict(
                    source=(
                        "constant_plus_real_triplet",
                        "constant_plus_imaginary_triplet",
                        "sharp_port",
                    )[pilot],
                    exposure=time,
                    launch=gates,
                    measured=readings.tolist(),
                    held_predictions=expected.tolist(),
                    complete_count_hashes=[row[2] for row in detections],
                    held_difference_standard_errors=errors.tolist(),
                    max_held_standardized_error=float(
                        max(abs(readings[held] - expected) / errors)
                    ),
                    max_abs_exact_held_miss=float(
                        max(abs(exact[held] - forecast @ exact[trained]))
                    ),
                )
            )
    _, unrotated, _ = regional_operators(untwirled)
    return dict(
        cells=geometry.cells,
        reference_sweeps=4,
        reference_rotations=len(rotations()),
        reference_acceptance=float(acceptance),
        min_star_density_eigenvalue=float(floor),
        shots_per_setting=SHOTS,
        settings_per_star=16,
        calibration_ports=trained,
        held_ports=held,
        slow_state_variation_rank=rank(variations),
        axial_coefficients=dict(
            scalar=float(scalar),
            mixed=float(mixed),
            transverse=float(transverse),
            anisotropy=float(anisotropy),
        ),
        scalar_variation_contrast=float(scalar - transverse - anisotropy / 3),
        max_axial_form_error=float(np.max(abs(compressed - predicted))),
        ninth_variation_singular_value=float(singular_values[8]),
        pure_triplet_observable_rank=rank(compressed[:, 1:, 1:]),
        full_one_load_observable_rank=rank(effects),
        slow_sector_held_operator_error=closure(compressed),
        full_sector_held_operator_error=closure(effects),
        single_service_order_slow_sector_held_error=closure(
            basis.T @ unrotated @ basis
        ),
        readings=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    encoded = json.dumps(result, indent=2) + "\n"
    if args.write:
        (HERE / "modular.json").write_text(encoded)
    else:
        print(encoded, end="")
