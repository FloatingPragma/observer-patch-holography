"""Test the integrated current response and its centre on the actual one-cell pair."""

from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, expm

from model.candidate_laws.interval_echo import currents
from model.executable_model import cell


STEPS = (512, 2048)
SHOTS = 8192
SEED = 216031
I3 = np.eye(3)
Z3 = np.zeros((3, 3))
I6 = np.eye(6)
NAMES = ("sham", "pure_weak_centre", "compensating_word", "diagonal_comparison")
LABELS = ((0, 0, 0), (0, 1, 0), (1, 1, 4), (1, 1, 1))


def coordinates(matrices):
    flat = np.asarray(matrices).reshape(len(matrices), -1)
    return np.column_stack((flat.real, flat.imag))


def ordered_word(addresses, coefficients, angle, steps):
    """Every half step applies precisely the registered pulse with negative signed strength."""
    half = [expm(0.5j * angle * x / steps * h) for x, h in zip(coefficients, addresses)]
    step = I6
    for u in half + half[::-1]:
        step = u @ step
    return np.linalg.matrix_power(step, steps)


def central_words(addresses, generators, coefficients, steps):
    angles = (2 * np.pi, 2 * np.pi, np.pi / 3)
    words = [ordered_word(addresses, c, angle, steps) for c, angle in zip(coefficients, angles)]
    errors = [np.linalg.norm(u - expm(1j * angle * g), 2)
              for u, angle, g in zip(words, angles, generators)]
    grid, selected, kernel = [], {}, []
    for a, b, c in product(range(3), range(2), range(6)):
        u = np.linalg.matrix_power(words[2], c) @ np.linalg.matrix_power(words[1], b) @ np.linalg.matrix_power(words[0], a)
        residue = (2 * a + c) % 6
        oracle = block_diag(np.exp(2j * np.pi * residue / 6) * I3, I3)
        grid.append(float(np.linalg.norm(u - oracle, 2)))
        if np.linalg.norm(u - I6, 2) < 0.1:
            kernel.append([a, b, c])
        if (a, b, c) in LABELS:
            selected[NAMES[LABELS.index((a, b, c))]] = u
    predicted_kernel = [list(k) for k in product(range(3), range(2), range(6))
                        if (2 * k[0] + k[2]) % 6 == 0]
    if kernel != predicted_kernel or max(grid) >= 0.05:
        raise ArithmeticError("The finite words did not resolve the proposed central actions.")
    half = ordered_word(addresses, coefficients[1], np.pi, steps)
    report = dict(base_word_operator_errors=errors, maximum_grid_operator_error=max(grid),
                  resolved_kernel=kernel, states_per_triplet=3, grid_size=len(grid),
                  weak_half_turn_distance_from_identity=float(np.linalg.norm(half - I6, 2)),
                  weak_half_turn_squared_error=float(np.linalg.norm(half @ half - I6, 2)))
    return selected, report


def run():
    car = cell()
    sector, source, h, _, bands, _, _, _, address_list, _ = currents._service(car, True, 1.0)
    addresses = np.array(address_list)
    first, second = addresses[:, :3, :3], addresses[:, 3:, 3:]
    form_error = max(np.max(abs(addresses[:, :3, 3:])), np.max(abs(addresses[:, 3:, :3])),
                     np.max(abs(second.real)), np.max(abs(second + second.transpose(0, 2, 1))),
                     np.max(abs(addresses - addresses.conj().transpose(0, 2, 1))))
    coords = coordinates(addresses)
    ranks = [int(np.linalg.matrix_rank(x, tol=1e-9))
             for x in (coordinates(first), coordinates(second), coords)]
    if ranks != [9, 3, 12] or form_error > 1e-10:
        raise ArithmeticError("The two triplet generators did not span the stated block algebras.")

    # These axes are selected by the native port-0 source and response.
    axis = bands[0].T @ source[:, 0]
    axis /= np.linalg.norm(axis)
    rotation = second[0] / np.linalg.norm(second[0], 2)
    generators = [block_diag(I3 / 3 - np.outer(axis, axis.conj()), Z3),
                  block_diag(Z3, rotation), block_diag(I3, Z3)]
    coefficients = np.linalg.lstsq(coords.T, coordinates(generators).T, rcond=1e-12)[0].T
    reconstruction = np.linalg.norm(coefficients @ coords - coordinates(generators), axis=1)
    if reconstruction.max() > 1e-10 or not np.allclose(np.linalg.eigvalsh(rotation), (-1, 0, 1)):
        raise ArithmeticError("The requested axes could not be calibrated on the address menu.")
    convergence = {}
    for steps in STEPS:
        selected, report = central_words(addresses, generators, coefficients, steps)
        convergence[str(steps)] = report

    basis = np.column_stack(bands)
    state = source[:, currents.SOURCE_PORT]
    kept = basis.T @ state
    spectator = state - basis @ kept
    free = expm(-16j * h)

    def readings(u):
        weights = abs(free @ (spectator + basis @ (u @ kept))) ** 2
        reading = np.r_[weights[sector.seam_states], weights[~sector.adjacent].sum()]
        if abs(reading.sum() - 1) > 1e-10:
            raise ArithmeticError("A terminal outcome was lost.")
        return reading / reading.sum()

    sham = readings(I6)
    ideal_diagonal = readings(block_diag(-I3, I3))
    event = ideal_diagonal - sham > 1e-10
    outcomes = {}
    for name, label, stream in zip(NAMES, LABELS, np.random.SeedSequence(SEED).spawn(len(NAMES))):
        probabilities = readings(selected[name])
        counts = np.random.default_rng(stream).multinomial(SHOTS, probabilities)
        event_count = int(counts[event].sum())
        baseline_count = outcomes["sham"]["event_count"] if name != "sham" else event_count
        measured = (event_count - baseline_count) / SHOTS
        pooled = (event_count + baseline_count) / (2 * SHOTS)
        margin = 4 * np.sqrt(2 * pooled * (1 - pooled) / SHOTS)
        outcomes[name] = dict(label=list(label), total_variation_from_sham=float(abs(probabilities - sham).sum() / 2),
                              probabilities=probabilities, counts=counts.tolist(),
                              event_count=event_count, event_probability=float(probabilities[event].sum()),
                              measured_increase_over_sham=measured, four_se_margin=margin,
                              resolved_increase=bool(measured > margin))
    return dict(source="Interval Echo, guarded native currents, full 66-state pair source",
                source_port=currents.SOURCE_PORT, free_duration=16,
                sphere_patch="one intact echosahedron", ranks=ranks,
                rank_tolerance=1e-9, address_basis_smallest_singular_value=np.linalg.svd(coords, compute_uv=False)[-1],
                maximum_block_form_error=form_error, axis_reconstruction_errors=reconstruction,
                axis_port_coefficients=coefficients, weak_axis_eigenvalues=np.linalg.eigvalsh(rotation),
                primitive_central_weights=[1, 0], cover_coordinate="primitive circle on the first triplet",
                order_of_central_labels=["SU3 cube roots", "SU2 signs", "primitive U1 sixth roots"],
                exact_pair_kernel="(omega^a I3, (-1)^b I2, omega^(-a)), a=0,1,2; b=0,1",
                convergence=convergence, shots_per_arm=SHOTS, random_seed=SEED,
                receiver_seams=[list(e) for e in car.seams], final_outcome="quiet",
                resolved_event_outcomes=np.flatnonzero(event).tolist(), terminal_measurements=outcomes)


def portable(x):
    if isinstance(x, dict):
        return {k: portable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [portable(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return round(float(x), 10) + 0.0
    return x


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the complete measurement to this JSON file.")
    args = parser.parse_args()
    receipt = portable(run())
    if args.output:
        args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    else:
        if receipt != json.loads(Path(__file__).with_name("result.json").read_text()):
            raise SystemExit("The archived current experiment differs.")
        print("Recovered the six-element native kernel on all 36 words.")
        for name, row in receipt["terminal_measurements"].items():
            print(name, row["event_count"], "/", SHOTS, row["total_variation_from_sham"])
