"""One negative arrival, its route rest-line receipt, and the later outside-head comparison."""

import argparse
from hashlib import sha256
import json
from math import ceil, floor, sqrt
from pathlib import Path

import numpy as np

from model.executable_model import federation
from ......fan_recoil.fresh_kick.propagation.measure import CELLS, generator, prepare
from ......fan_recoil.fresh_kick.propagation.pole_clock.measure import spectral
from ......fan_recoil.self_start.measure import DRAW_PAD, bound
from ......replay import compare
from ..local_read import PAULI
from ..measure import inputs
from ..spectrum import lines, readers, spectral_filter
from ..weak_probe import angular_input
from ..wheel import LAND, RIM, preparation, sides, wheel_basis
from .instrument import SHARED, exterior
from .pole_handoff.check import collector_replay, primary_layout
from .record_banks import current_bank, fixed_witness, frame, sample, witness_bank


HERE = Path(__file__).resolve().parent
GUIDE = HERE.parents[4] / "fan_recoil" / "fresh_kick" / "propagation"
SHOTS = 2**48


def collect_columns(target):
    """Invert the registered two-port rotations and phase plates on all inputs."""
    size = len(target)
    columns = np.eye(size, dtype=complex)
    tails = np.sqrt(np.cumsum(target[::-1] ** 2))[::-1]
    for j in range(size - 2, -1, -1):
        cosine, sine = target[j] / tails[j], tails[j + 1] / tails[j]
        # For a positive real target the forward primitive includes diag(1,-i)
        # after exp(+i*theta*X). Use its adjoint in reverse contact order.
        two = np.array([[cosine, sine], [-1j * sine, 1j * cosine]])
        columns[j : j + 2] = two @ columns[j : j + 2]
    np.testing.assert_allclose(
        columns.conj().T @ columns, np.eye(size), atol=4e-14, rtol=0
    )
    np.testing.assert_allclose(columns[0], target, atol=3e-14, rtol=0)
    return columns


def route_read(rng):
    saved_route = (GUIDE / "results.json").read_bytes()
    saved_clock = (GUIDE / "pole_clock" / "records.json").read_bytes()
    route = json.loads(saved_route)
    pole = json.loads(saved_clock)
    ports = route["ports"]
    board = federation(2)
    assert len(ports) == 3 * CELLS and ports[:4] == [1045, 1052, 1047, 1032]
    degrees = np.asarray(board.degree)[ports]
    np.testing.assert_array_equal(degrees, np.tile([6, 5, 6], CELLS))
    for a, b in zip(ports[:-1], ports[1:], strict=True):
        assert b in board.neighbours[a]
    h = generator(degrees)
    energies, eig, _ = spectral(h)
    e0 = 10 - 2 * sqrt(17)
    f = np.tile([1, 1 - e0 / 4, 1], CELLS)
    f /= np.linalg.norm(f)
    assert np.linalg.norm(h @ f - e0 * f) < 2e-13
    np.testing.assert_allclose(prepare(f.astype(complex)), f, atol=2e-13, rtol=0)
    collector = collect_columns(f)
    filter_data = pole["filter"]
    stages = filter_data["stages"]
    size = filter_data["coherent_wait_appointments_per_stage"]
    step = filter_data["appointment_spacing"]
    center = filter_data["center_measured_in_appointment_units"]

    def response(energy):
        x = (energy - center) * step / (2 * np.pi)
        return (
            np.exp(-1j * np.pi * (size - 1) * x) * np.sinc(size * x) / np.sinc(x)
        ) ** stages

    # The coherent acceptance, not just its squared modulus, precedes collection.
    response_all = response(energies)
    filtered = (eig * response_all) @ eig.T
    actual_row = collector[0] @ filtered
    expected_row = response(e0) * f
    # The protected analytic identity in spectral coordinates is exact. This
    # direct matrix replay allows the phase error of a long coherent wait.
    assert np.linalg.norm(actual_row - expected_row) < 3e-10
    ground_gain = float(abs(response(e0)) ** 2)
    returned = ground_gain * float(f[0] ** 2)
    off_pole = float(abs(eig[0, 1:]) ** 2 @ abs(response_all[1:]) ** 2)
    accepted = returned + off_pole
    assert 1 - filter_data["pole_attenuation_upper"] <= ground_gain <= 1 + 2e-15
    assert off_pole <= filter_data["off_pole_response_upper"]
    # A spatially odd input on the two end ports has zero overlap with the
    # even ground line. This independent preparation uses the same contacts.
    odd = np.zeros(len(f), dtype=complex)
    odd[[0, -1]] = np.array([1, -1]) / sqrt(2)
    np.testing.assert_allclose(prepare(odd), odd, atol=3e-13, rtol=0)
    assert abs(actual_row @ odd) < 2e-13
    odd_filter_mass = float(np.linalg.norm(filtered @ odd) ** 2)
    # Removing every coherent appointment leaves precisely the same rank-one
    # spatial projector, now with unit efficiency on the rest eigenmode.
    only = collector[0].conj()[:, None] @ collector[0][None, :]
    np.testing.assert_allclose(only, np.outer(f, f), atol=2e-14, rtol=0)
    np.testing.assert_allclose(h @ only, e0 * only, atol=3e-14, rtol=0)
    bare_rng = np.random.RandomState(10451753)
    bare = {}
    for key, vector in (
        ("sharp_1045", np.eye(1, len(f), 0)[0]),
        ("native_swap_1052", np.eye(1, len(f), 1)[0]),
        ("rest_packet", f),
        ("odd_two_ends", odd),
    ):
        prob = float(abs(collector[0] @ vector) ** 2)
        assert abs(prob - abs(np.vdot(f, vector)) ** 2) < 4e-14
        bare[key] = {
            "return_probability": prob,
            "return_and_refusal_counts": sample(
                [prob, 1 - prob], SHOTS, bare_rng
            ).tolist(),
        }
    calibrations = {}
    for name, pf, pc in (
        ("rest_packet", ground_gain, ground_gain),
        ("sharp_arrival", accepted, returned),
        ("odd_two_ends", odd_filter_mass, 0),
    ):
        calibrations[name] = sample([pc, pf - pc, 1 - pf], SHOTS, rng).tolist()
    return (
        (accepted, returned),
        ports,
        {
            "guide_receipt_sha256": sha256(saved_route).hexdigest(),
            "calibration_sha256": sha256(saved_clock).hexdigest(),
            "route_words": pole["signed_two_word_gate"],
            "route_ports": len(ports),
            "rest_frequency_in_neutral_clock_units": pole[
                "rest_frequency_in_neutral_clock_units"
            ],
            "inverse_two_contact_pulses": len(ports) - 1,
            "collector_coefficient_squared_given_sharp_capture": returned,
            "certified_coefficient_squared_interval": [
                float((1 - filter_data["pole_attenuation_upper"]) * f[0] ** 2 - 1e-15),
                float(f[0] ** 2 + 1e-15),
            ],
            "sharp_filter_acceptance_given_capture": accepted,
            "sharp_filtered_excited_weight_given_capture": off_pole,
            "off_pole_acceptance_after_collection": 0,
            "collector_without_filter": {
                "coherent_filter_stages": 0,
                "originals_per_calibration_input": SHOTS,
                "independent_inputs": bare,
                "off_pole_acceptance": 0,
                "full_filter_vs_bypass_maximum_effect_difference": filter_data[
                    "pole_attenuation_upper"
                ],
            },
            "calibration_shots_given_input": SHOTS,
            "calibration_bins": [
                "origin_returned",
                "collector_refused",
                "filter_refused",
            ],
            "calibration_counts": calibrations,
        },
    )


def run():
    rng = np.random.RandomState(104510521)
    collector, ports, receipt = route_read(rng)
    route, basis = preparation(), wheel_basis()
    handoff = primary_layout()
    joined = collector_replay(handoff, ports, receipt["route_words"])
    assert joined["returned_collector_word_matches_before_handoff"]
    assert handoff["wheel_input_words"] == route["wheel_input_words"]
    # The restored origin has zero reversed-prefix heads. Both route words
    # pass the handoff without any spectator-dependent phase.
    effects, short, long, currents, _ = exterior(basis, route)
    h, _, _, rows = lines(basis, 0.5)
    index = next(i for i, row in enumerate(rows) if (1, 0) in row[2])
    selected = spectral_filter(h, rows, index)
    v, eps = angular_input(basis, h)
    u = v[:, 0] * v[0, 0].conj() / abs(v[0, 0])
    v = np.column_stack((u, u.conj()))
    np.testing.assert_allclose(v @ v.conj().T, selected, rtol=0, atol=2e-13)
    k = v[[basis.index((LAND, word)) for word in route["wheel_input_words"]]].conj().T
    transfer = np.kron(k, np.eye(2))
    side_read = readers(basis, currents)[sides(RIM).index(SHARED)]
    ops, design, tom = frame(effects, v, eps)
    w, coeff, guard = fixed_witness(
        design, transfer, ops, tom["calibration_operator_norm_guard"]
    )
    saved = json.loads((HERE / "records.json").read_text())
    np.testing.assert_allclose(
        coeff, saved["witness_coefficients_outside_by_reference"], atol=3e-11, rtol=0
    )
    direct_w, direct_coeff, direct_guard = fixed_witness(np.eye(4), np.eye(4), PAULI, 0)
    cases = []
    for time, bank, pad, t in inputs()[1]:

        def read(matrix, options, c, witness, error, randomized):
            return witness_bank(
                t,
                matrix,
                options,
                c,
                witness,
                pad,
                error,
                randomized,
                rng,
                collector=collector,
                shots=SHOTS,
                comparisons=6,
                receiver="face" if options is PAULI else "outside",
            )

        raw = read(np.eye(4), PAULI, direct_coeff, direct_w, direct_guard, False)
        outer = read(transfer, ops, coeff, w, guard, False)
        # Enrollment precedes either analyzer choice. Use the first setting's
        # acceptance counts directly, including every original refusal.
        passed = sum(outer["settings"][0]["counts"][:4])
        qlo, qhi = receipt["certified_coefficient_squared_interval"]
        acceptance_pad = qhi * (4 * pad + 4e-10) + (qhi - qlo) * t[0, 0] + 4 * DRAW_PAD
        acceptance = bound(passed, SHOTS, acceptance_pad)
        ext_current = current_bank(
            bank,
            v @ k,
            short,
            long,
            side_read,
            pad,
            rng,
            collector=collector,
            shots=2**46,
        )
        case = {
            "arrival_appointment": time,
            "sharp_capture_per_original": float(t[0, 0]),
            "collector_returns_per_original": float(t[0, 0] * collector[1]),
            "joint_heralds_per_million_99pct_first_exterior_setting": [
                floor(1e10 * acceptance[0]) / 1e4,
                ceil(1e10 * acceptance[1]) / 1e4,
            ],
            "face_reference_after_collection": raw,
            "wheel_exterior_reference_after_collection": outer,
            "wheel_current_after_collection": ext_current,
        }
        if time == 30.8125:
            case["face_scrambled_reference"] = read(
                np.eye(4), PAULI, direct_coeff, direct_w, direct_guard, True
            )
            case["exterior_scrambled_reference"] = read(
                transfer, ops, coeff, w, guard, True
            )
            assert all(
                z["simultaneous_99pct"][0] > 0
                for z in (
                    case["face_scrambled_reference"],
                    case["exterior_scrambled_reference"],
                )
            )
            assert (
                raw["simultaneous_99pct"][1] < 0 and outer["simultaneous_99pct"][1] < 0
            )
            assert (
                ext_current["reports"]["plus"]["outside_current_simultaneous_99pct"][1]
                < 0
            )
            assert (
                ext_current["reports"]["minus"]["outside_current_simultaneous_99pct"][0]
                > 0
            )
        else:
            assert (
                raw["simultaneous_99pct"][0] > 0 and outer["simultaneous_99pct"][0] > 0
            )
        cases.append(case)
    return {
        "route_collector": receipt,
        "handoff": {
            "commissioned_background": handoff["commissioned_background"],
            "native_address_schedule": [
                list(op) for op in handoff["native_address_schedule"]
            ],
            "exactly_replayed_route_ports_for_both_signed_words": joined[
                "physical_route_ports_certified"
            ],
            "recovered_wheel_input_words_in_positive_chart": list(
                handoff["wheel_input_words"]
            ),
            "join_receipt_sha256": sha256(
                (HERE / "pole_handoff" / "results.json").read_bytes()
            ).hexdigest(),
        },
        "wheel_energy_filter_rejection_waits": len(rows) - 1,
        "appointments": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    rec = run()
    path = HERE / "collected.json"
    if args.write:
        path.write_text(json.dumps(rec, indent=2, allow_nan=False) + "\n")
    else:
        compare(rec, json.loads(path.read_text()), "pole_and_exterior")
    for case in rec["appointments"]:
        print(
            case["arrival_appointment"],
            case["wheel_exterior_reference_after_collection"]["simultaneous_99pct"],
        )


if __name__ == "__main__":
    main()
