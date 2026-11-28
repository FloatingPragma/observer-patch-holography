"""Replay the outside current and source-reference reads on the captured material."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from ......replay import compare
from ..apron import instrument as outside_hamiltonian
from ..measure import READS, inputs
from ..spectrum import lines, readers, spectral_filter
from ..weak_probe import angular_input
from ..wheel import LAND, preparation, sides, RIM, wheel_basis
from .instrument import SHARED, exterior
from .record_banks import current_bank, fixed_witness, frame, witness_bank


HERE = Path(__file__).resolve().parent


def run():
    route, basis = preparation(), wheel_basis()
    effects, short, long, currents, setup = exterior(basis, route)
    h, _, _, energies = lines(basis, 0.5)
    at = next(i for i, row in enumerate(energies) if (1, 0) in row[2])
    passed = spectral_filter(h, energies, at)
    v, error = angular_input(basis, h)
    # Anchor an otherwise arbitrary eigenvector phase, for portable coefficients.
    assert abs(v[0, 0]) > 0.1
    u = v[:, 0] * v[0, 0].conj() / abs(v[0, 0])
    v = np.column_stack((u, u.conj()))
    np.testing.assert_allclose(v @ v.conj().T, passed, atol=2e-13, rtol=0)
    k = v[[basis.index((LAND, word)) for word in route["wheel_input_words"]]].conj().T
    embed = np.eye(len(basis))[
        :, [basis.index((LAND, w)) for w in route["wheel_input_words"]]
    ]
    side = sides(RIM).index(SHARED)
    for j in currents:
        np.testing.assert_allclose(
            v.conj().T @ (currents[side] - j) @ v, 0, atol=3e-14, rtol=0
        )
    side_read = readers(basis, currents)[side]
    transfer = np.kron(k, np.eye(2))
    ops, matrix, apparatus = frame(effects, v, error)
    w, coeff, guard = fixed_witness(
        matrix, transfer, ops, apparatus["calibration_operator_norm_guard"]
    )
    rng = np.random.RandomState(10531045)
    axis, appointments = inputs()
    cases = []
    for wait, bank, pad, moments in appointments:
        selected = current_bank(bank, v @ k, short, long, side_read, pad, rng)
        if wait == 30.8125:
            assert (
                selected["reports"]["plus"]["outside_current_simultaneous_99pct"][1] < 0
            )
            assert (
                selected["reports"]["minus"]["outside_current_simultaneous_99pct"][0]
                > 0
            )
            assert (
                selected["reports"]["unread"]["outside_current_simultaneous_99pct"][1]
                < 0
            )
        witness = witness_bank(moments, transfer, ops, coeff, w, pad, guard, False, rng)
        if wait == 27.0625:
            assert witness["simultaneous_99pct"][0] > 0
        else:
            assert witness["simultaneous_99pct"][1] < 0
        cases.append(
            {
                "arrival_appointment": wait,
                "energy_heralded": selected,
                "without_energy_selection": current_bank(
                    bank, embed, short, long, side_read, pad, rng
                ),
                "source_and_outside": witness,
            }
        )
        if wait == 30.8125:
            mixed = witness_bank(
                moments, transfer, ops, coeff, w, pad, guard, True, rng
            )
            assert mixed["simultaneous_99pct"][0] > 0
            cases[-1]["reference_randomized"] = mixed

    # At unchanged elapsed time, halving only the adjoining face rate gives
    # a quarter-cycle in each enabled face block. The outside flip probability
    # and the two-setting response both drop by one half on every word.
    closed, ports, changed, _ = outside_hamiltonian(basis, route, 0)
    half, *_ = outside_hamiltonian(basis, route, 0.25)
    involution = np.eye(len(half)) - (half - closed) / 0.25
    quarter = (
        (np.eye(len(half)) + involution) - 1j * (np.eye(len(half)) - involution)
    ) / 2
    active = quarter[changed][:, ports]
    np.testing.assert_array_equal(active.conj().T @ active, effects["bare"] / 2)
    # On a closed face the outside record always stays put. The fair setting
    # then has agreement probability 1/2, for every wheel preparation.
    np.testing.assert_array_equal(changed[ports], np.zeros(len(ports), dtype=bool))
    return {
        "source_word_reads_sha256": sha256(READS.read_bytes()).hexdigest(),
        "written_words_and_waits": setup,
        "reference_report_axis": axis.tolist(),
        "current_scale_times_twice_agreement_probability_minus_one": 4,
        "wheel_energy_filter_rejection_waits": len(energies) - 1,
        "half_face_rate_same_wait_contrast_factor": 0.5,
        "tomography": apparatus,
        "witness_coefficients_outside_by_reference": coeff.tolist(),
        "appointments": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    receipt = HERE / "records.json"
    if args.write:
        receipt.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare(result, json.loads(receipt.read_text()), "exterior_meter")
    for case in result["appointments"]:
        print(
            case["arrival_appointment"],
            case["source_and_outside"]["simultaneous_99pct"],
        )


if __name__ == "__main__":
    main()
