"""Neutral-side process tomography of the charged magnetic doublet."""

import argparse
import json
from math import atan2, log, sqrt
from pathlib import Path

import numpy as np
from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet.measure import (
    charged_source,
)
from model.experiments.interval_echo_lepton_register import measure as rails
from .aperture import compact_instrument
from model.experiments.interval_echo_gauge.arrival_ladder.measure import (
    aperture_map,
    layout,
)

HERE = Path(__file__).resolve().parent
PAULI = [
    np.array([[0, 1], [1, 0]], complex),
    np.array([[0, -1j], [1j, 0]], complex),
    np.diag([1.0, -1.0]),
]


def stations():
    board = federation(0)
    owner = layout()["source"] // 12
    ports = [12 * owner + i for i in range(12)]
    assert set(ports).isdisjoint(layout()["offline"])
    assert all((ports[a], ports[b]) in board.seams for a, b in rails.CELL_SEAMS)
    return owner, ports


def preparation(target, owner):
    """Reduce each desired output magnitude along one physical path, then set phases."""
    target = np.asarray(target, complex)
    np.testing.assert_allclose(np.vdot(target, target), 1, atol=1e-12)
    order = list(rails.path())
    unitary = np.eye(12, dtype=complex)
    for i, (a, b) in enumerate(zip(order, order[1:])):
        angle = atan2(np.linalg.norm(target[order[i + 1 :]]), abs(target[a]))
        rails.rotate(unitary, a, b, angle, 1, owner)
    for i, coefficient in enumerate(target):
        rails.phase(unitary, i, np.angle(coefficient), 1, owner)
    np.testing.assert_allclose(unitary[:, order[0]], target, atol=2e-12)
    return unitary


def six_sources(columns, owner):
    words = []
    for axis, operator in zip("XYZ", PAULI, strict=True):
        _, basis = np.linalg.eigh(operator)
        for label, eigenstate in [("+", basis[:, 1]), ("-", basis[:, 0])]:
            target = columns @ eigenstate
            pulse = preparation(target, owner)
            words.append((axis + label, pulse))
    return words


def run_channel(inputs, charged, neutral, operation, forget, rng):
    """Each setting has its own cohort; sum retained neutral and charged records."""
    records = []
    empirical = np.zeros((6, 6))
    predicted = np.zeros_like(empirical)
    source_rail = rails.path()[0]
    for j, (launch, source) in enumerate(inputs):
        initial = source[:, source_rail]
        state = neutral @ initial
        density = np.outer(state, state.conj())
        if forget:
            density = np.diag(density.diagonal())
        density = operation @ density @ operation.conj().T
        missed = charged @ initial
        charged_density = np.outer(missed, missed.conj())
        if forget:
            charged_density = np.diag(charged_density.diagonal())
        banks = []
        for i, (label, pulse) in enumerate(inputs):
            # The inverse state synthesis is a matched terminal vacancy test.
            # Resolve both rail alternatives by the common target and donor.
            n = (pulse.conj().T @ density @ pulse).diagonal().real
            c = (pulse.conj().T @ charged_density @ pulse).diagonal().real
            probabilities = np.clip(np.r_[n, c], 0, 1)
            np.testing.assert_allclose(probabilities.sum(), 1, atol=1e-11)
            counts = rails.count_bank(rng, probabilities)
            predicted[i, j] = n[source_rail]
            empirical[i, j] = counts[source_rail] / rails.SHOTS
            banks.append(dict(projector=label, terminal_counts=counts.tolist()))
        records.append(dict(preparation=launch, reads=banks))

    def infer(p):
        # Vector of expected I,X,Y,Z for each source from independent arm tests.
        even = (p[::2] + p[1::2]).mean(axis=0)
        obs = np.vstack((even, p[::2] - p[1::2]))
        return np.column_stack((obs.mean(axis=1), ((obs[:, ::2] - obs[:, 1::2]) / 2)))

    pred, seen = infer(predicted), infer(empirical)
    return dict(
        records=records,
        predicted_pauli_transfer=pred.tolist(),
        sampled_pauli_transfer=seen.tolist(),
        predicted_identity_entanglement_fidelity=float(np.trace(pred) / 4),
        sampled_identity_entanglement_fidelity=float(np.trace(seen) / 4),
    )


def energy_filtered_round_trip(owner, ports, columns, turns, gates, rng):
    """Retain both charge reads and repeat the same compact energy selection."""
    from model.experiments.interval_echo_lepton_register.pole_meter.measure import setup

    board, integer_projection, denominator, *_ = setup()
    assert board.n == federation(0).n and owner == ports[0] // 12
    p = integer_projection[np.ix_(ports, ports)] / denominator
    herald = compact_instrument(p, owner)
    accept, refuse = herald["accepted"], herald["refused"]
    weights, directions = np.linalg.eigh(columns.conj().T @ p @ columns)
    initial = columns @ directions
    bank = []
    for i, weight in enumerate(weights):
        state = preparation(initial[:, i], owner)[:, rails.path()[0]]
        admitted = accept @ state
        refused = refuse @ state
        for name, turn in [
            ("park", np.eye(12)),
            ("first", turns[0]),
            ("second", turns[1]),
        ]:
            first_miss = gates["charged_to_charged"] @ admitted
            neutral = turn @ gates["charged_to_neutral"] @ admitted
            residual = gates["neutral_to_neutral"] @ neutral
            returned = gates["neutral_to_charged"] @ neutral
            recovered = accept @ returned
            probabilities = np.array(
                [
                    np.linalg.norm(refused) ** 2,
                    np.vdot(first_miss, first_miss).real,
                    np.vdot(residual, residual).real,
                    np.vdot(recovered, recovered).real,
                    np.linalg.norm(refuse @ returned) ** 2,
                ]
            )
            np.testing.assert_allclose(probabilities.sum(), 1, atol=2e-12)
            bank.append(
                dict(
                    initial_principal_mode=i,
                    neutral_turn=name,
                    original_probabilities=probabilities.tolist(),
                    counts=rails.count_bank(rng, probabilities).tolist(),
                    recovered_given_source=float(probabilities[3] / weight),
                )
            )
    return dict(
        initial_compact_admission=weights.tolist(),
        common_flag_ports=list(map(int, herald["flag_ports"])),
        maximum_donor_residue=float(herald["max_donor_residue"]),
        flag_routes=herald["routes"],
        material_basis_rotations=int(herald["rotation_count"]),
        banks=bank,
        outcomes=[
            "source refused",
            "charged at first read",
            "neutral at second read",
            "returned charged in compact line",
            "returned charged outside compact line",
        ],
    )


def measure():
    owner, ports = stations()
    columns, projector, (r, s) = charged_source(owner)
    inputs = six_sources(columns, owner)
    driven = aperture_map(ports)
    sham = aperture_map(ports, angle=0)
    unit = np.eye(len(ports))
    rng = np.random.default_rng(95047)
    settings = {}
    for name, gates, operation, forget in [
        ("coherent", driven, unit, False),
        ("address_recorded", driven, unit, True),
        ("neutral_half_turn", driven, r, False),
        ("no_drive", sham, unit, False),
    ]:
        settings[name] = run_channel(
            inputs,
            gates["charged_to_charged"],
            gates["charged_to_neutral"],
            operation,
            forget,
            rng,
        )
    # Keep the prepared code fixed to isolate sensitivity to the receiver law.
    zero_contrast = aperture_map(ports, kappa=0)
    changed_law = run_channel(
        inputs, zero_contrast["charged_to_charged"], zero_contrast["charged_to_neutral"],
        unit, False, np.random.default_rng(95048),
    )
    # Theory check derived independently from the calibrated rotation on the code.
    logical = columns.conj().T @ r @ columns
    paulis = [np.eye(2), *PAULI]
    turn = np.array(
        [
            [np.trace(x @ logical @ y @ logical.conj().T).real / 2 for y in paulis]
            for x in paulis
        ]
    )
    np.testing.assert_allclose(
        settings["neutral_half_turn"]["predicted_pauli_transfer"], turn, atol=1e-11
    )
    dephased_fidelity = float(np.sum(projector.diagonal().real ** 2) / 4)
    np.testing.assert_allclose(
        settings["address_recorded"]["predicted_identity_entanglement_fidelity"],
        dephased_fidelity,
        atol=1e-12,
    )
    raw = rails.addresses(rails.CELL_SEAMS[0])
    a, _, plus, minus = raw
    family = rails.reader(raw, 2 / 9, 1, owner).conj().T[:, [a, minus, plus]]
    k = driven["charged_to_neutral"]
    family_effect = family.conj().T @ k.conj().T @ k @ family
    return dict(
        schema=1,
        owner=owner,
        source_rail=rails.path()[0],
        workstation_target=layout()["a"],
        donor=layout()["donor"],
        port_order=ports,
        shots_per_read=rails.SHOTS,
        complete_terminal_order=[
            "neutral target/donor=2,1 then 12 local vacancy addresses",
            "charged target/donor=1,2 then 12 local vacancy addresses",
        ],
        inferred_channel_order=["I", "X", "Y", "Z"],
        settings=settings,
        zero_contrast_receiver=changed_law,
        energy_selected_round_trip=energy_filtered_round_trip(
            owner, ports, columns, (r, s), driven, rng
        ),
        simultaneous_99_process_entry_radius=2
        * sqrt(log(2 * 6 * 6 * len(settings) / 0.01) / (2 * rails.SHOTS))
        + 1e-9,
        dephased_identity_fidelity_formula=dephased_fidelity,
        family_neutral_success=np.real(np.diag(family_effect)).tolist(),
        family_effect_distance_from_identity=float(
            np.max(abs(family_effect - np.eye(3)))
        ),
        calibrated_turn_matrix=turn.round(12).tolist(),
        maximum_order_sum_on_code=float(
            np.linalg.norm((r @ s + s @ r) @ k @ columns, ord=2)
        ),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    out = json.dumps(measure(), indent=2) + "\n"
    if p.parse_args().write:
        (HERE / "results.json").write_text(out)
    else:
        print(out, end="")


if __name__ == "__main__":
    main()
