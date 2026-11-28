"""Replay joint signed records on the finite three-arm threshold preparations."""

import argparse
from hashlib import sha256
from itertools import pairwise, product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta as beta_distribution

from ... import apparatus, measure as source
from ...recoil_search import threshold
from ..measure import conjugation

HERE = Path(__file__).resolve().parent
THRESHOLD = threshold.HERE / "threshold.json"
SEED = 740011
BOUNDARY_STARTS = 2**46
CUT_STARTS = 2**41
TAIL = 0.01 / 128


def calibrated_echo():
    law = apparatus.LAW
    echo = law.echo()
    labels = {}
    for raw, (a, b) in enumerate(law.SRC):
        outputs = np.flatnonzero(echo[:, raw])
        assert len(outputs) == 2 and min(outputs) >= law.N_SRC
        observed = set()
        for j in outputs:
            mem, stock = divmod(int(j), law.N_SRC)
            _, contrast = law.LABELS[mem - 1]
            observed.add((sum(law.SRC[stock]) - 2, contrast))
        assert observed == {(a + b - 2, a - b)}
        np.testing.assert_allclose(np.sum(abs(echo[outputs, raw]) ** 2), 1, atol=1e-14)
        labels[a, b] = observed.pop()
    assert len(set(labels.values())) == len(law.SRC)
    assert labels[0, 0] == (-2, 0) and labels[2, 2] == (2, 0)
    return labels


def raw_record(labels, seam, address, charge):
    a, b = sorted(seam)
    return labels[1 + charge * (address == a), 1 + charge * (address == b)]


def idle_boundary(board, fixture):
    owner = fixture["tip"] // 12
    boundary = [edge for edge in board.external if any(p // 12 == owner for p in edge)]
    assert len(boundary) == 3
    boundary_ports = set(sum((list(edge) for edge in boundary), []))
    assert fixture["tip"] not in boundary_ports and not boundary_ports.intersection(
        fixture["donor"]
    )
    # Every evolving written-memory edge is internal to one owner. Cross-owner
    # collars are prepared blank or held in the zero-write block.
    assert all(a // 12 == b // 12 for a, b in fixture["marked"])
    np.testing.assert_array_equal(apparatus.blocks()[0, 1], np.zeros((18, 18)))
    law = apparatus.LAW
    blank = np.eye(law.DIM)[:, law.SRC.index((1, 1))]
    idle_write = law.echo() @ blank
    for term in law.local_term_parts():
        np.testing.assert_allclose(term @ blank, 0, atol=1e-14)
        np.testing.assert_allclose(term @ idle_write, 0, atol=1e-14)
    np.testing.assert_allclose(law.echo() @ idle_write, blank, atol=1e-14)
    spoke = min(
        tuple(sorted(edge)) for edge in fixture["blank"] if fixture["tip"] in edge
    )
    assert all(p // 12 == owner for p in spoke) and not boundary_ports.intersection(
        spoke
    )
    return owner, boundary, boundary_ports, spoke


def position_probabilities(weighted):
    return (abs(weighted.reshape(-1, 64)) ** 2).sum(axis=1) / 3


def draw(probabilities, starts, rng):
    assert min(probabilities) >= -1e-13 and abs(sum(probabilities) - 1.0) < 1e-12
    valid = np.maximum(probabilities, 0)
    return rng.multinomial(starts, valid / valid.sum()).tolist()


def bank(probabilities, starts, rng):
    outcomes = {name: float(p) for name, p in probabilities.items() if p > 0}
    outcomes["clock refusal"] = (1 / 3) - sum(outcomes.values())
    outcomes["source refusal"] = 2 / 3
    counts = draw(list(outcomes.values()), starts, rng)
    return dict(
        computed_probabilities=outcomes,
        synthetic_counts=dict(zip(outcomes, counts, strict=True)),
    )


def bin_interval(success, total):
    if not total:
        return [0.0, 1.0]
    low = beta_distribution.ppf(TAIL, success, total - success + 1) if success else 0.0
    high = (
        beta_distribution.isf(TAIL, success + 1, total - success)
        if success < total
        else 1.0
    )
    return [float(low), float(high)]


def boundary_bank(labels, fixture, spoke, boundary_ports, ports, weights, charge, rng):
    record = {}
    for port, mass in zip(ports, weights, strict=True):
        deposit = raw_record(labels, spoke, port, charge)
        signature = "changed" if port in boundary_ports else "idle"
        key = f"{deposit[0]},{deposit[1]}/{signature}"
        record[key] = record.get(key, 0.0) + mass
    deposit = raw_record(labels, spoke, fixture["tip"], charge)
    key = f"{deposit[0]},{deposit[1]}"
    arrival = record.get(f"{key}/idle", 0.0)
    assert arrival > 0 and record.get(f"{key}/changed", 0.0) == 0
    data = bank(record, BOUNDARY_STARTS, rng)
    observed = data["synthetic_counts"].get(f"{key}/idle", 0)
    data.update(
        admitted_tip_reading=list(deposit),
        signed_tip_arrival_per_original_start=float(arrival),
        signed_tip_arrivals=observed,
        boundary_idle_given_signed_tip=1.0,
        boundary_change_given_signed_tip_upper99=bin_interval(0, observed)[1],
    )
    return data


def arm_bank(labels, cut, ports, weights, charge, rng):
    responses = [(0, 0), (charge, charge), (charge, -charge)]
    joint = np.zeros((3,) * len(cut))
    for port, mass in zip(ports, weights, strict=True):
        measured = [raw_record(labels, seam, port, charge) for seam in cut]
        joint[tuple(responses.index(record) for record in measured)] += mass
    admission = float(joint.sum())
    law = joint / admission
    marginal = np.stack(
        [law.sum(axis=tuple(j for j in range(3) if j != arm)) for arm in range(3)]
    )
    independent = np.einsum("i,j,k->ijk", *marginal)
    patterns = list(product(range(3), repeat=3))
    multiple = np.array([sum(x != 0 for x in p) > 1 for p in patterns])
    multi = float(independent.ravel()[multiple].sum())
    assert np.all(law.ravel()[multiple] == 0) and multi > 0
    for arm in range(3):
        axes = tuple(j for j in range(3) if j != arm)
        np.testing.assert_allclose(
            law.sum(axis=axes), independent.sum(axis=axes), atol=5e-16
        )
    keys = [
        "/".join(f"{responses[i][0]},{responses[i][1]}" for i in pattern)
        for pattern in patterns
    ]
    native = bank(dict(zip(keys, joint.ravel(), strict=True)), CUT_STARTS, rng)
    control = bank(
        dict(zip(keys, admission * independent.ravel(), strict=True)), CUT_STARTS, rng
    )
    return dict(
        admission=admission,
        local_conditional_probabilities=marginal.tolist(),
        native=native,
        independent_local_control=control,
        native_multiple_deposits_given_clock=0.0,
        independent_multiple_deposits_given_clock=multi,
        native_vs_independent_conditional_total_variation=float(
            np.sum(abs(law - independent)) / 2
        ),
        native_vs_independent_original_start_total_variation=float(
            admission * np.sum(abs(law - independent)) / 2
        ),
    )


def spectra(length_case, deposited):
    _, ports, legs, h, actions, q, small = threshold.detour(length_case)
    chart = apparatus.fixture()
    blanks = {
        frozenset(edge)
        for leg in legs
        for edge in pairwise([chart["hub"], *leg, chart["tip"]])
    }
    assert all(frozenset(edge) in blanks for edge in deposited["physical_cuts"])
    energies, v = eigh(small, driver="evd")
    alpha = source.BASE[:, 0] @ (q[source.ONE] @ v)
    modes = q @ v
    centres = [mode["energy"] for mode in deposited["lowest_and_next_harmonic"]]
    first = energies.searchsorted(16.0)
    second = energies.searchsorted(17 - np.cos(1.5 * np.pi / deposited["length"]))
    np.testing.assert_allclose(energies[[first, second]], centres, rtol=0, atol=9e-13)
    first, second = (
        alpha
        * threshold.clock(
            energies, dict(center=e, tick=threshold.DT, number_delays=threshold.DELAYS)
        )
        for e in centres
    )
    cycle = actions[source.CY]
    assert (h[cycle, :][:, cycle] != h).nnz == 0
    return ports, energies, modes, first, second


def run():
    rng = np.random.default_rng(SEED)
    chart = apparatus.fixture()
    board = apparatus.winding.flight.federation(0)
    labels = calibrated_echo()
    assert conjugation() == 99
    assert (apparatus.native_graph(-1) != apparatus.native_graph(1)).nnz == 0
    owner, boundary, boundary_ports, spoke = idle_boundary(board, chart)
    negative, positive = labels[0, 1], labels[2, 1]
    assert negative != labels[1, 1] and positive != labels[1, 1]
    sensitivity = [
        [list(raw_record(labels, seam, site, q)) for site, q in product(seam, (-1, 1))]
        for seam in boundary
    ]
    assert all(
        labels[1, 1] not in [tuple(x) for x in readings] for readings in sensitivity
    )
    data = json.loads(THRESHOLD.read_text())
    rows = []
    for length_case, deposited in enumerate(data["lengths"]):
        ports, energies, modes, first, second = spectra(length_case, deposited)
        native_weights = position_probabilities(modes @ first)
        cut = deposited["physical_cuts"]
        assert len(set(sum(cut, []))) == 6
        boundary_rows = {
            str(charge): boundary_bank(
                labels, chart, spoke, boundary_ports, ports, native_weights, charge, rng
            )
            for charge in (-1, 1)
        }
        for reading in boundary_rows.values():
            expected = deposited["lowest_and_next_harmonic"][0]
            ref = (
                expected["original_admission"] * expected["mid_and_far_conditional"][1]
            )
            np.testing.assert_allclose(
                reading["signed_tip_arrival_per_original_start"], ref, rtol=1e-6, atol=0
            )
        middle = [
            dict(
                preparation="lowest line",
                time=0,
                kappa=8,
                **arm_bank(labels, cut, ports, native_weights, -1, rng),
            )
        ]
        mix = (first + second) / 2
        for event in deposited["free_beat"][:3]:
            time = event["time"]
            evolution = np.exp(-1j * (energies - energies[0]) * time)
            weights = position_probabilities(modes @ (mix * evolution))
            np.testing.assert_allclose(
                weights.sum(),
                deposited["coherent_source_admission"],
                rtol=2e-12,
                atol=0,
            )
            middle.append(
                dict(
                    preparation="two tones",
                    time=time,
                    kappa=8,
                    **arm_bank(labels, cut, ports, weights, -1, rng),
                )
            )
        if length_case == 0:
            h4 = apparatus.native_graph(-1, 4)
            operation = -1j * time * h4
            changed = expm_multiply(
                operation, modes @ mix, traceA=operation.diagonal().sum()
            )
            result = arm_bank(
                labels, cut, ports, position_probabilities(changed), -1, rng
            )
            clean = middle[-1]["native"]["computed_probabilities"]
            other = result["native"]["computed_probabilities"]
            delta = (
                sum(
                    abs(clean.get(k, 0.0) - other.get(k, 0.0))
                    for k in clean.keys() | other.keys()
                )
                / 2
            )
            assert delta > 1e-10
            middle.append(
                dict(
                    preparation="two tones",
                    time=time,
                    kappa=4,
                    native_vs_changed_original_start_total_variation=float(delta),
                    **result,
                )
            )
        rows.append(
            dict(
                length=deposited["length"],
                cuts=cut,
                boundary=boundary_rows,
                middle=middle,
            )
        )
    return dict(
        seed=SEED,
        simultaneous_bound_tail=TAIL,
        source_threshold_sha256=sha256(THRESHOLD.read_bytes()).hexdigest(),
        starts_per_boundary_bank=BOUNDARY_STARTS,
        starts_per_cut_bank=CUT_STARTS,
        receiving_owner=owner,
        receiving_tip=chart["tip"],
        receiving_spoke=list(spoke),
        physical_boundary=[list(x) for x in boundary],
        boundary_sensitivity_reads=sensitivity,
        source_double_deposit_control=list(labels[0, 0]),
        signed_meter_calibration={f"{a},{b}": list(z) for (a, b), z in labels.items()},
        lengths=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as file:
            file.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    for row in result["lengths"]:
        print(
            row["length"],
            row["boundary"]["-1"]["signed_tip_arrivals"],
            row["middle"][0]["independent_multiple_deposits_given_clock"],
        )
