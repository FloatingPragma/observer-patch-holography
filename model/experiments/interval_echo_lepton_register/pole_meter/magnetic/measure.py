"""Measure the compact 35 line on the magnetic doublet used by the neutral receiver."""

import argparse
import json
from math import atan2, log, pi, sqrt
from pathlib import Path
import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.executable_model import rotations
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet import (
    measure as mag,
)
from model.experiments.interval_echo_gauge.arrival_ladder.spin_transfer import (
    measure as transfer,
)
from model.experiments.interval_echo_lepton_register import measure as rails
from .. import measure as pole
from .certificate import enclose

HERE = Path(__file__).resolve().parent
STEPS = 4096
SHOTS = 2**25


def bank(rng, p):
    rounded = np.round(p, 12)
    if min(rounded) < 0 or abs(1 - sum(rounded)) > 2e-10:
        raise ValueError("Incomplete terminal probabilities")
    return rng.multinomial(SHOTS, rounded / rounded.sum()).tolist()


def compact_analyzer(projector, owner):
    # The accepted two addresses retain the qubit coherently. Only the other
    # ten ports are queried, until the first bright vacancy read.
    # Fix the unread complement instrument across eigensolver bases. Take
    # the first independent projected coordinate columns in port order.
    basis = []
    for column in projector.T:
        remainder = column.copy()
        for v in basis:
            remainder -= np.dot(v, remainder) * v
        if np.linalg.norm(remainder) > 1e-8:
            basis.append(remainder / np.linalg.norm(remainder))
        if len(basis) == 2:
            break
    target = np.column_stack(basis)
    word = np.eye(12, dtype=complex)
    for col in range(2):
        for port in range(11, col, -1):
            col_value, other = (word @ target)[[col, port], col].real
            if abs(other) > 1e-12:
                rails.rotate(word, col, port, -atan2(other, col_value), 1, owner)
    np.testing.assert_allclose(word[:2].conj().T @ word[:2], projector, atol=1e-12)
    return word


def complementary_rail_checks(board, ports):
    """Route each vacancy read to an external double and unwrite its echo."""
    external = [(v, w) for v in ports for w in board.neighbours[v] if w not in ports]
    gateway, reference = min(external)
    local_gate = ports.index(gateway)
    law = mag.LAW
    echo = law.echo()
    bright = np.tile([word == (1, 1) for word in law.SRC], law.DIM // 9)
    sharp = [law.SRC.index((v, 2)) for v in (0, 1)]
    columns = np.eye(law.DIM)[:, sharp]
    for flag, expected in [(bright, [1, 0]), (~bright, [0, 1])]:
        recovered = echo.conj().T @ (flag[:, None] * (echo @ columns))
        np.testing.assert_allclose(recovered, columns * expected, atol=2e-13)
    checks = []
    for port in range(2, 12):
        route = rails.route(port, local_gate)
        swaps = list(zip(route, route[1:]))
        transported = list(range(12))
        for a, b in swaps:
            transported[a], transported[b] = transported[b], transported[a]
        assert transported[local_gate] == port
        for a, b in reversed(swaps):
            transported[a], transported[b] = transported[b], transported[a]
        assert transported == list(range(12))
        checks.append([ports[v] for v in route])
    return dict(
        double_reference_port=reference,
        checked_routes=checks,
        vacancy_and_filled_bright_amplitudes=[1, 0],
    )


def principal_reader(target, owner):
    """Map both complex left singular vectors onto sharp rails, retaining phase."""
    word = np.eye(12, dtype=complex)
    for col in range(2):
        for port in range(11, col, -1):
            a, b = (word @ target)[[col, port], col]
            if abs(b) < 1e-13:
                continue
            rails.phase(word, col, -np.angle(a), 1, owner)
            rails.phase(word, port, -np.angle(b), 1, owner)
            rails.rotate(word, col, port, -atan2(abs(b), abs(a)), 1, owner)
        rails.phase(word, col, -np.angle((word @ target)[col, col]), 1, owner)
    np.testing.assert_allclose(word @ target, np.eye(12)[:, :2], atol=1e-12)
    return word


def equalized_transfer(columns, proj, polar, analyzer, inputs, names, owner, rng):
    from model.experiments.interval_echo_gauge.arrival_ladder.spin_transfer.aperture import (
        compact_instrument,
    )

    donor = compact_instrument(proj, owner)
    common_admission = donor["accepted"]
    weights, v = np.linalg.eigh(columns.conj().T @ proj @ columns)
    target = proj @ columns @ v @ np.diag(weights**-0.5)
    read = principal_reader(target, owner)
    attenuation = np.eye(12, dtype=complex)
    rails.rotate(attenuation, 1, 2, np.arccos(sqrt(weights[0] / weights[1])), 1, owner)
    # The compact input has support on principal rails 0 and 1. Check just
    # rail 2 once, then undo the frame; both successful outputs are unread.
    before_check = attenuation @ read @ common_admission
    after_check = before_check.copy()
    after_check[2, :] = 0
    transfer_matrix = read.conj().T @ after_check
    restricted = transfer_matrix @ columns
    np.testing.assert_allclose(restricted, sqrt(weights[0]) * polar, atol=2e-12)
    # The reference has no intervention. Its successful joint state is
    # (identity tensor polar) times the prepared Bell state.
    bell = restricted.reshape(-1, order="F") / sqrt(2)
    ideal = polar.reshape(-1, order="F") / sqrt(2)
    fidelity = abs(np.vdot(ideal, bell)) ** 2 / np.vdot(bell, bell).real
    banks = []
    donor_banks = []
    for name, initial in zip(names, inputs.T, strict=True):
        rejected = abs((analyzer @ initial)[2:]) ** 2
        later = abs((before_check @ initial)[2]) ** 2
        kept = transfer_matrix @ initial
        donor_refused = donor["refused"] @ initial
        donor_banks.append(
            dict(
                input=name,
                outcomes=bank(
                    rng,
                    np.array(
                        [
                            np.vdot(donor_refused, donor_refused).real,
                            later,
                            np.vdot(kept, kept).real,
                        ]
                    ),
                ),
            )
        )
        banks.append(
            dict(
                input=name,
                outcomes=bank(rng, np.r_[rejected, later, np.vdot(kept, kept).real]),
            )
        )
    return dict(
        common_success=float(weights[0]),
        retained_amplitude_on_second_rail=float(sqrt(weights[0] / weights[1])),
        conditional_reference_fidelity=float(fidelity),
        complete_banks=banks,
        common_donor_banks=donor_banks,
        common_donor_ports=donor["flag_ports"],
        common_donor_residue=donor["max_donor_residue"],
    )


def encoding():
    owner, ports = transfer.stations()
    assert owner == 6
    columns, _, turns = mag.charged_source(owner)
    board, num, den, _, ref, _ = pole.setup()
    sub = num[np.ix_(ports, ports)]
    projector = sub / den
    effect = columns.conj().T @ projector @ columns
    weights, frame = np.linalg.eigh(effect)
    polar = projector @ columns @ frame @ np.diag(weights**-0.5) @ frame.conj().T
    np.testing.assert_allclose(polar @ polar.conj().T, projector, atol=2e-12)
    analyzer = compact_analyzer(projector, owner)
    # Pull actual preparations from the independent neutral-side process.
    starts = transfer.six_sources(columns, owner)
    for name, v in zip(["minimum", "maximum"], frame.T, strict=True):
        starts.append((name, transfer.preparation(columns @ v, owner)))
    inputs = np.column_stack([u[:, rails.path()[0]] for _, u in starts])
    ref_tip = pole.storage_record()["launch_ports"][1]
    source = ports[rails.path()[0]]
    parent = shortest_path(
        board.adjacency(),
        indices=source,
        directed=False,
        unweighted=True,
        return_predecessors=True,
    )[1]
    route = [ref_tip]
    while route[-1] != source:
        route.append(int(parent[route[-1]]))
    return (
        owner,
        ports,
        board,
        sub,
        den,
        columns,
        turns,
        inputs,
        [name for name, _ in starts],
        polar,
        analyzer,
        ref,
        route[::-1],
    )


def measure():
    (
        owner,
        ports,
        board,
        num,
        den,
        columns,
        (r, s),
        inputs,
        names,
        polar,
        analyzer,
        ref,
        route,
    ) = encoding()
    proj = num / den
    h, e0, _, _, _ = pole.spectrum(1)
    np.testing.assert_allclose(h @ ref, e0 * ref, atol=2e-12)
    energies, v = np.linalg.eigh(h)
    full = np.zeros((board.n, len(names)), complex)
    full[ports, :] = inputs
    weights = abs(v.T @ full) ** 2
    times = np.arange(STEPS) * pole.STEP
    amp = np.exp(-1j * np.outer(times, energies - e0)) @ weights
    delays = np.r_[STEPS, 2 * (STEPS - np.arange(1, STEPS))] / STEPS**2
    rng = np.random.default_rng(350620)
    shot = sqrt(2 * log(2 * len(names) / 0.01) / SHOTS)
    tail = 1 / (STEPS * np.sin(pole.STEP * 0.2 / 2)) ** 2
    exact = (inputs.conj() * (proj @ inputs)).sum(axis=0).real
    ramsey = {}
    for name, prob, z in zip(names, exact, amp.T, strict=True):
        all_out = delays @ pole.outcomes(z)
        counts = bank(rng, all_out)
        mean = (counts[0] - counts[1]) / SHOTS
        ramsey[name] = dict(
            counts=counts,
            exact_pole_weight=float(prob),
            expected_quadrature=float(all_out[0] - all_out[1]),
            observed_quadrature=mean,
            simultaneous_interval=[mean - shot - tail, mean + shot],
        )
    coherent = columns.conj().T @ proj @ columns

    def mean(name):
        return ramsey[name]["observed_quadrature"]

    measured = (mean("Z+") + mean("Z-")) / 2 * np.eye(2, dtype=complex)
    for pauli, axis in zip(transfer.PAULI, "XYZ", strict=True):
        measured += (mean(axis + "+") - mean(axis + "-")) / 2 * pauli

    def matrix_form(a):
        return dict(real=a.real.tolist(), imaginary=a.imag.tolist())

    turn_sets = {}
    for label, u in [
        ("idle", np.eye(12)),
        ("first", r),
        ("second", s),
        ("first_after_second", r @ s),
        ("second_after_first", s @ r),
    ]:
        remaining = polar.conj().T @ u.conj().T @ proj @ u @ polar
        rows = []
        for name, initial in zip(names, inputs.T, strict=True):
            first = analyzer @ initial
            kept = analyzer[:2].conj().T @ first[:2]
            last = analyzer @ u @ kept
            joint = np.r_[
                abs(first[2:]) ** 2,
                abs(last[2:]) ** 2,
                np.vdot(last[:2], last[:2]).real,
            ]
            rows.append(
                dict(
                    input=name,
                    complete_counts=bank(rng, joint),
                    both_pass=float(joint[-1]),
                )
            )
        turn_sets[label] = dict(
            conditional_retention_extrema=np.linalg.eigvalsh(remaining).tolist(),
            banks=rows,
        )
    assert np.max(abs(r @ s + s @ r)) < 2e-11
    links = mag.links(1)
    # A magnetic turn is a permutation multiplied by nonzero phases. Every
    # port can be moved to the same fixed gateway, where a compact mode vanishes.
    gate = next(i for i, p in enumerate(ports) if board.degree[p] == 6)
    perms = rotations()
    assert np.all(proj[gate] == 0)
    chosen = []
    for i in range(12):
        g = next(p for p in perms if p[i] == gate)
        phase = mag.compensated_turn(g, links)[i]
        chosen.append(
            dict(
                port=i,
                rotation_index=perms.index(g),
                gateway=gate,
                phase_units=int(phase),
            )
        )
    all_retained = []
    for g in perms:
        phase = mag.compensated_turn(g, links)
        op = np.eye(12)[:, list(g)] @ np.diag(np.exp(2j * pi * phase / 20))
        kept = polar.conj().T @ op.conj().T @ proj @ op @ polar
        all_retained.append(np.linalg.eigvalsh(kept).tolist())
    return dict(
        owner=owner,
        reference_owner=pole.storage_record()["home_cells"][1],
        source_to_reference_route=route,
        compact_complement_checks=complementary_rail_checks(board, ports),
        steps=STEPS,
        step=pole.STEP,
        starts_per_bank=SHOTS,
        frequency_reference=e0,
        shot_radius_99=shot,
        off_pole_leakage=tail,
        exact_effect=matrix_form(coherent),
        measured_effect=matrix_form(measured),
        principal_weights=np.linalg.eigvalsh(coherent).tolist(),
        angle_degrees=np.degrees(
            np.arccos(np.sqrt(np.linalg.eigvalsh(coherent)))
        ).tolist(),
        interval_certificate=enclose(columns, num, den),
        ramsey_banks=ramsey,
        paired_herald_outcomes="first refusal at rails 2..11; next refusal at rails 2..11 after turn; both pass",
        unbiased_transfer=equalized_transfer(
            columns, proj, polar, analyzer, inputs, names, owner, rng
        ),
        turn_banks=turn_sets,
        all_rotation_conditional_retention=all_retained,
        transitive_gateway_witness=chosen,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    doc = json.dumps(measure(), indent=2) + "\n"
    if p.parse_args().write:
        (HERE / "results.json").write_text(doc)
    else:
        print(doc, end="")


if __name__ == "__main__":
    main()
