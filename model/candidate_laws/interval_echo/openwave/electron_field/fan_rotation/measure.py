"""Compare controlled angular reads of the archived 14→24 arrival."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from runpy import run_path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
from scipy.linalg import eigh, eigvalsh
from scipy.sparse.csgraph import connected_components

HERE = Path(__file__).resolve().parent
FIELD = HERE.parent
TIMES = (27.0625, 30.8125)
PHASE = np.pi / 2
INPUTS = {
    "first_ring_odd": np.array([1, 0]),
    "last_ring_odd": np.array([0, 1]),
    "circular_plus": np.array([1, 1j]) / np.sqrt(2),
    "circular_minus": np.array([1, -1j]) / np.sqrt(2),
}


def effect(bra_columns, released):
    amplitudes = bra_columns.conj().T @ released
    return amplitudes.conj().T @ amplitudes


def probabilities(gram):
    probabilities = {
        name: float(np.vdot(v, gram @ v).real) for name, v in INPUTS.items()
    }
    probabilities["equal_mixture"] = float(np.trace(gram).real / 2)
    return {name: float(np.clip(value, 0, 1)) for name, value in probabilities.items()}


def flight(modes, source, times=TIMES):
    energy, vectors = modes
    return [
        vectors @ (np.exp(-1j * t * energy[:, None]) * (vectors.conj().T @ source))
        for t in times
    ]


def physical_pose(ports, local_edges, board=None):
    from model.executable_model import cell, federation
    from model.executable_model.symmetry import rotations

    right = ports["right"]
    gate = right[0] % 12
    seams = set(cell().seams)
    live = {tuple(sorted((right[a] % 12, right[b] % 12))) for a, b in local_edges}
    assert len(live) == 9 and live <= seams
    if board is None:
        board = federation(0)
    assert tuple(sorted((ports["left"][0], right[0]))) in board.external
    nontrivial_at_gate = []
    matching = []
    for r in rotations():
        pose = {tuple(sorted((r[a], r[b]))) for a, b in live}
        if r[gate] == gate and pose == live:
            matching.append(r)
        elif r[gate] == gate:
            nontrivial_at_gate.append(
                {
                    "permutation": list(r),
                    "written_seams_still_at_same_address": len(live & pose),
                }
            )
    assert matching == [tuple(range(12))] and len(nontrivial_at_gate) == 4
    return {
        "proper_turns": len(rotations()),
        "move_gateway_aperture": len(rotations())
        - len(nontrivial_at_gate)
        - len(matching),
        "preserve_gateway_but_change_written_support": nontrivial_at_gate,
        "preserve_aperture_and_written_support": len(matching),
    }


def phase_read(space, receiver, sink, arrivals):
    first_ring = space.projector(4, -1)
    np.testing.assert_array_equal(first_ring @ first_ring, first_ring)
    # exp(-i pi P)=I-2P for the addressed native odd-ring projector.
    ring = np.eye(len(space.basis)) - 2 * first_ring
    outer_bits = [space.bits[6 + a, 6 + b] for a, b in ((1, 2), (4, 5))]
    pulse = np.array(
        [(-1) ** sum(bool(word & bit) for bit in outer_bits) for _, word in space.basis]
    )
    head = np.diag(pulse)
    rz, rx = (
        array.astype(int) for array in (np.diag([-1, 1]), np.array([[0, 1], [1, 0]]))
    )
    identity = np.eye(len(space.basis))
    np.testing.assert_array_equal(ring @ ring, identity)
    np.testing.assert_array_equal(head @ head, identity)
    for columns in (receiver, sink):
        np.testing.assert_array_equal(
            ring @ columns, columns @ np.kron(np.eye(columns.shape[1] // 2), rz)
        )
        head_action = np.kron(np.eye(columns.shape[1] // 2), rx)
        # Overall head-pulse phase is set by the two raw head signs in the chart.
        sign = int(np.sign(np.vdot(head @ columns, columns @ head_action)))
        np.testing.assert_array_equal(head @ columns, sign * columns @ head_action)
        np.testing.assert_array_equal(columns.T @ (ring @ head + head @ ring), 0)

    bank_operators = {
        "repeat_native_ring_vs_hold": (ring @ ring, identity),
        "native_ring_and_supplied_heads_reverse_order": (ring @ head, head @ ring),
        "phase_read_positive_control": (identity, identity),
        "phase_read_negative_control": (-identity, identity),
    }
    records = []
    for t, full, columns in zip(TIMES, arrivals, (receiver, sink), strict=True):
        original = effect(columns, full)
        banks = {}
        for name, (arm, reference) in bank_operators.items():
            bright = effect(columns, (arm + reference) @ full / 2)
            dark = effect(columns, (arm - reference) @ full / 2)
            np.testing.assert_allclose(bright + dark, original, atol=2e-14, rtol=0)
            refused = np.eye(2) - bright - dark
            np.testing.assert_allclose(
                bright + dark + refused, np.eye(2), atol=2e-14, rtol=0
            )
            banks[name] = {
                "bright": probabilities(bright),
                "dark": probabilities(dark),
                "joint_detector_refused": probabilities(refused),
            }
        records.append(
            {"time": t, "includes_source_vacuum": columns is sink, "banks": banks}
        )
    return {
        "receiver_outer_rim_chart_bonds": [[7, 8], [10, 11]],
        "isolated_native_ring_duration_in_pi_units": 1,
        "supplied_head_phase_in_pi_units_each": 1,
        "native_ring_logical_action": rz.tolist(),
        "supplied_head_logical_action": (sign * rx).tolist(),
        "launched_with_ideal_supplied_preparation": True,
        "source_refusal_in_this_ensemble": 0,
        "records": records,
    }


def calibrate_paths(space, triangle, rim, h, phase_bits):
    """Calibrate the relative phase using actual closed matter-plus-ring paths."""
    records = {sense: [] for sense in (1, -1)}
    ring_mask, _ = space.orient(triangle)
    for sense in (1, -1):
        oriented = triangle[::sense]
        for x, word in space.basis:
            if x not in oriented:
                continue
            start = oriented.index(x)
            a, b, c = oriented[start:] + oriented[:start]
            contacts = [tuple(sorted(e)) for e in ((a, b), (b, c), (c, a))]
            assert rim in contacts
            point = (x, word)
            material = 0
            for contact in contacts:
                target = space.material_step(point, contact)
                if target == point:
                    break
                i, j = space.ids[target], space.ids[point]
                assert h[i, j] == -4
                material += phase_bits[i, j]
                point = target
            else:
                target = (a, word ^ ring_mask)
                assert point == target
                i, j = space.ids[x, word], space.ids[target]
                assert h[i, j] == -0.5
                field = int(phase_bits[i, j])
                assert material == sense * space.sign and field == -material
                records[sense].append((int(material), field))
    assert all(records.values())
    return [
        {
            "sense": sense,
            "number_of_closed_paths": len(reads),
            "material_phase_in_quarter_turns": reads[0][0],
            "closing_ring_phase_in_quarter_turns": reads[0][1],
        }
        for sense, reads in records.items()
    ]


def stationary_response(h, modes, matter_bits, ring_bits, heads):
    """Resolve the ground's quadratic response independently of its zero slope."""
    e, vectors = modes
    ground = vectors[:, 0]
    derivative = 1j * h * matter_bits
    second = -h * matter_bits**2
    coupling = vectors[:, 1:].T @ derivative @ ground
    curvature = float(
        ground @ second @ ground - 2 * sum(abs(coupling) ** 2 / (e[1:] - e[0]))
    )
    step = 1 / 1024
    shifts = {}
    for name, exponents in (("material", matter_bits), ("field", ring_bits)):
        levels = [
            eigvalsh(
                h * np.exp(1j * direction * step * exponents), subset_by_index=[0, 0]
            )[0]
            for direction in (1, -1)
        ]
        np.testing.assert_allclose(levels[0], levels[1], atol=2e-13, rtol=0)
        shifts[name] = [float(levels[0] - e[0]), float(levels[1] - e[0])]
    numerical = sum(shifts["material"]) / step**2
    np.testing.assert_allclose(numerical, curvature, atol=1e-7, rtol=5e-5)
    assert curvature > 1e-7
    for angle in (PHASE, step):
        m = h * np.exp(1j * angle * matter_bits)
        f = h * np.exp(-1j * angle * ring_bits)
        gauge = np.exp(-1j * angle * heads)
        np.testing.assert_allclose(
            m, gauge[:, None] * f * gauge[None, :].conj(), atol=6e-15, rtol=0
        )
    return {
        "radian_step": step,
        "ground_linear_slopes_at_zero_exact": {"material": 0, "field": 0, "both": 0},
        "ground_energy_increases_at_plus_minus_finite_probe": shifts,
        "ground_material_curvature_from_variation": curvature,
        "ground_material_curvature_from_finite_probe": numerical,
        "opposite_material_and_field_probes_are_isospectral": True,
    }


def flux_read(space, face, physical, h, modes, source, arrival_columns, arrivals):
    rim = tuple(face[1:])
    bit = space.bits[rim]
    heads = np.array([int(bool(word & bit)) for _, word in space.basis])
    phase_bits = heads[None, :] - heads[:, None]
    paths = calibrate_paths(space, face, rim, h, phase_bits)
    sites = np.array([x for x, _ in space.basis])
    moves = sites[:, None] != sites[None, :]
    matter_bits = phase_bits * moves
    ring_bits = phase_bits * ~moves
    np.testing.assert_array_equal(matter_bits + ring_bits, phase_bits)
    matter_current = 1j * h * matter_bits
    ring_current = 1j * h * ring_bits
    head_operator = np.diag(heads)
    np.testing.assert_array_equal(
        matter_current + ring_current, 1j * (h @ head_operator - head_operator @ h)
    )

    energy, vectors = modes
    source_modes = vectors.T @ source
    delta = energy[:, None] - energy[None, :]
    currents = [
        vectors.T @ current @ vectors for current in (matter_current, ring_current)
    ]
    for current in currents:
        np.testing.assert_allclose(np.diag(current), 0, atol=6e-15, rtol=0)

    fluxes = {}
    for tag, selected in (
        ("matter", matter_bits),
        ("field", ring_bits),
        ("both", phase_bits),
    ):
        phase_records = {}
        for sense in (1, -1):
            gauge = np.exp(-1j * sense * PHASE * heads)
            affected = h * np.exp(1j * sense * PHASE * selected)
            measurements = flight(eigh(affected), source)
            if tag == "both":
                np.testing.assert_allclose(
                    affected,
                    gauge[:, None] * h * gauge[None, :].conj(),
                    atol=6e-15,
                    rtol=0,
                )
                rephased = flight(eigh(affected), gauge[:, None] * source)
                for native, control in zip(arrivals, rephased, strict=True):
                    np.testing.assert_allclose(
                        control, gauge[:, None] * native, atol=2e-12, rtol=0
                    )
            phase_records[str(sense)] = [
                probabilities(effect(read, launched))
                for read, launched in zip(arrival_columns, measurements, strict=True)
            ]
        fluxes[tag] = phase_records

    flow = []
    for t, released in zip(TIMES, arrivals, strict=True):
        integral = t * np.exp(0.5j * t * delta) * np.sinc(t * delta / (2 * np.pi))
        integrated = [
            source_modes.T @ (current * integral) @ source_modes for current in currents
        ]
        change = (
            released.conj().T @ head_operator @ released
            - source.T @ head_operator @ source
        )
        np.testing.assert_allclose(sum(integrated), change, atol=1e-12, rtol=0)
        for matrix in integrated:
            np.testing.assert_allclose(matrix, matrix.conj().T, atol=3e-14, rtol=0)
        flow.append(
            {
                name: {
                    key: float(np.vdot(prep, matrix @ prep).real)
                    for key, prep in INPUTS.items()
                }
                for name, matrix in zip(
                    (
                        "integrated_matter_current",
                        "integrated_field_current",
                        "head_occupation_change",
                    ),
                    (*integrated, change),
                    strict=True,
                )
            }
        )
    return {
        "physical_face_in_orientation_order": [physical[x] for x in face],
        "phase_link_in_orientation_order": [physical[x] for x in rim],
        "material_tour_and_closing_ring_controls": paths,
        "zero_flux_current_expectation": "zero for every real energy eigenvector",
        "native_stationary_ground": stationary_response(
            h, modes, matter_bits, ring_bits, heads
        ),
        "flux_arrival_probabilities": fluxes,
        "native_head_balance": flow,
    }


def captured_currents(space, h, field, source, receiver):
    """Compare currents within the prepared and captured code, keeping spectator words."""
    sites = np.array([x for x, _ in space.basis])
    moves = sites[:, None] != sites[None, :]
    for offset, columns in ((0, source), (6, receiver)):
        for a, b in field.edges:
            bit = space.bits[a + offset, b + offset]
            heads = np.array([int(bool(word & bit)) for _, word in space.basis])
            current = 1j * h * (heads[None, :] - heads[:, None])
            for select in (moves, ~moves):
                np.testing.assert_array_equal(
                    columns.T @ (current * select) @ columns, 0
                )
    return {
        "bonds_tested_at_source_and_receiver_each": len(field.edges),
        "receiving_doublet_includes_all_sender_words": receiver.shape[1] == 8,
        "incoming_projected_material_and_field_current_exact": 0,
        "captured_projected_material_and_field_current_exact": 0,
    }


def run_sector(field, heads, grades, sign, ports):
    from search import JoinedSpace, vacuum_code

    space = JoinedSpace(field, (heads, heads), sign)
    h = space.hamiltonian([6, *grades[0], 6, *grades[1]])
    assert connected_components(h, directed=False, return_labels=False) == 1
    modes = eigh(h)
    assert min(np.diff(modes[0])) > 1e-9
    vacuum, code = vacuum_code(field, heads)
    source = np.zeros((len(h), 2))
    sink = np.zeros_like(source)
    receiver = np.zeros((len(h), 4, 2))
    source[space.portal[0]] = np.einsum("ik,j->ijk", code, vacuum)
    sink[space.portal[1]] = np.einsum("i,jk->ijk", vacuum, code)
    receiver[space.portal[1]] = np.einsum("ia,jk->ijak", np.eye(4), code)
    receiver = receiver.reshape(len(h), 8)
    remote = (receiver, sink)
    arrivals = flight(modes, source)
    physical = ports["left"] + ports["right"]
    return {
        "charge_sign": sign,
        "guide_dimension": len(h),
        "minimum_native_spectral_gap_between_levels": float(min(np.diff(modes[0]))),
        "native_arrival": [
            probabilities(effect(read, full))
            for read, full in zip(remote, arrivals, strict=True)
        ],
        "projected_currents": captured_currents(space, h, field, source, receiver),
        "phase_read": phase_read(space, receiver, sink, arrivals),
        "faces": [
            flux_read(
                space,
                tuple(off + j for j in face),
                physical,
                h,
                modes,
                source,
                remote,
                arrivals,
            )
            for off in (0, 6)
            for face in field.rings
        ],
    }


def run():
    manifest = dict(
        line.split("  ", maxsplit=1)[::-1]
        for line in (FIELD / "records.sha256").read_text().splitlines()
    )
    root = FIELD.parents[4]
    sys.path.insert(0, str(root))
    with TemporaryDirectory() as work, ZipFile(FIELD / "records.zip") as archive:
        directory = Path(work)
        assert sorted(archive.namelist()) == sorted(manifest)
        assert all(
            sha256(archive.read(name)).hexdigest() == digest
            for name, digest in manifest.items()
        )
        (directory / root.name).symlink_to(root, target_is_directory=True)
        archive.extractall(directory)
        package = directory / "experiments" / "interval_echo_collective_spin"
        sys.path.insert(0, str(package))
        from search import fan

        inherited = json.loads((package / "results.json").read_text())["headline"][
            "best_worst_qubit_success"
        ]
        reference, physical = inherited["release"], inherited["actual_ports"]
        assert (
            reference["sign"] == 1
            and reference["head_sectors"][0] == reference["head_sectors"][1]
        )
        assert (
            reference["best_receiver_only_time"],
            reference["best_worst_time"],
        ) == TIMES
        heads = tuple(reference["head_sectors"][0])
        field = fan(4)
        conjugate = tuple(
            sum(i in edge for edge in field.edges) - head
            for i, head in enumerate(heads)
        )
        positive, negative = [
            run_sector(field, words, inherited["degrees"], sign, physical)
            for words, sign in ((heads, 1), (conjugate, -1))
        ]
        np.testing.assert_allclose(
            positive["native_arrival"][0]["equal_mixture"],
            reference["receiver_only_mean_joint_success_at_that_time"],
            atol=1e-12,
            rtol=0,
        )
        for left, right in zip(positive["faces"], negative["faces"], strict=True):
            for field_action, left_action in left["flux_arrival_probabilities"].items():
                right_action = right["flux_arrival_probabilities"][field_action]
                for sense in (1, -1):
                    for i in range(len(TIMES)):
                        np.testing.assert_allclose(
                            list(left_action[str(sense)][i].values()),
                            list(right_action[str(-sense)][i].values()),
                            atol=2e-12,
                            rtol=0,
                        )
        return {
            "source_and_receiver_geometry": physical,
            "release_times": list(TIMES),
            "supplied_phase_magnitude_in_pi_units": PHASE / np.pi,
            "proper_receiver_pose_census": physical_pose(physical, field.edges),
            "charge_conjugation_reverses_calibrated_flux": True,
            "sectors": [positive, negative],
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run()
    if args.output:
        with args.output.open("x") as f:
            f.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(
            record, json.loads((HERE / "results.json").read_text()), "four_face_angular"
        )
        print("Verified the four-face angular controls.")


if __name__ == "__main__":
    main()
