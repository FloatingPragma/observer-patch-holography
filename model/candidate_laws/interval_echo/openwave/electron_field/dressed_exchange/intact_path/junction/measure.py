"""Audit the complete head and reference records of a microscopic Y exchange."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import pairwise, product
import json
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[2]
START = (1045, 1028)
FIRST_LEG = (1045, 1052, 1047, 1032)
SECOND_LEG = (1028, 1044, 1045, 1010)
COMPLETION = (1032, 1047, 1052, 1045, 1044, 1028)
ORIGIN_RETURN = (1010, 1045)


def edge(u, v):
    return tuple(sorted((u, v)))


def routes():
    # Pairs retain the directed route for verification; shutters address
    # undirected contacts. The first two legs are identical in both arms.
    prefix = (*pairwise(FIRST_LEG), *pairwise(SECOND_LEG))
    arms = {
        "crossed": (*prefix, *pairwise(COMPLETION), *pairwise(ORIGIN_RETURN)),
        "uncrossed": (
            *prefix,
            *pairwise(reversed(SECOND_LEG)),
            *pairwise(reversed(FIRST_LEG)),
        ),
    }
    assert Counter(edge(*e) for e in arms["crossed"]) == Counter(
        edge(*e) for e in arms["uncrossed"]
    )
    usage = Counter(edge(*e) for e in arms["crossed"])
    assert len(usage) == 6 and set(usage.values()) == {2}
    assert all((b, a) in arms["crossed"] for a, b in arms["crossed"])
    return arms


def geometry(board, cell):
    offset = 1044
    internal = [(offset + a, offset + b) for a, b in cell.seams]
    face_a = (1045, 1049, 1053)
    face_b = min(f for f in board.faces() if START[1] in f)
    leaf = sorted(edge(*e) for e in pairwise((*face_b, face_b[0])))
    tree = {edge(*e) for e in routes()["crossed"]}
    original = {edge(*e) for e in pairwise((*face_a, face_a[0]))}
    assert not (tree & (original | set(leaf))) and not (original & set(leaf))
    bridges = sorted(tree.intersection(board.external))
    assert len(bridges) == 3 and all(any(p // 12 == 87 for p in e) for e in bridges)
    assert len(tree - set(bridges)) == 3 and tree - set(bridges) <= set(internal)
    assert set(leaf + internal + bridges) <= set(board.seams)
    assert len({p // 12 for p in (1045, 1010, 1028, 1032)}) == 4
    return (
        internal + bridges + leaf,
        {q: face for q, face in zip(START, (face_a, face_b), strict=True)},
        bridges,
    )


def states(edges, bridges, faces, inherited, sign):
    bits = {e: 1 << i for i, e in enumerate(edges)}
    initial = inherited["proper_pose_checkpoints"][0][
        "head_words_at_five_successive_axial_poses"
    ][0]
    other = inherited["proper_pose_checkpoints"][1][
        "head_words_at_five_successive_axial_poses"
    ][0]
    assert other == [w ^ ((1 << 30) - 1) for w in initial]
    # Ordered excess bridge heads. 1045 is the head after its actual arrival.
    heads = {(1010, 1045): 1045, (1028, 1044): 1028, (1032, 1047): 1047}
    assert set(bridges) == set(heads)
    bridge_word = sum(bits[e] for e in bridges if heads[e] == e[0])
    face_b = faces[START[1]]
    arrows = list(pairwise((*face_b, face_b[0])))
    b_mask = sum(bits[edge(*e)] for e in arrows)
    b_directed = sum(bits[edge(a, b)] for a, b in arrows if a < b)
    words = [
        a | bridge_word | b
        for a, b in product(initial, (b_directed, b_directed ^ b_mask))
    ]
    completed_mask = (1 << len(edges)) - 1
    if sign == -1:
        words = [w ^ completed_mask for w in words]
    assert len(set(words)) == 4
    return words, bits


def gauss(positions, word, edges, sign):
    counts = Counter(e[0] if word >> i & 1 else e[1] for i, e in enumerate(edges))
    for port in positions:
        counts[port] -= sign
    return {p: n for p, n in counts.items() if n}


def micro_run(edges, faces, bridges, inherited, sign):
    words, bits = states(edges, bridges, faces, inherited, sign)
    held = {
        p: sum(bits[edge(*e)] for e in pairwise((*f, f[0]))) for p, f in faces.items()
    }
    assert words[0] ^ words[2] == held[START[0]]
    assert words[0] ^ words[1] == held[START[1]]
    assert (
        len({tuple(sorted(gauss(START, w, edges, sign).items())) for w in words}) == 1
    )
    matrix = np.eye(4, dtype=int)
    outcomes = {}
    for name, schedule in routes().items():
        pictures = []
        encoding = np.zeros_like(matrix)
        for j, initial in enumerate(words):
            coordinates = list(START)
            word = initial
            original_record = gauss(coordinates, word, edges, sign)
            path = [(tuple(coordinates), word)]
            factor = 1
            for source, target in schedule:
                assert source in coordinates and target not in coordinates
                allowed = (
                    (source < target) == bool(word & bits[edge(source, target)])
                ) == (sign == 1)
                assert allowed
                mover = coordinates.index(source)
                spectator = coordinates[1 - mover]
                factor *= (
                    -1 if min(source, target) < spectator < max(source, target) else 1
                )
                coordinates[mover] = target
                word ^= bits[edge(source, target)]
                assert gauss(coordinates, word, edges, sign) == original_record
                assert coordinates[0] // 12 != coordinates[1] // 12
                assert all(not ((word ^ initial) & mask) for mask in held.values())
                path.append((tuple(coordinates), word))
            assert tuple(coordinates) == (START[::-1] if name == "crossed" else START)
            assert word == initial
            assert factor == (-1 if name == "crossed" else 1)
            encoding[words.index(word), j] = 1
            pictures.append(path)
        np.testing.assert_array_equal(encoding, matrix)
        assert all(
            [p for p, _ in pic] == [p for p, _ in pictures[0]] for pic in pictures
        )
        outcomes[name] = {
            "directed_contact_order": schedule,
            "native_endpoint_field_code_matrix": encoding.tolist(),
            "calibrated_origin_addresses_by_prefix": [p for p, _ in pictures[0]],
            "complete_head_words_by_prefix_and_input": [
                [w for _, w in pic] for pic in pictures
            ],
            "material_cells_distinct_throughout": True,
            "all_incident_port_gauss_records_retained": True,
            "odd_worldline_sign_control": factor,
        }
    return {
        "charge_sign": sign,
        "complete_ordered_head_word_inputs": words,
        "routes": outcomes,
    }


def reference_register(file):
    data = json.loads(file.read_text())
    cohort = next(
        x
        for x in data["cases"][-1]["stock_arms"]
        if x["preparation"] == "clock_bypassed"
    )
    np.testing.assert_allclose(
        [cohort["original_clock_admission"], cohort["fresh_write_per_original"]],
        1,
        atol=2e-15,
        rtol=0,
    )
    flight = next(
        x
        for x in cohort["flights"]
        if x["appointment"] == 30.8125 and x["material_bridge_rate"] == 4
    )
    gate = next(
        x
        for x in flight["pair_readings"]
        if x["port_face_and_words"] == [1045, 0, 334, 365]
    )
    t = np.array(gate["receiver_reference_pauli_moments"])
    pad = flight["observable_absolute_error_upper_per_pair"]
    p = t[0, 0]
    assert p > 0 and gate["partial_transpose_minimum_interval"][-1] < 0
    bell = float((t[0, 0] + t[1, 1] - t[2, 2] + t[3, 3]) / 4)
    # Contraction of the second, locally supplied Bell pair with the two
    # crossed Bell effects reduces to one quarter of the original effect.
    return {
        "arrival_rate_per_original": float(p),
        "arrival_and_original_pair_probability_error_upper": pad,
        "joint_bell_test_per_original_by_pairing": [bell, bell / 4],
        "joint_bell_test_given_face_gate_by_pairing": [bell / p, bell / (4 * p)],
        "pairing_order": ["original_references", "crossed_references"],
        "pairing_calibration_with_two_supplied_bell_inputs": {
            "identity_field_code": [1, 0.25],
            "imposed_marker_swap": [0.25, 1],
        },
        "original_A_and_A_reference_pt_minimum_interval_per_launch": gate[
            "partial_transpose_minimum_interval"
        ],
        "original_B_and_B_reference_pt_minimum_per_launch": float(-p / 2),
        "crossed_marker_reference_marginals_factor": True,
        "bright_dark_comparison_order": [
            "clock_refused",
            "fresh_write_refused",
            "ideal_second_source_refused",
            "arrival_refused",
            "accepted_but_final_charge_or_field_refused",
            "bright",
            "dark",
        ],
        "native_phase_comparison_per_original": [0, 0, 0, float(1 - p), 0, float(p), 0],
        "supplied_negative_phase_control_per_original": [
            0,
            0,
            0,
            float(1 - p),
            0,
            0,
            float(p),
        ],
    }


def incoming_head_checks(archive, FanPath):
    cases = []
    for sign in (1, -1):
        apparatus = FanPath(archive, 4, sign, asymmetric=True)
        entry, hub = apparatus.bridges[-1]
        assert (entry, hub) == (1010, 1045)
        complement = 0 if sign == 1 else sum(apparatus.field.bits.values())
        words = {334 ^ complement, 365 ^ complement}
        end = apparatus.blocks[3, 0]
        count = 0
        for i in range(end.rows.start, end.rows.stop):
            port, packed = apparatus.basis[i]
            if packed[-1] not in words:
                continue
            assert hub == port
            j = apparatus.ids[entry, packed]
            bridge = apparatus.bridge_term.getrow(i)
            assert bridge.nnz == 2 and bridge[0, j] == -4 and bridge[0, i] == 4
            count += 1
        assert count > 0
        cases.append(
            {
                "charge_sign": sign,
                "all_accepted_spectator_and_face_words_reverse_native_incoming_bridge": count,
            }
        )
    return cases


def run():
    sys.path.insert(0, str(HERE.parents[7]))
    from model.executable_model import cell, federation
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.fixture import (
        FanPath,
        archived_fixture,
    )

    reference = HERE.parent / "face_turn" / "cell_frame.json"
    inherited = json.loads(reference.read_text())
    source = FIELD / "fan_path" / "local_launch" / "word_reads.json"
    board = federation(2)
    edges, faces, bridges = geometry(board, cell())
    np.testing.assert_allclose(
        expm(-0.5j * np.pi * np.array([[1, -1], [-1, 1]])),
        [[0, 1], [1, 0]],
        atol=2e-15,
        rtol=0,
    )
    with archived_fixture() as archive:
        incoming = incoming_head_checks(archive, FanPath)
    return {
        "source_frame_receipt_sha256": sha256(reference.read_bytes()).hexdigest(),
        "arrival_pauli_receipt_sha256": sha256(source.read_bytes()).hexdigest(),
        "original_load_addresses_by_source": START,
        "physical_tree_arms_from_central_hub": [
            SECOND_LEG[2:],
            tuple(reversed(SECOND_LEG[:3])),
            FIRST_LEG,
        ],
        "physical_marker_faces_by_source": faces,
        "complete_head_bit_order": edges,
        "full_native_contact_pulse_exposure_in_pi_units": "1/8",
        "sectors": [
            micro_run(edges, faces, bridges, inherited, sign) for sign in (1, -1)
        ],
        "original_arrival_incoming_bridge_head_checks": incoming,
        "paired_reference_and_phase_records": reference_register(source),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(
            result, json.loads((HERE / "results.json").read_text()), "junction_exchange"
        )
        print("Verified the complete 1045 Y exchange and retained-reference records.")


if __name__ == "__main__":
    main()
