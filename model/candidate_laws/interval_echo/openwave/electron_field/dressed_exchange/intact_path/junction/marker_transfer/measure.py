"""Replay an exact native-pulse exchange of both material and face orientations."""

import argparse
from collections import Counter, deque
from hashlib import sha256
from itertools import pairwise
import json
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[3]
TARGET = (0, 2, 1, 3)
# Pulse indices are the sorted active material links, then faces 1045 and
# 1028. Compiled by a simultaneous four-input bidirectional graph search.
SEPARATED = (
    0,
    2,
    1,
    4,
    3,
    6,
    4,
    7,
    6,
    9,
    12,
    7,
    6,
    4,
    3,
    9,
    13,
    6,
    4,
    3,
    4,
    6,
    8,
    10,
    5,
    0,
    6,
    4,
    5,
    10,
    8,
)
# A separate Schreier-generator compilation swaps the two markers with
# material origins returned. Appending the material exchange gives a control.
MARKERS_ONLY = (
    2,
    4,
    2,
    9,
    6,
    9,
    2,
    4,
    2,
    11,
    2,
    1,
    3,
    13,
    7,
    4,
    6,
    4,
    7,
    3,
    11,
    2,
    1,
    3,
    13,
    3,
    2,
    4,
    2,
    9,
    6,
    9,
    2,
    4,
    2,
)
# Every odd-use native gate is disabled on the full four-word code after
# the indicated short path. Its conjugate is a native identity on the code.
DISABLE = {1: (0, 4), 2: (0, 4), 3: (0, 4), 12: (7, 9), 13: (2, 3)}


def edge(a, b):
    return tuple(sorted((a, b)))


class Junction:
    def __init__(self, inherited, board):
        self.edges = list(map(tuple, inherited["complete_head_bit_order"]))
        self.bits = {joint: 1 << i for i, joint in enumerate(self.edges)}
        self.start = tuple(inherited["original_load_addresses_by_source"])
        self.faces = [
            inherited["physical_marker_faces_by_source"][str(port)]
            for port in self.start
        ]
        raw_route = inherited["sectors"][0]["routes"]["crossed"][
            "directed_contact_order"
        ]
        tree = {edge(*bond) for bond in raw_route}
        supports = {edge(*joint) for f in self.faces for joint in pairwise((*f, f[0]))}
        self.links = sorted(tree | supports)
        assert len(tree) == len(supports) == 6 and not tree & supports
        assert set(self.links) <= set(board.seams)
        assert all(tuple(sorted(f)) in board.faces() for f in self.faces)
        self.material_only = tuple(self.links.index(edge(*e)) for e in raw_route)
        self.rings = []
        for face in self.faces:
            darts = list(pairwise((*face, face[0])))
            self.rings.append(
                (
                    sum(self.bits[edge(*joint)] for joint in darts),
                    sum(self.bits[a, b] for a, b in darts if a < b),
                )
            )
        self.words = [
            sec["complete_ordered_head_word_inputs"] for sec in inherited["sectors"]
        ]
        self.complement = (1 << len(self.edges)) - 1
        assert self.words[1] == [w ^ self.complement for w in self.words[0]]

    def state_record(self, coordinates, word, sign):
        counts = Counter(
            e[0] if word >> i & 1 else e[1] for i, e in enumerate(self.edges)
        )
        for port in coordinates:
            counts[port] -= sign
        return {p: n for p, n in counts.items() if n}

    def act(self, state, op, sign):
        a, b, word = state
        if op >= len(self.links):
            mask, direction = self.rings[op - len(self.links)]
            return (
                (a, b, word ^ mask)
                if word & mask in (direction, direction ^ mask)
                else state
            )
        u, v = self.links[op]
        if (u in (a, b)) == (v in (a, b)):
            return state
        source, target = (u, v) if u in (a, b) else (v, u)
        bit = self.bits[u, v]
        if (bool(word & bit) == (source == u)) != (sign == 1):
            return state
        return (target, b, word ^ bit) if source == a else (a, target, word ^ bit)

    def reachable(self):
        states = [(*self.start, w) for w in self.words[0]]
        queue, seen = deque(states), set(states)
        while queue:
            state = queue.popleft()
            for op in range(len(self.links) + len(self.rings)):
                target = self.act(state, op, 1)
                if target not in seen:
                    seen.add(target)
                    queue.append(target)
        return len(seen)

    def contiguous_records(self, trace):
        results = [[], []]
        for states in trace:
            for who in range(2):
                direction_bit = 2 if who == 0 else 1
                first, second = [
                    (i, i ^ direction_bit) for i in range(4) if not i & direction_bit
                ]
                delta_first = states[first[0]][2] ^ states[first[1]][2]
                delta_second = states[second[0]][2] ^ states[second[1]][2]
                found = None
                for j, (mask, orientation) in enumerate(self.rings):
                    ports = {state[who] for state in states}
                    if (
                        delta_first == delta_second == mask
                        and all(
                            word & mask in (orientation, orientation ^ mask)
                            for _, _, word in states
                        )
                        and len(ports) == 1
                        and ports <= set(self.faces[j])
                    ):
                        found = j
                results[who].append(found)
        return results

    def replay(self, pulse, words, sign):
        start = tuple((*self.start, w) for w in words)
        trace = [start]
        sectors = [self.state_record(s[:2], s[2], sign) for s in start]
        parity = [1] * len(start)
        for op in pulse:
            current = tuple(self.act(state, op, sign) for state in trace[-1])
            for j, (before, after) in enumerate(zip(trace[-1], current, strict=True)):
                assert self.act(after, op, sign) == before
                assert self.state_record(after[:2], after[2], sign) == sectors[j]
                # In the convention of ascending physical address, this is
                # the externally stipulated odd-worldline comparator.
                parity[j] *= (
                    1 if (before[0] < before[1]) == (after[0] < after[1]) else -1
                )
            trace.append(current)
        return trace, parity


def output_record(lab, pulses, words, sign, target, require_separated):
    history, parity = lab.replay(pulses, words, sign)
    end = history[-1]
    assert end == tuple((*lab.start[::-1], words[j]) for j in target)
    assert parity == [-1] * 4
    collisions = [
        [
            step
            for step, rows in enumerate(history)
            if rows[j][0] // 12 == rows[j][1] // 12
        ]
        for j in range(4)
    ]
    if require_separated:
        assert not any(collisions)
    swaps = sum(op < len(lab.links) for op in pulses)
    amplitude = np.eye(4, dtype=int)[:, target]
    return {
        "pulse_indices": pulses,
        "material_pulses": swaps,
        "ring_pulses": len(pulses) - swaps,
        "full_ordered_material_and_head_records": history,
        "same_cell_occupations_at_pulse_boundaries_by_input": collisions,
        "native_terminal_marker_amplitude_matrix": amplitude.tolist(),
        "supplied_odd_worldline_comparator_terminal_amplitude": (-amplitude).tolist(),
        "local_attached_face_read_by_input_reference_and_pulse_boundary": lab.contiguous_records(
            history
        ),
    }


def matching_reference(lab):
    uncrossed = tuple(
        x for gate, route in DISABLE.items() for x in (*route, gate, *reversed(route))
    )
    for sign, words in zip((1, -1), lab.words, strict=True):
        for gate, route in DISABLE.items():
            displaced, _ = lab.replay(route, words, sign)
            assert all(lab.act(state, gate, sign) == state for state in displaced[-1])
    first, second = Counter(SEPARATED), Counter(uncrossed)
    assert (
        {k for k, v in first.items() if v % 2}
        == set(DISABLE)
        == {k for k, v in second.items() if v % 2}
    )
    joint = first | second
    paths = []
    for pulse, counts in ((SEPARATED, first), (uncrossed, second)):
        padding = tuple(k for k in sorted(joint) for _ in range(joint[k] - counts[k]))
        assert all(a == b for a, b in zip(padding[::2], padding[1::2], strict=True))
        aligned = (*pulse, *padding)
        assert Counter(aligned) == joint
        paths.append(aligned)
    return paths


def ring_control(lab, words, sign):
    results = {}
    for cut in (12, 13):
        history, _ = lab.replay(tuple(k for k in SEPARATED if k != cut), words, sign)
        samples = history[-1]
        results[str(cut)] = {
            "endpoint_occupied_addresses_and_full_word_by_input": samples,
            "expected_swapped_joint_code_columns": [
                i
                for i, actual in enumerate(samples)
                if actual == (*lab.start[::-1], words[TARGET[i]])
            ],
        }
    return results


def source_moments(path):
    document = json.loads(path.read_text())
    cohort = next(
        a
        for a in document["cases"][-1]["stock_arms"]
        if a["preparation"] == "clock_bypassed"
    )
    flight = next(
        a
        for a in cohort["flights"]
        if a["appointment"] == 30.8125 and a["material_bridge_rate"] == 4
    )
    read = next(
        a
        for a in flight["pair_readings"]
        if a["port_face_and_words"] == [1045, 0, 334, 365]
    )
    return np.array(read["receiver_reference_pauli_moments"]), flight[
        "observable_absolute_error_upper_per_pair"
    ]


def local_marginal(density, keep):
    other = tuple(i for i in range(4) if i not in keep)
    order = keep + other
    routed = density.reshape((2,) * 8).transpose(*order, *(4 + i for i in order))
    return np.trace(routed.reshape(4, 4, 4, 4), axis1=1, axis2=3)


def receiver_records(moments):
    pauli = [
        np.eye(2),
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.diag([1, -1]),
    ]
    a = sum(
        moments[i, j] * np.kron(x, y) / 4
        for i, x in enumerate(pauli)
        for j, y in enumerate(pauli)
    )
    phi = np.array([1, 0, 0, 1]) / np.sqrt(2)
    b = np.outer(phi, phi)
    total = np.kron(a, b)
    swap = (0, 1, 2, 3)
    crossed = (2, 1, 0, 3)
    original_bells = np.kron(phi, phi)
    crossed_bells = original_bells.reshape((2,) * 4).transpose(0, 3, 2, 1).ravel()
    p = float(moments[0, 0])
    outcomes = {}
    for name, operation in (
        ("material_only", swap),
        ("both_material_and_markers", crossed),
    ):
        output = (
            total.reshape((2,) * 8)
            .transpose(*operation, *(i + 4 for i in operation))
            .reshape(16, 16)
        )
        bells = [
            float(np.vdot(test, output @ test).real)
            for test in (original_bells, crossed_bells)
        ]
        pair = {}
        for label, legs in (
            ("1045_with_its_original_reference", (0, 1)),
            ("1045_with_reference_from_1028", (0, 3)),
            ("1028_with_its_original_reference", (2, 3)),
            ("1028_with_reference_from_1045", (2, 1)),
        ):
            reduced = local_marginal(output, legs)
            pt = reduced.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
            least = float(np.linalg.eigvalsh(pt)[0])
            pair[label] = {
                "reference_receiver_pt_min_per_original": least,
                "receiver_reference_pauli_moments_per_original": [
                    [float(np.trace(reduced @ np.kron(x, y)).real) for y in pauli]
                    for x in pauli
                ],
            }
        outcomes[name] = {
            "original_then_crossed_bell_pass_probabilities_per_original": bells,
            "original_then_crossed_bell_pass_probabilities_given_arrival": [
                b / p for b in bells
            ],
            "pairings": pair,
        }
    np.testing.assert_allclose(
        outcomes["material_only"][
            "original_then_crossed_bell_pass_probabilities_per_original"
        ],
        outcomes["both_material_and_markers"][
            "original_then_crossed_bell_pass_probabilities_per_original"
        ][::-1],
        atol=2e-17,
    )
    return outcomes


def phase_records(moments, pad):
    p = float(moments[0, 0])
    swap = np.eye(4)[:, TARGET]
    equal = np.eye(4)[:, 0]
    odd_marker = np.array([0, 1, -1, 0]) / np.sqrt(2)
    for state, symmetry in ((equal, 1), (odd_marker, -1)):
        np.testing.assert_array_equal(swap @ state, symmetry * state)
        for external in (1, -1):
            ports = [
                (np.eye(4) + phase * external * swap) @ state / 2 for phase in (1, -1)
            ]
            np.testing.assert_allclose(
                [float(np.vdot(a, a)) for a in ports],
                [int(external * symmetry == s) for s in (1, -1)],
                atol=2e-15,
                rtol=0,
            )
    records = []
    for marker in (0, 1):
        selected = float((p + ((-1) ** marker) * moments[3, 0]) / 2)
        records.append(
            {
                "matching_sharp_marker_selection": [marker, marker],
                "outcome_order": [
                    "ideal_source_refusal",
                    "no_face_arrival",
                    "other_first_marker_word",
                    "final_charge_or_field_refused",
                    "bright",
                    "dark",
                ],
                "native_compared_with_uncrossed_reference": [
                    0,
                    1 - p,
                    p - selected,
                    0,
                    selected,
                    0,
                ],
                "supplied_negative_arm_phase_control": [
                    0,
                    1 - p,
                    p - selected,
                    0,
                    0,
                    selected,
                ],
                "absolute_error_bound_each_event": 2 * pad,
            }
        )
    return {
        "identity_echo_uses_only_native_positive_full_pulses": True,
        "unresolved_original_reference_ensemble": {
            "source_or_fresh_write_refusal": 0,
            "terminal_face_refusal": 1 - p,
            "other_final_refusal": 0,
            "native_bright_dark": [3 * p / 4, p / 4],
            "supplied_negative_phase_bright_dark": [p / 4, 3 * p / 4],
            "absolute_error_bound_each_event": pad,
        },
        "equal_sharp_word_controls": records,
        "known_marker_input_bright_dark_given_admission": {
            "equal_word_for_native_plus_swap": [1, 0],
            "antisymmetric_marker_pair_for_native_plus_swap": [0, 1],
            "equal_word_for_scalar_minus_swap_control": [0, 1],
            "antisymmetric_marker_pair_for_scalar_minus_swap_control": [1, 0],
        },
    }


def run():
    sys.path.insert(0, str(HERE.parents[8]))
    from model.executable_model import federation

    original = HERE.parent / "results.json"
    given = json.loads(original.read_text())
    lab = Junction(given, federation(2))
    face_receipt = FIELD / "fan_path" / "local_launch" / "word_reads.json"
    moments, pad = source_moments(face_receipt)
    # Both isolated terms have the same calibrated relative eigenspace phase.
    two_words = np.array([[1, -1], [-1, 1]])
    active = expm(-0.5j * np.pi * two_words)
    np.testing.assert_allclose(active, [[0, 1], [1, 0]], atol=2e-15, rtol=0)
    matched_cross, matched_return = matching_reference(lab)
    sectors = []
    for sign, words in zip((1, -1), lab.words, strict=True):
        marker_trace, parity = lab.replay(MARKERS_ONLY, words, sign)
        assert marker_trace[-1] == tuple((*lab.start, words[j]) for j in TARGET)
        assert parity == [1] * 4
        ref, sign_ref = lab.replay(matched_return, words, sign)
        assert ref[-1] == ref[0] and sign_ref == [1] * 4
        assert all(a // 12 != b // 12 for rows in ref for a, b, _ in rows)
        sectors.append(
            {
                "charge_sign": sign,
                "collision_free_joint_exchange": output_record(
                    lab, SEPARATED, words, sign, TARGET, True
                ),
                "split_marker_then_material_exchange": output_record(
                    lab, MARKERS_ONLY + lab.material_only, words, sign, TARGET, False
                ),
                "prior_material_only_control": output_record(
                    lab, lab.material_only, words, sign, (0, 1, 2, 3), True
                ),
                "crossed_arm_with_exposure_matched_native_padding": output_record(
                    lab, matched_cross, words, sign, TARGET, True
                ),
                "uncrossed_native_same_exposure_control": {
                    "pulse_indices": matched_return,
                    "full_ordered_material_and_head_records": ref,
                    "terminal_marker_amplitude_matrix": np.eye(4, dtype=int).tolist(),
                    "all_cells_distinct": True,
                },
                "omit_one_of_the_two_addressed_face_pulses_from_the_fixed_schedule": ring_control(
                    lab, words, sign
                ),
            }
        )
    return {
        "junction_receipt_sha256": sha256(original.read_bytes()).hexdigest(),
        "arrival_receipt_sha256": sha256(face_receipt.read_bytes()).hexdigest(),
        "material_addresses_order": lab.links,
        "local_face_ring_order_after_material_pulses": lab.faces,
        "complete_head_bit_order": lab.edges,
        "sharp_origin_lift_reachable_microrecords": lab.reachable(),
        "sectors": sectors,
        "face_arrival_per_original_launch": float(moments[0, 0]),
        "face_arrival_and_pauli_absolute_error_upper": pad,
        "terminal_observer_records": receiver_records(moments),
        "matched_phase_runs": phase_records(moments, pad),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(
            result, json.loads((HERE / "results.json").read_text()), "marker_transfer"
        )
        print(
            "Verified both exact marker swaps, collision-free routing, and the two reference banks."
        )


if __name__ == "__main__":
    main()
