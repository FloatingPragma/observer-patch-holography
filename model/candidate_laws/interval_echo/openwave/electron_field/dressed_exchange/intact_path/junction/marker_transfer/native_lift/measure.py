"""Bound and replay native partial-pulse phase comparisons on the 1045 arrival."""

import argparse
from collections import Counter, defaultdict, deque
from fractions import Fraction
from hashlib import sha256
from itertools import pairwise
import json
from pathlib import Path
from runpy import run_path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
TRANSFER = HERE.parent
INTACT_PATH = HERE.parents[2]
FIELD = HERE.parents[4]


def read(path):
    return json.loads(path.read_text())


def explore(seeds, permutation):
    seen = set(seeds)
    waiting = deque(seen)
    while waiting:
        for target in permutation(waiting.popleft()):
            if target not in seen:
                seen.add(target)
                waiting.append(target)
    return seen


def partition(seeds, permutation):
    orbits = []
    for seed in seeds:
        if all(seed not in orbit for orbit in orbits):
            orbits.append(explore([seed], permutation))
    assert set(seeds) <= set.union(*orbits)
    return orbits


def rational(value):
    value = Fraction(value)
    return [value.numerator, value.denominator]


def phase_floor(orbits, code, density, estimated_error):
    block = np.zeros((len(code), len(code)))
    for orbit in orbits:
        mask = np.array([int(word in orbit) for word in code])
        block += np.outer(mask, mask) / len(orbit)
    p = float(np.trace(block @ density).real)
    sharp = [
        rational(Fraction(1, len(next(o for o in orbits if w in o)))) for w in code
    ]
    plus = Fraction(1, len(code)) * sum(
        Fraction(sum(w in o for w in code) ** 2, len(o)) for o in orbits
    )
    assert 0 < plus < Fraction(1, 2)
    return {
        "orbit_sizes": [len(o) for o in orbits],
        "sharp_input_inclusive_bright_lower": sharp,
        "equal_coherent_sum_inclusive_bright_lower": rational(plus),
        "equal_coherent_sum_negative_code_amplitude_upper_against_unit_return": rational(
            1 - 2 * plus
        ),
        "archived_per_original_launch_bright_lower_estimate": p,
        "archived_bright_lower_estimate_absolute_error_upper": estimated_error + 2e-15,
    }


def half_pulses(seed, schedule, involution):
    # A half swap is ((1-i) I + (1+i) P)/2. Retain exact Gaussian numerators.
    current = {seed: (1, 0)}
    for step, gate in enumerate(schedule, start=1):
        following = defaultdict(lambda: [0, 0])
        for word, (real, imaginary) in current.items():
            here = following[word]
            here[0] += real + imaginary
            here[1] += imaginary - real
            there = following[involution(word, gate)]
            there[0] += real - imaginary
            there[1] += imaginary + real
        current = {
            word: (real, imaginary)
            for word, (real, imaginary) in following.items()
            if real or imaginary
        }
        assert (
            sum(x for x, _ in current.values()),
            sum(y for _, y in current.values()),
        ) == (1 << step, 0)
    norm = 1 << (2 * len(schedule))
    assert sum(x * x + y * y for x, y in current.values()) == norm
    assert sum(x for x, _ in current.values()) == (1 << len(schedule))
    assert sum(y for _, y in current.values()) == 0
    return current


def combine(first, second, phase):
    joined = {}
    for word in first.keys() | second.keys():
        a, b = first.get(word, (0, 0))
        c, d = second.get(word, (0, 0))
        joined[word] = (a + phase * c, b + phase * d)
    return joined


def exact_effect(columns, keep, denominator):
    size = len(columns)
    real = [[0] * size for _ in columns]
    imaginary = [[0] * size for _ in columns]
    for i, a in enumerate(columns):
        for j, b in enumerate(columns):
            for word in keep & a.keys() & b.keys():
                x, y = a[word]
                u, v = b[word]
                real[i][j] += x * u + y * v
                imaginary[i][j] += x * v - y * u
    normalized = np.array(
        [
            [
                complex(a / denominator, b / denominator)
                for a, b in zip(row_re, row_im, strict=True)
            ]
            for row_re, row_im in zip(real, imaginary, strict=True)
        ]
    )
    return normalized, {
        "common_denominator": denominator,
        "gram_real_integer": real,
        "gram_imaginary_integer": imaginary,
    }


def junction(parent_ops, parent_input, arrived, pad):
    from model.executable_model import federation

    old = read(TRANSFER / "results.json")
    lab = parent_ops["Junction"](parent_input, federation(2))
    active_count = len(lab.links) + len(lab.rings)
    assert active_count == 14
    cohorts = []
    density = np.kron(arrived, np.eye(2) / 2)
    for signed, archive in zip((1, -1), old["sectors"], strict=True):
        crossed = tuple(
            archive["crossed_arm_with_exposure_matched_native_padding"]["pulse_indices"]
        )
        reference = tuple(
            archive["uncrossed_native_same_exposure_control"]["pulse_indices"]
        )
        assert Counter(crossed) == Counter(reference) and len(crossed) == 37

        def physical(state, op):
            labels, word = state
            u, v, updated = lab.act((*labels, word), op, signed)
            return (tuple(sorted((u, v))), updated)

        starting = tuple(
            (tuple(sorted(lab.start)), word) for word in lab.words[(1 - signed) // 2]
        )

        def neighbors(state):
            for op in range(active_count):
                target = physical(state, op)
                assert physical(target, op) == state
                assert lab.state_record(*target, signed) == lab.state_record(
                    *state, signed
                )
                yield target

        orbits = partition(starting, neighbors)
        assert len(orbits) == 1 and len(orbits[0]) > len(starting)
        floor = phase_floor(orbits, starting, density, 2 * pad)
        support = set.union(*orbits)
        allowed = {
            (positions, word)
            for positions, word in support
            if positions == starting[0][0]
            and all(word & mask in (arrow, arrow ^ mask) for mask, arrow in lab.rings)
        }
        assert allowed == set(starting)
        split = []
        exact_banks = []
        cumulative = np.zeros((4, 4), dtype=complex)
        for phase in (1, -1):
            columns = [
                combine(
                    half_pulses(seed, crossed, physical),
                    half_pulses(seed, reference, physical),
                    phase,
                )
                for seed in starting
            ]
            for accepts, detector in ((True, allowed), (False, support - allowed)):
                effect, raw = exact_effect(
                    columns, detector, 1 << (2 * len(crossed) + 2)
                )
                cumulative += effect
                probability = float(np.trace(effect @ density).real)
                split.append(probability)
                exact_banks.append(
                    {
                        "path_port_phase": phase,
                        "charge_and_field_accepts": accepts,
                        "conditional_bank_on_each_sharp_input": [
                            rational(
                                Fraction(
                                    raw["gram_real_integer"][j][j],
                                    raw["common_denominator"],
                                )
                            )
                            for j in range(4)
                        ],
                        **raw,
                    }
                )
        np.testing.assert_allclose(cumulative, np.eye(4), atol=2e-15, rtol=0)
        for key in ("gram_real_integer", "gram_imaginary_integer"):
            assert exact_banks[0][key] == exact_banks[2][key]
        assert sum(split[:2]) != sum(split[2:])
        mu = float(np.trace(arrived).real)
        assert abs(sum(split) - mu) < 2e-15
        assert (
            split[0] + split[1]
            >= floor["archived_per_original_launch_bright_lower_estimate"]
        )
        cohorts.append(
            {
                "charge_sign": signed,
                "registered_arms": {
                    "intended_crossed": crossed,
                    "uncrossed_reference": reference,
                },
                "same_positive_exposure_on_each_native_address": True,
                "terminal_local_flags_select_exactly_four_whole_codewords": True,
                "orbit_floor": floor,
                "exact_sharp_and_coherent_banks": exact_banks,
                "original_launch_outcomes_order": [
                    "source_refused",
                    "face_arrival_refused",
                    "bright_code_accept",
                    "bright_code_refuse",
                    "dark_code_accept",
                    "dark_code_refuse",
                ],
                "native_original_launch_probabilities": [0.0, 1 - mu, *split],
                "planted_minus_headline_probabilities": [
                    0.0,
                    1 - mu,
                    *split[2:],
                    *split[:2],
                ],
                "original_launch_probability_error_upper_each": 2 * pad + 2e-15,
            }
        )
    assert (
        cohorts[0]["exact_sharp_and_coherent_banks"]
        == cohorts[1]["exact_sharp_and_coherent_banks"]
    )
    return cohorts


def bits_on_faces(seams, faces):
    addresses = {joint: 1 << i for i, joint in enumerate(seams)}
    masks = []
    for a, b, c in faces:
        directed = list(pairwise((a, b, c, a)))
        bits = [addresses[tuple(sorted((u, v)))] for u, v in directed]
        masks.append(
            (
                sum(bits),
                sum(bit for bit, (u, v) in zip(bits, directed, strict=True) if u < v),
            )
        )
    return masks


def completed_cell(frame, arrived, pad):
    from model.executable_model import cell
    from model.executable_model.symmetry import rotations

    shape = cell()
    seams = sorted(shape.seams)
    assert list(shape.seams) == seams
    rings = bits_on_faces(seams, shape.faces())
    material = tuple((u, v, 1 << i) for i, (u, v) in enumerate(seams))
    fixed_port = frame["fixed_charge_port"]
    hub = fixed_port % 12
    physical_owner = fixed_port - hub

    def gauss(state, q):
        load, word = state
        values = [0] * shape.n
        values[load] -= q
        for u, v, bit in material:
            values[u if word & bit else v] += 1
        return tuple(values)

    def all_native(state, q, held):
        load, word = state
        for mask, arrow in rings:
            if word & mask in (arrow, arrow ^ mask):
                yield (load, word ^ mask)
        if not held:
            for u, v, bit in material:
                if load in (u, v) and (bool(word & bit) == (load == u)) == (q == 1):
                    yield (v if load == u else u, word ^ bit)

    results = []
    for pose in frame["proper_pose_checkpoints"]:
        q = pose["material_sign"]
        code = tuple(
            (hub, w) for w in pose["head_words_at_five_successive_axial_poses"][0]
        )
        observed = gauss(code[0], q)
        assert all(gauss(item, q) == observed for item in code)
        table = []
        orders = Counter()
        for rotate in rotations():
            changed = [
                physical_owner + v
                for v in range(shape.n)
                if observed[v] != observed[rotate[v]]
            ]
            power, order = rotate, 1
            while power != tuple(range(shape.n)):
                power = tuple(rotate[v] for v in power)
                order += 1
            if not changed:
                orders[order] += 1
            table.append(
                {
                    "physical_rotation": [physical_owner + v for v in rotate],
                    "order": order,
                    "changed_sharp_gauss_ports": changed,
                }
            )
        assert orders == {1: 1, 5: 4} and len(table) == 60
        floors = []
        for held in (True, False):
            orbits = partition(code, lambda s: all_native(s, q, held))
            assert len(orbits) == 1
            assert all(
                gauss(state, q) == observed
                and all(
                    state in all_native(other, q, held)
                    for other in all_native(state, q, held)
                )
                for state in orbits[0]
            )
            floors.append(
                {
                    "hub_material_held": held,
                    **phase_floor(orbits, code, arrived, 2 * pad),
                }
            )
        results.append(
            {
                "charge_sign": q,
                "physical_hub": fixed_port,
                "sharp_gauss_record": observed,
                "code_orbits": floors,
                "proper_gauss_preserving_rotations_by_order": dict(orders),
                "all_spatial_rotations": table,
            }
        )
    return results


def run():
    prerequisites = {
        "junction": TRANSFER.parent / "results.json",
        "marker_transfer": TRANSFER / "results.json",
        "arrival": FIELD / "fan_path" / "local_launch" / "word_reads.json",
        "whole_cell": INTACT_PATH / "face_turn" / "cell_frame.json",
    }
    inherited = run_path(str(TRANSFER / "measure.py"))
    moments, pad = inherited["source_moments"](prerequisites["arrival"])
    pauli = (
        np.eye(2),
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.diag([1, -1]),
    )
    arrived = sum(value * p for value, p in zip(moments[:, 0], pauli, strict=True)) / 2
    return {
        "input_sha256": {
            name: sha256(path.read_bytes()).hexdigest()
            for name, path in prerequisites.items()
        },
        "junction_half_pulse_comparison": junction(
            inherited, read(prerequisites["junction"]), arrived, pad
        ),
        "one_charge_completed_cell_rotation_check": completed_cell(
            read(prerequisites["whole_cell"]), arrived, pad
        ),
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
        compare(result, read(HERE / "results.json"), "native_lift")
        print(
            "Verified the invariant phase floor, all partial-pulse outcomes and the physical Gauss test."
        )


if __name__ == "__main__":
    sys.path.insert(0, str(HERE.parents[9]))
    main()
