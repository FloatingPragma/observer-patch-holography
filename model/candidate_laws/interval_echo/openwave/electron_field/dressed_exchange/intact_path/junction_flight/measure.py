"""Resolve two-load return and exchange interference on the actual 1045 tree."""

import argparse
from collections import deque
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
from random import Random

import numpy as np
from scipy.linalg import expm

from ....fan_path.fixture import archived_fixture
from ....fan_recoil.self_start.measure import bound, draw
from ....replay import compare
from model.executable_model import cell


HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[2]
EDGES = [
    (1045, 1010),
    (1045, 1044),
    (1044, 1028),
    (1045, 1052),
    (1052, 1047),
    (1047, 1032),
]
EXCESS_HEAD = [1045, 1044, 1028, 1045, 1052, 1047]
START = (1028, 1045)
RATE = 4
TRIALS = 2**34
GUARD = 2e-10


def enabled(pair, word, sign, cut, distinct_owners):
    for bit, (u, v) in enumerate(EDGES):
        if (cut and bit == 0) or ((u in pair) == (v in pair)):
            continue
        at, dest = (u, v) if u in pair else (v, u)
        other = next(p for p in pair if p != at)
        head = v if word & (1 << bit) else u
        if (head == at) != (sign == 1) or (
            distinct_owners and dest // 12 == other // 12
        ):
            continue
        yield tuple(sorted((other, dest))), word ^ (1 << bit), at, dest, other


def gauss(pair, word, sign, ports):
    heads = [side[(word >> j) & 1] for j, side in enumerate(EDGES)]
    return tuple(heads.count(p) - sign * int(p in pair) for p in ports)


def generators(degrees, sign, cut, distinct_owners):
    initial_word = sum(
        int(v == h) << bit
        for bit, ((_, v), h) in enumerate(zip(EDGES, EXCESS_HEAD, strict=True))
    )
    if sign == -1:
        initial_word ^= 63
    first = (START, initial_word)
    reached = {first}
    todo = deque([first])
    flux = gauss(*first, sign, degrees)
    while todo:
        pair, w = todo.popleft()
        assert gauss(pair, w, sign, degrees) == flux
        for dest, target, *_ in enabled(pair, w, sign, cut, distinct_owners):
            if (dest, target) not in reached:
                reached.add((dest, target))
                todo.append((dest, target))
    states = sorted(reached)
    ids = {pair: i for i, (pair, _) in enumerate(states)}
    assert len(ids) == len(
        states
    )  # checked cut law: the address pair fixes all six heads
    n = len(states)
    h = np.zeros((n, n), dtype=int)
    fermion = h.copy()
    labeled = np.zeros((2 * n, 2 * n), dtype=int)
    opposite_spin = labeled.copy()
    for i, (pair, w) in enumerate(states):
        neighbors = list(enabled(pair, w, sign, cut, distinct_owners))
        diagonal = sum(8 * (6 - degrees[p]) for p in pair) + RATE * len(neighbors)
        h[i, i] = fermion[i, i] = diagonal
        for assignment in (0, 1):
            labeled[2 * i + assignment, 2 * i + assignment] = diagonal
            opposite_spin[2 * i + assignment, 2 * i + assignment] = diagonal
        for dest, target, at, new, other in neighbors:
            j = ids[dest]
            assert states[j][1] == target
            # Ordered spinless Fock creation/annihilation at actual site labels.
            sign_mat = -1 if min(at, new) < other < max(at, new) else 1
            h[i, j] = -RATE
            fermion[i, j] = -RATE * sign_mat
            for assignment in (0, 1):
                initial_up = pair[assignment]
                final_up = new if initial_up == at else initial_up
                next_assignment = dest.index(final_up)
                row, col = 2 * j + next_assignment, 2 * i + assignment
                labeled[row, col] = -RATE
                # Independent one-up/one-down Fock space, ordered spin orbitals.
                spin = int(initial_up != at)
                old_id, new_id, other_id = (
                    2 * at + spin,
                    2 * new + spin,
                    2 * other + 1 - spin,
                )
                pauli_sign = (
                    -1 if min(old_id, new_id) < other_id < max(old_id, new_id) else 1
                )
                opposite_spin[row, col] = -RATE * pauli_sign
    symmetric = np.kron(np.eye(n, dtype=int), [[1], [1]])
    alternating = np.kron(np.eye(n, dtype=int), [[1], [-1]])
    for a in (h, fermion, labeled, opposite_spin):
        np.testing.assert_array_equal(a, a.T)
    np.testing.assert_array_equal(labeled @ symmetric, symmetric @ h)
    np.testing.assert_array_equal(labeled @ alternating, alternating @ fermion)
    np.testing.assert_array_equal(
        opposite_spin @ alternating, alternating @ h
    )  # singlet
    np.testing.assert_array_equal(
        opposite_spin @ symmetric, symmetric @ fermion
    )  # triplet
    return states, ids[START], h, fermion, labeled


def exchange_distance(labeled, first):
    source, goal = 2 * first, 2 * first + 1
    distance = {source: 0}
    todo = deque([source])
    while todo:
        at = todo.popleft()
        for dest in np.flatnonzero(labeled[at] < 0):
            if dest not in distance:
                distance[dest] = distance[at] + 1
                todo.append(dest)
    return distance.get(goal)


def return_series(h, fermion, first, distance):
    # Obtain exact Taylor coefficients of |<first|exp(-iuH)|first>|^2.
    def moments(a):
        a = a.astype(object)
        column = np.eye(len(a), dtype=int)[:, first].astype(object)
        out = []
        for _ in range(distance + 1):
            out.append(int(column[first]))
            column = a @ column
        return out

    def probability(m, order):
        return sum(
            Fraction(
                (-1) ** (k + order // 2) * m[k] * m[order - k],
                factorial(k) * factorial(order - k),
            )
            for k in range(order + 1)
        )

    n, f = moments(h), moments(fermion)
    assert f[:distance] == n[:distance] and distance % 2 == 0
    coefficient = probability(f, distance) - probability(n, distance)
    assert coefficient
    return str(coefficient)


def wave(a, first, u):
    e, v = np.linalg.eigh(a)
    assert np.linalg.norm(v.T @ v - np.eye(len(a)), ord=2) < 1e-13
    assert np.linalg.norm(a @ v - v * e, ord=2) < 1e-12
    result = v @ (np.exp(-1j * u * e) * v[first])
    np.testing.assert_allclose(result, expm(-1j * u * a)[:, first], atol=3e-13, rtol=0)
    return result


def physical_setup():
    with archived_fixture() as (model, _):
        board = model.federation(2)
        ports = sorted({p for edge in EDGES for p in edge})
        for u, v in EDGES:
            assert v in board.neighbours[u]
        degrees = {p: int(board.degree[p]) for p in ports}
    shape = cell()
    frame = json.loads((HERE.parent / "face_turn/cell_frame.json").read_text())
    for chart in frame["proper_pose_checkpoints"]:
        sign = chart["material_sign"]
        words = chart["head_words_at_five_successive_axial_poses"][0]
        for edge, required in zip(EDGES, EXCESS_HEAD, strict=True):
            u, v = edge
            if u // 12 == v // 12:
                assert u // 12 == 1045 // 12
                normalized = tuple(sorted((u % 12, v % 12)))
                bit = tuple(shape.seams).index(normalized)
                for word in words:
                    true_head = (
                        normalized[0] if word >> bit & 1 else normalized[1]
                    ) + 12 * (u // 12)
                    assert (true_head == required) == (sign == 1)
    archived = json.loads((FIELD / "fan_path/local_launch/word_reads.json").read_text())
    arm = next(
        a
        for a in archived["cases"][-1]["stock_arms"]
        if a["preparation"] == "clock_bypassed"
    )
    data = next(
        f
        for f in arm["flights"]
        if f["appointment"] == 30.8125 and f["material_bridge_rate"] == 4
    )
    gate = next(
        r
        for r in data["pair_readings"]
        if r["port_face_and_words"] == [1045, 0, 334, 365]
    )
    return (
        degrees,
        gate["reference_joint_success"],
        data["observable_absolute_error_upper_per_pair"],
    )


def run():
    degrees, p, pad = physical_setup()
    rng = Random(10281045)
    cases = []
    for guarded in (False, True):
        for cut in (False, True):
            states, first, native, fermion, labeled = generators(
                degrees, 1, cut, guarded
            )
            complement, nf, nh, nferm, nl = generators(degrees, -1, cut, guarded)
            assert first == nf and [(pair, 63 ^ w) for pair, w in states] == complement
            for a, b in zip((native, fermion, labeled), (nh, nferm, nl), strict=True):
                np.testing.assert_array_equal(a, b)
            distance = exchange_distance(labeled, first)
            assert (distance is None) == cut
            if cut:
                # In path order the origin order cannot reverse. The signed
                # gauge also verifies equality of ALL address laws at every time.
                gauge = np.zeros(len(native), dtype=int)
                gauge[first] = 1
                queue = deque([first])
                while queue:
                    i = queue.popleft()
                    for j in np.flatnonzero(native[i] < 0):
                        want = gauge[i] * (fermion[i, j] // native[i, j])
                        if not gauge[j]:
                            gauge[j] = want
                            queue.append(j)
                        assert gauge[j] == want
                assert np.all(gauge)
                np.testing.assert_array_equal(fermion, gauge[:, None] * native * gauge)
            coefficient = (
                None if cut else return_series(native, fermion, first, distance)
            )
            reads = []
            separates = False
            for sixteenths in (0, 1, 2, 6, 16):
                u = np.pi * sixteenths / 16
                a, b = (wave(h, first, u) for h in (native, fermion))
                tagged = wave(labeled, 2 * first, u).reshape(-1, 2)
                np.testing.assert_allclose(
                    a, tagged[:, 0] + tagged[:, 1], rtol=0, atol=2e-13
                )
                np.testing.assert_allclose(
                    b, tagged[:, 0] - tagged[:, 1], rtol=0, atol=2e-13
                )
                records = []
                for pred in (
                    abs(a) ** 2,
                    abs(b) ** 2,
                    np.sum(abs(tagged) ** 2, axis=1),
                ):
                    counts = draw([*(p * pred), 1 - p], TRIALS, rng)
                    records.append(
                        {
                            "address_probabilities_given_both_inputs": pred.tolist(),
                            "simulated_counts_in_address_order_then_refusal": counts,
                            "return_per_original": p * pred[first],
                            "return_simultaneous_99pct_interval_per_original": bound(
                                counts[first], TRIALS, pad + GUARD
                            ),
                        }
                    )
                diff = abs(b[first]) ** 2 - abs(a[first]) ** 2
                if cut:
                    np.testing.assert_allclose(
                        abs(a) ** 2, abs(b) ** 2, rtol=0, atol=2e-13
                    )
                ai, bi = (
                    x["return_simultaneous_99pct_interval_per_original"]
                    for x in records[:2]
                )
                separates |= min(ai[1], bi[1]) < max(ai[0], bi[0])
                reads.append(
                    {
                        "wait_in_pi_over_16": sixteenths,
                        "native_polarized_and_labelled": records,
                        "swapped_origin_return_probability_in_tagged_run": abs(
                            tagged[first, 1]
                        )
                        ** 2,
                        "polarized_minus_native_return_given_both": diff,
                    }
                )
            assert separates != cut
            cases.append(
                {
                    "distinct_owner_guard": guarded,
                    "spare_arm_cut": cut,
                    "shortest_odd_exchange_walk": distance,
                    "exact_first_polarized_minus_native_return_coefficient": coefficient,
                    "address_pair_and_excess_head_word": [
                        [list(pair), w] for pair, w in states
                    ],
                    "physical_return_bin_index": first,
                    "instruments": reads,
                }
            )
    return {
        "design_sha256": sha256((HERE / "DESIGN.md").read_bytes()).hexdigest(),
        "ports_degrees": [[v, d] for v, d in degrees.items()],
        "edges_and_starting_excess_head": [
            [*edge, h] for edge, h in zip(EDGES, EXCESS_HEAD, strict=True)
        ],
        "first_arrival_per_original": p,
        "second_material_input_and_spectator_completion_supplied": True,
        "conditional_probability_guard": GUARD,
        "simulated_original_launches_per_bank": TRIALS,
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(data, indent=1) + "\n")
    else:
        compare(
            data, json.loads((HERE / "records.json").read_text()), "junction_flight"
        )
        print(
            "Verified both Y-junction flights and cuts, both charge charts, spin sectors and count banks."
        )


if __name__ == "__main__":
    main()
