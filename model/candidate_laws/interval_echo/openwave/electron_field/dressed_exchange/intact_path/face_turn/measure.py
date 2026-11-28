"""Check a full material orbit of the physical 1045 rim with the face qubit live."""

import argparse
from collections import Counter, deque
from functools import lru_cache
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
PAIR = (334, 365)


def edge(a, b):
    return tuple(sorted((a, b)))


class LocalTurn:
    def __init__(self, path):
        self.ports = path.cells[-1]
        assert self.ports[0] == 1045 and len(self.ports) == 6
        self.q = path.sign
        self.edges = [*path.field.edges, (1, 5)]
        self.bits = {e: 1 << i for i, e in enumerate(self.edges)}
        self.old_mask = (1 << len(path.field.edges)) - 1
        self.extra = self.bits[1, 5]
        self.faces = [*path.field.rings, (0, 5, 1)]
        assert set(path.field.edges) | {(1, 5)} == set(self.edges)
        assert (1, 5) not in path.field.edges
        assert edge(self.ports[1], self.ports[5]) in path.board.seams
        actual = {frozenset(f) for f in path.board.faces()}
        assert all(frozenset(self.ports[i] for i in f) in actual for f in self.faces)
        self.rings = []
        for f in self.faces:
            darts = list(pairwise((*f, f[0])))
            mask = sum(self.bits[edge(a, b)] for a, b in darts)
            periodic = sum(self.bits[edge(a, b)] for a, b in darts if a < b)
            self.rings.append((mask, periodic))
        assert self.rings[0][0] == PAIR[0] ^ PAIR[1]
        start = PAIR if self.q == 1 else tuple(w ^ self.old_mask for w in PAIR)
        assert set(start) <= set(path.words_at[0])
        self.chart_order_pair = start

    def gauss(self, x, w):
        degrees = [0] * len(self.ports)
        for (a, b), bit in self.bits.items():
            degrees[a if w & bit else b] += 1
        degrees[x] -= self.q
        return degrees

    def ring(self, x, w, f):
        mask, cyclic = self.rings[f]
        if w & mask not in (cyclic, cyclic ^ mask):
            return w
        target = w ^ mask
        assert self.gauss(x, target) == self.gauss(x, w)
        return target

    def allowed_hop(self, a, b, w):
        bond = self.bits[edge(a, b)]
        if ((a < b) == bool(w & bond)) != (self.q == 1):
            return None
        target = w ^ bond
        assert self.gauss(b, target) == self.gauss(a, w)
        return target


def orbit_search(lab, route, source):
    """Enumerate the product graph of ring pulses and prescribed material hops."""
    initial = (0, source)
    queue = deque([initial])
    parents = {initial: None}
    frontiers = [0] * len(route)
    first = None
    same = tuple(source)
    reversed_record = same[::-1]
    while queue:
        key = queue.popleft()
        k, words = key
        frontiers[k] += 1
        if k == len(route) - 1 and words in (same, reversed_record) and first is None:
            first = key
        neighbours = [
            ((k, tuple(lab.ring(route[k], w, f) for w in words)), ("ring", f))
            for f in range(len(lab.rings))
        ]
        if k < len(route) - 1:
            a, b = route[k : k + 2]
            new = tuple(lab.allowed_hop(a, b, w) for w in words)
            if None not in new:
                neighbours.append(((k + 1, new), ("hop", a, b)))
        for nxt, action in neighbours:
            if nxt not in parents:
                parents[nxt] = (key, action)
                queue.append(nxt)
    actions = None
    if first is not None:
        target = first
        actions = []
        while parents[target] is not None:
            target, action = parents[target]
            actions.append(action)
        actions = actions[::-1]
        if first[1] != source:
            assert tuple(lab.ring(0, w, 0) for w in first[1]) == source
            actions.append(("ring", 0))
    return (frontiers, actions)


def full_rim_route(site, sense):
    rim = [1 + (site - 1 + i * sense) % 5 for i in range(6)]
    return (0, *rim, 0)


def exhaustive(lab):
    conditions = []
    for site, sense, added_head in product(range(1, 6), (1, -1), (0, 1)):
        route = full_rim_route(site, sense)
        source = tuple(w | (added_head * lab.extra) for w in lab.chart_order_pair)
        assert lab.gauss(0, source[0]) == lab.gauss(0, source[1])
        full, both = orbit_search(lab, route, source)
        one = [orbit_search(lab, route, (w,)) for w in source]
        conditions.append(
            {
                "added_bit": added_head,
                "sense": sense,
                "physical_material_route": [lab.ports[p] for p in route],
                "common_two_word_states_by_number_of_hops": full,
                "common_two_word_restoring_schedule": both,
                "single_word_states_by_number_of_hops": [x[0] for x in one],
                "single_word_restoring_schedules": [x[1] for x in one],
            }
        )
    return conditions


def undirected_pulse(action):
    return ("hop", *edge(*action[1:])) if action[0] == "hop" else tuple(action)


def apply_pulse(lab, state, action):
    x, word, winding = state
    if action[0] == "ring":
        return x, lab.ring(x, word, action[1]), winding
    a, b = action[1:]
    if x not in (a, b):
        return state
    y = b if x == a else a
    updated = lab.allowed_hop(x, y, word)
    if updated is None:
        return state
    change = 0 if not x or not y else (1 if (y - x) % 5 == 1 else -1)
    return y, updated, winding + change


def matched_reference(lab, orbit, words):
    """Search contact orders with the same exposure and zero rim displacement."""
    initial = tuple((0, w, 0) for w in words)
    counts = Counter(map(undirected_pulse, orbit))
    ops = tuple(sorted(counts))

    @lru_cache(None)
    def search(rest, pair):
        if not any(rest):
            return () if pair == initial else None
        for j, count in enumerate(rest):
            if not count:
                continue
            op = ops[j]
            moved = tuple(apply_pulse(lab, state, op) for state in pair)
            tail = search(rest[:j] + (count - 1,) + rest[j + 1 :], moved)
            if tail is not None:
                return (op, *tail)
        return None

    schedule = search(tuple(counts[op] for op in ops), initial)
    return schedule, search.cache_info().currsize


def replay_route(lab, schedule, words):
    trace = [tuple((0, w, 0) for w in words)]
    for action in schedule:
        trace.append(tuple(apply_pulse(lab, record, action) for record in trace[-1]))
    return trace


def local_attachment(lab, trace, schedule):
    attached = []
    for first, second in trace:
        x, a, _ = first
        y, b, _ = second
        assert x == y and a != b
        eligible = [
            j
            for j, (mask, cyclic) in enumerate(lab.rings)
            if a ^ b == mask
            and a & mask in (cyclic, cyclic ^ mask)
            and x in lab.faces[j]
        ]
        attached.append(eligible[0] if eligible else None)
    material_steps = [j for j, op in enumerate(schedule) if op[0] == "hop"]
    assert all(attached[j] is None or attached[j + 1] is None for j in material_steps)
    return {
        "attached_face_by_prefix_including_zero": attached,
        "material_steps_with_one_attached_face_before_and_after": 0,
        "material_hops": len(material_steps),
    }


def phase_rows(path):
    saved = json.loads(path.read_text())
    arms = saved["cases"][-1]["stock_arms"]
    bypass = next(arm for arm in arms if arm["preparation"] == "clock_bypassed")
    np.testing.assert_allclose(
        [bypass["original_clock_admission"], bypass["fresh_write_per_original"]],
        1,
        atol=5e-15,
        rtol=0,
    )
    flight = next(
        x
        for x in bypass["flights"]
        if x["material_bridge_rate"] == 4 and x["appointment"] == 30.8125
    )
    pair = next(
        row
        for row in flight["pair_readings"]
        if row["port_face_and_words"] == [1045, 0, *PAIR]
    )
    p = pair["reference_joint_success"]
    reference = pair["partial_transpose_minimum_interval"]
    assert 0 < p < 1 and reference[-1] < 0
    return {
        "arrival_pair_probability_per_original_launch": p,
        "arrival_probability_absolute_error_upper": flight[
            "observable_absolute_error_upper_per_pair"
        ],
        "source_reference_partial_transpose_minimum_interval_per_original_launch": reference,
        "ordered_outcomes": [
            "clock_refused",
            "fresh_write_refused",
            "arrival_refused",
            "arrived_but_winding_or_field_return_refused",
            "bright",
            "dark",
        ],
        "native_phase_comparison": [0, 0, 1 - p, 0, p, 0],
        "supplied_minus_reference_phase_control": [0, 0, 1 - p, 0, 0, p],
        "reference_partial_transpose_interval_on_native_bright_output": reference,
    }


def presentation(lab, cases):
    successes = sorted(
        (len(c["common_two_word_restoring_schedule"]), j, c)
        for j, c in enumerate(cases)
        if c["common_two_word_restoring_schedule"] is not None
    )
    attempts = []
    for _, j, case in successes:
        source = tuple(w | case["added_bit"] * lab.extra for w in lab.chart_order_pair)
        orbit = case["common_two_word_restoring_schedule"]
        control, visited = matched_reference(lab, orbit, source)
        attempts.append(
            {
                "zero_based_case": j,
                "reference_search_states": visited,
                "matched_zero_winding_reference_exists": control is not None,
            }
        )
        if control is None:
            continue
        assert Counter(map(undirected_pulse, orbit)) == Counter(control)
        active = replay_route(lab, orbit, source)
        ref = replay_route(lab, control, source)
        assert active[-1] == tuple((0, w, case["sense"] * 5) for w in source)
        assert ref[-1] == tuple((0, w, 0) for w in source)
        steps = [
            x[0][0]
            for i, x in enumerate(active)
            if i == 0 or x[0][0] != active[i - 1][0][0]
        ]
        assert [lab.ports[p] for p in steps] == case["physical_material_route"]
        assert all(a[0] == b[0] for a, b in active)
        return {
            "reference_search_attempts": attempts,
            "selected_route": case["physical_material_route"],
            "added_bit": case["added_bit"],
            "physical_fan_faces": [[lab.ports[i] for i in f] for f in lab.faces],
            "orbit_pulses_in_chart_order": orbit,
            "zero_winding_pulses_in_chart_order": control,
            "each_step_in_physical_port_head_word_and_winding_for_both_words": {
                label: [
                    [(lab.ports[p], w, steps) for p, w, steps in pair] for pair in trace
                ]
                for label, trace in (("orbit", active), ("matched_reference", ref))
            },
            "native_accepted_amplitude_by_arm": [[[1, 0], [0, 1]]] * 2,
            "attachment_audit": local_attachment(lab, active, orbit),
        }
    raise ValueError("No phase-matched route for this signed control")


def mechanical_fixed_sector(lab):
    from model.executable_model.symmetry import rotations

    physical = tuple(port % 12 for port in lab.ports)
    rows = []
    for bit in (0, 1):
        source = lab.chart_order_pair[0] | bit * lab.extra
        fixed = lab.gauss(0, source)
        for steps in range(1, 5):
            mapping = (0, *(1 + (p - 1 + steps) % 5 for p in range(1, 6)))
            assert any(
                tuple(r[p] for p in physical) == tuple(physical[j] for j in mapping)
                for r in rotations()
            )
            rotated = 0
            for (a, b), bitmask in lab.bits.items():
                head = mapping[a if source & bitmask else b]
                joint = edge(mapping[a], mapping[b])
                if head == joint[0]:
                    rotated |= lab.bits[joint]
            changed = lab.gauss(0, rotated)
            assert changed != fixed
            rows.append(
                {
                    "opened_head_bit": bit,
                    "step_angle_in_fifths_of_a_turn": steps,
                    "changed_physical_gauss_sites": [
                        lab.ports[i] for i in range(6) if changed[i] != fixed[i]
                    ],
                }
            )
    return rows


def run():
    sys.path.insert(0, str(HERE.parents[7]))
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.fixture import (
        FanPath,
        archived_fixture,
    )

    two_word_laplacian = np.array([[1, -1], [-1, 1]])
    for rate, duration in ((4, np.pi / 8), (0.5, np.pi)):
        # The zero energy line is fixed and the eigenvalue 2*rate reverses.
        np.testing.assert_allclose(
            expm(-1j * duration * rate * two_word_laplacian),
            [[0, 1], [1, 0]],
            atol=2e-15,
            rtol=0,
        )
    sectors = []
    with archived_fixture() as archived:
        for q in (1, -1):
            lab = LocalTurn(FanPath(archived, 4, q, asymmetric=True))
            conditions = exhaustive(lab)
            sectors.append(
                {
                    "charge_sign": q,
                    "ordered_input_words_ignoring_added_bit": lab.chart_order_pair,
                    "physical_chart_ports": lab.ports,
                    "all_trajectories": conditions,
                    "matched_winding_comparison": presentation(lab, conditions),
                    "rigid_pose_fixed_exterior_control": mechanical_fixed_sector(lab),
                }
            )
    source = FIELD / "fan_path" / "local_launch" / "word_reads.json"
    return {
        "sectors": sectors,
        "arrival_receipt_sha256": sha256(source.read_bytes()).hexdigest(),
        "interferometer_per_original_launch": phase_rows(source),
        "orbit_arrivals_with_original_four_faces_only": 0,
        "full_recovered_field_arrivals_with_added_rim_bond_but_no_fifth_face_pulse": 0,
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
        compare(result, json.loads((HERE / "results.json").read_text()), "face_turn")
        print("Verified every oriented 1045 rim route and its retained face reference.")


if __name__ == "__main__":
    main()
