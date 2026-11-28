"""Check the phase effects and full-cell counts without the measurement code."""

from functools import lru_cache
from itertools import pairwise
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text())


def check_junction(calibration, prepared, case):
    edges = tuple(map(tuple, calibration["complete_head_bit_order"]))
    link_index = {edge: i for i, edge in enumerate(edges)}
    pulse_links = tuple(map(tuple, calibration["material_addresses_order"]))
    faces = calibration["local_face_ring_order_after_material_pulses"]
    pulse_faces = tuple(
        tuple((link_index[tuple(sorted((a, b)))], a) for a, b in pairwise((*f, f[0])))
        for f in faces
    )
    occupancy = frozenset(prepared["original_load_addresses_by_source"])
    q = case["charge_sign"]
    ordered_words = next(
        r["complete_ordered_head_word_inputs"]
        for r in prepared["sectors"]
        if r["charge_sign"] == q
    )
    seeds = tuple(
        (
            occupancy,
            tuple(
                edge[0] if word >> i & 1 else edge[1] for i, edge in enumerate(edges)
            ),
        )
        for word in ordered_words
    )

    def flip_heads(heads, changed):
        directions = list(heads)
        for i in changed:
            directions[i] = sum(edges[i]) - heads[i]
        return tuple(directions)

    def transform(state, gate):
        sites, heads = state
        if gate < len(pulse_links):
            edge = pulse_links[gate]
            i = link_index[edge]
            a, b = edge if q == 1 else edge[::-1]
            source = a if heads[i] == edge[0] else b
            target = sum(edge) - source
            if source in sites and target not in sites:
                return (sites - {source}) | {target}, flip_heads(heads, [i])
        else:
            ring = pulse_faces[gate - len(pulse_links)]
            agreed = sum(heads[i] == u for i, u in ring)
            if agreed in (0, 3):
                return sites, flip_heads(heads, [i for i, _ in ring])
        return state

    gates = range(len(pulse_links) + len(pulse_faces))
    states = list(seeds)
    index = {state: i for i, state in enumerate(states)}
    for state in states:
        for gate in gates:
            target = transform(state, gate)
            if target not in index:
                index[target] = len(states)
                states.append(target)
    size = len(states)
    assert case["orbit_floor"]["orbit_sizes"] == [size]
    assert size == 120
    identity = np.eye(size)
    steps = [
        expm(
            -0.25j
            * np.pi
            * (identity - identity[[index[transform(s, gate)] for s in states]])
        )
        for gate in gates
    ]
    arms = []
    for name in ("intended_crossed", "uncrossed_reference"):
        columns = identity[:, : len(seeds)]
        for gate in case["registered_arms"][name]:
            columns = steps[gate] @ columns
        arms.append(columns)
    bright = (arms[0] + arms[1]) / 2
    dark = (arms[0] - arms[1]) / 2
    assert np.max(np.abs(bright.sum(axis=0) - 1)) < 1e-12
    assert np.max(np.abs(dark.sum(axis=0))) < 1e-12

    banks = (bright[:4], bright[4:], dark[:4], dark[4:])
    for amplitude, expected in zip(
        banks, case["exact_sharp_and_coherent_banks"], strict=True
    ):
        computed = amplitude.conj().T @ amplitude
        numerator = np.asarray(
            expected["gram_real_integer"], dtype=float
        ) + 1j * np.asarray(expected["gram_imaginary_integer"], dtype=float)
        np.testing.assert_allclose(
            computed, numerator / expected["common_denominator"], rtol=0, atol=3e-14
        )


def check_frame(cases):
    from model.executable_model import cell
    from model.executable_model.symmetry import rotations

    shape = cell()
    edges = tuple(shape.seams)
    # Count the whole sharp signed Gauss sector directly, independent of the
    # dynamically reachable word enumeration in the producer.
    residual_degrees = [
        [sum(port in seam for seam in edges[first:]) for port in range(12)]
        for first in range(len(edges) + 1)
    ]

    @lru_cache(None)
    def count_step(first, need):
        if first == len(edges):
            return int(not any(need))
        hits = 0
        a, b = edges[first]
        for head in (a, b):
            if not need[head]:
                continue
            after = need[:head] + (need[head] - 1,) + need[head + 1 :]
            if (
                after[a] <= residual_degrees[first + 1][a]
                and after[b] <= residual_degrees[first + 1][b]
            ):
                hits += count_step(first + 1, after)
        return hits

    normal, complemented = cases
    counts = normal["sharp_gauss_record"]
    assert [5 - n for n in counts] == complemented["sharp_gauss_record"]
    expected = []
    for charge in range(12):
        require = tuple(n + int(port == charge) for port, n in enumerate(counts))
        expected.append(count_step(0, require))
    for comparison in cases:
        hub = comparison["physical_hub"] % 12
        held, mobile = [row["orbit_sizes"][0] for row in comparison["code_orbits"]]
        assert expected[hub] == held and sum(expected) == mobile
        directions = {
            pose["order"]
            for pose in comparison["all_spatial_rotations"]
            if not pose["changed_sharp_gauss_ports"]
        }
        assert directions == {1, 5}
        record = comparison["sharp_gauss_record"]
        count = 0
        for rotation in rotations():
            if tuple(record[v] for v in rotation) != tuple(record):
                continue
            count += 1
            assert all(
                rotation[rotation[v]] != v for v in range(12) if rotation[v] != v
            )
        assert count == 5
    return held, mobile


def main():
    result = load(HERE / "results.json")
    calibration = load(HERE.parent / "results.json")
    prepared = load(HERE.parents[1] / "results.json")
    for case in result["junction_half_pulse_comparison"]:
        check_junction(calibration, prepared, case)
    held, mobile = check_frame(result["one_charge_completed_cell_rotation_check"])
    print(
        f"Independent dense-phase and Gauss counting checks passed on both signs ({held} held, {mobile} mobile)."
    )


if __name__ == "__main__":
    sys.path.insert(0, str(HERE.parents[9]))
    main()
