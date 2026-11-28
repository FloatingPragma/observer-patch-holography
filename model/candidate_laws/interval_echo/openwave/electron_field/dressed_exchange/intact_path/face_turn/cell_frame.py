"""Compile discrete axial poses of the whole 1045 field with native face swaps."""

import argparse
from functools import lru_cache
from itertools import pairwise
import json
from pathlib import Path
from runpy import run_path
import sys

HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[2]
# Found in a simultaneous two-word breadth-first search from the least
# spectator completion. The certificate below checks the compiled action.
FIRST_POSE = (2, 9, 12, 4, 3, 15, 7, 5, 0, 13, 18)


def edge(a, b):
    return tuple(sorted((a, b)))


def complete(heads, extra, goal):
    missing = tuple(a - b for a, b in zip(goal, heads, strict=True))
    assert min(missing) >= 0

    @lru_cache(None)
    def count(i, rest):
        if i == len(extra):
            return (int(not any(rest)), 0)
        options = []
        for bit, p in enumerate(extra[i]):
            if rest[p] == 0:
                continue
            updated = rest[:p] + (rest[p] - 1,) + rest[p + 1 :]
            n, word = count(i + 1, updated)
            if n:
                options.append((n, word | bit << i))
        return sum(n for n, _ in options), min((w for _, w in options), default=0)

    return count(0, missing)


def pose_certificate(cell, rotation, local, added_bit, raw):
    all_bits = {edge(a, b): 1 << i for i, (a, b) in enumerate(cell.seams)}
    face_charts = [(0, j, j + 1) for j in range(1, 5)]
    first_nine = sorted(
        {edge(a, b) for f in face_charts for a, b in pairwise((*f, f[0]))}
    )
    chart_edges = [(local[a], local[b]) for a, b in (*first_nine, (1, 5))]
    active = {edge(a, b) for a, b in chart_edges}
    extra = [e for e in cell.seams if e not in active]
    assert len(extra) == 20
    start = [w | added_bit << 9 for w in raw]
    known = [0] * cell.n
    for i, (a, b) in enumerate(chart_edges):
        known[a if start[0] >> i & 1 else b] += 1
    opposite = next(i for i in range(cell.n) if i not in local and rotation[i] == i)
    lower = set(range(cell.n)) - set(local) - {opposite}
    assert len(lower) == 5 and {rotation[i] for i in lower} == lower
    # The upper rim needs common indegree 2: one site has no assigned head
    # and only two unused seams. Its other sites already require at least 2.
    upper_options = set.intersection(
        *(set(range(known[i], known[i] + 3)) for i in local[1:])
    )
    assert upper_options == {2} and known[local[0]] == 3
    assert [
        i for i in range(6) if (len(cell.seams) - known[local[0]] - i) % 5 == 0
    ] == [2]
    goal = [2] * cell.n
    for i in lower | {local[0]}:
        goal[i] = 3
    assert sum(goal) == len(cell.seams)
    completions, first = complete(known, extra, goal)
    assert completions
    spectator = sum(all_bits[e] for i, e in enumerate(extra) if not (first >> i & 1))
    seeds = tuple(
        spectator
        | sum(
            all_bits[edge(a, b)]
            for i, (a, b) in enumerate(chart_edges)
            if (a < b) == bool(w & (1 << i))
        )
        for w in start
    )

    def heads(w):
        counts = [0] * cell.n
        for (a, b), bit in all_bits.items():
            counts[a if w & bit else b] += 1
        return counts

    assert heads(seeds[0]) == heads(seeds[1]) == goal
    assert [goal[rotation[i]] for i in range(cell.n)] == goal
    return seeds, extra, first, completions, goal, all_bits, heads


def apply_pair(state, pulse, rings):
    mask, direction = rings[pulse]
    return tuple(
        w ^ mask if w & mask in (direction, direction ^ mask) else w for w in state
    )


def run():
    sys.path.insert(0, str(HERE.parents[7]))
    from model.executable_model import cell
    from model.executable_model.symmetry import rotations
    import numpy as np
    from scipy.linalg import expm

    comparison = json.loads((HERE / "results.json").read_text())["sectors"][0]
    actual = comparison["physical_chart_ports"]
    offset = 12 * (actual[0] // 12)
    local = tuple(p - offset for p in actual)
    added_bit = comparison["matched_winding_comparison"]["added_bit"]
    assert added_bit == 0
    shape = cell()
    turn = next(
        r
        for r in rotations()
        if r[local[0]] == local[0]
        and all(r[local[i]] == local[1 + i % 5] for i in range(1, 6))
    )
    raw = comparison["ordered_input_words_ignoring_added_bit"]
    start, extra, spectator, count, indegrees, bits, heads = pose_certificate(
        shape, turn, local, added_bit, raw
    )
    triangles = shape.faces()
    face_turn = [triangles.index(tuple(sorted(turn[v] for v in f))) for f in triangles]
    rings = []
    for triangle in triangles:
        arrows = list(pairwise((*triangle, triangle[0])))
        mask = sum(bits[edge(a, b)] for a, b in arrows)
        direction = sum(bits[edge(a, b)] for a, b in arrows if a < b)
        rings.append((mask, direction))

    def rotate(word):
        return sum(
            bits[edge(turn[a], turn[b])]
            for (a, b), bit in bits.items()
            if (turn[a] < turn[b]) == bool(word & bit)
        )

    cycle = []
    stage = FIRST_POSE
    for _ in range(5):
        cycle.append(stage)
        stage = tuple(face_turn[j] for j in stage)
    assert stage == FIRST_POSE
    active = tuple(j for part in cycle for j in part)
    # Each reference slot echoes the same face at equal total absolute
    # exposure, using a supplied inverse evolution for its second half.
    generator = np.array([[1, -1], [-1, 1]]) / 2
    forward = expm(-1j * np.pi * generator / 2)
    reverse = expm(1j * np.pi * generator / 2)
    np.testing.assert_allclose(reverse @ forward, np.eye(2), atol=2e-15, rtol=0)
    all_head_bits = (1 << len(bits)) - 1
    checkpoints = []
    for sign in (1, -1):
        source = start if sign == 1 else tuple(w ^ all_head_bits for w in start)
        target = source
        pair = source
        progress = [source]
        for moves in cycle:
            for pulse in moves:
                pair = apply_pair(pair, pulse, rings)
                assert heads(pair[0]) == heads(pair[1]) == heads(source[0])
            target = tuple(map(rotate, target))
            assert pair == target
            progress.append(pair)
        assert pair == source
        checkpoints.append(
            {
                "material_sign": sign,
                "head_words_at_five_successive_axial_poses": progress,
            }
        )

    phase = run_path(str(HERE / "measure.py"))["phase_rows"](
        FIELD / "fan_path" / "local_launch" / "word_reads.json"
    )
    return {
        "fixed_charge_port": actual[0],
        "additional_sharp_spectator_heads": [
            [offset + a, offset + b] for a, b in extra
        ],
        "full_positive_chart_spectator_word": spectator,
        "invariant_head_completions": count,
        "positive_full_cell_head_counts_in_increasing_port_order": indegrees,
        "physical_proper_72_degree_port_permutation_in_increasing_port_order": [
            offset + i for i in turn
        ],
        "addressed_physical_face_swaps_at_each_pose": [
            [[offset + p for p in triangles[j]] for j in part] for part in cycle
        ],
        "zero_turn_reference_each_addressed_face_echo_in_pi_units": ["1/2", "-1/2"],
        "reference_uses_supplied_inverse_ring_evolution": True,
        "proper_pose_checkpoints": checkpoints,
        "native_full_face_swaps": len(active),
        "reference_face_echoes": len(active),
        "native_conditional_phase": 1,
        "phase_bank_per_original_launch": phase,
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(data, indent=2) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(data, json.loads((HERE / "cell_frame.json").read_text()), "cell_frame")
        print("Verified the compiled five axial poses and their matched reference arm.")


if __name__ == "__main__":
    main()
