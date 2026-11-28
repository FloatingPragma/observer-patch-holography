"""Replay the jointly readable bridge heads on two through four written fans."""

import argparse
from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import pairwise
import json
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .measure import HERE, TOL


EXPOSURE = 30.8125


def local_terms(path):
    """Describe local contact and ring rules in physical addresses."""
    by_port = defaultdict(list)
    seam_cells = []
    for i, cell in enumerate(path.cells):
        for bit, (local_left, local_right) in enumerate(path.field.edges):
            a, b = cell[local_left], cell[local_right]
            by_port[a].append((b, i, 1 << bit, True))
            by_port[b].append((a, i, 1 << bit, False))
            seam_cells.append((i, bit, a, b))
    bridges = defaultdict(list)
    for i, (a, b) in enumerate(path.bridges):
        bridges[a].append((b, i))
        bridges[b].append((a, i))
    bits = {edge: 1 << i for i, edge in enumerate(path.field.edges)}
    faces = []
    for triangle in path.field.rings:
        oriented = list(pairwise((*triangle, triangle[0])))
        faces.append(
            (
                sum(bits[tuple(sorted(edge))] for edge in oriented),
                sum(bits[u, v] for u, v in oriented if u < v),
            )
        )
    return by_port, bridges, faces, seam_cells


def change(items, index, value):
    return (*items[:index], value, *items[index + 1 :])


def legal_moves(state, sign, by_port, bridges, faces):
    """Generate contact and ring moves while retaining actual bridge records."""
    x, words, heads = state
    for destination, i, mask, first in by_port[x]:
        head_is_at_x = bool(words[i] & mask) == first
        if head_is_at_x == (sign == 1):
            yield (destination, change(words, i, words[i] ^ mask), heads), 8
    for destination, i in bridges[x]:
        if (heads[i] == x) == (sign == 1):
            other = destination if heads[i] == x else x
            yield (destination, words, change(heads, i, other)), 8
    for i, word in enumerate(words):
        for mask, clockwise in faces:
            if (word & mask) in (clockwise, mask ^ clockwise):
                yield (x, change(words, i, word ^ mask), heads), 1


def joint_support(path, source, sign, terms):
    """Check the compact generator against a search that stores all bridge heads."""
    by_port, bridges, faces, seam_cells = terms
    compact = {key: i for i, key in enumerate(path.basis)}
    assert len(compact) == len(path.basis)
    initial = np.flatnonzero(np.any(source != 0, axis=1))
    first_hub = path.cells[0][0]
    assert all(path.basis[i][0] == first_hub for i in initial)
    start_heads = tuple(bridge[0 if sign == 1 else 1] for bridge in path.bridges)
    queue = deque((*path.basis[i], start_heads) for i in initial)
    reached = {int(i): start_heads for i in initial}
    h = path.h.tocsr()
    counts = Counter()
    while queue:
        state = queue.popleft()
        x, words, heads = state
        i = compact[x, words]
        expected = Counter()
        for (dest, next_words, next_heads), doubled_rate in legal_moves(
            state, sign, by_port, bridges, faces
        ):
            j = compact[dest, next_words]
            expected[j] += doubled_rate
            if j in reached:
                assert reached[j] == next_heads
            else:
                reached[j] = next_heads
                queue.append((dest, next_words, next_heads))
            kind = (
                "bridge" if heads != next_heads else "ring" if x == dest else "in_fan"
            )
            counts[kind] += 1
        actual = h.getrow(i)
        assert {
            int(j): -2 * float(amp)
            for j, amp in zip(actual.indices, actual.data, strict=True)
            if j != i and amp != 0
        } == dict(expected)
        assert 2 * h[i, i] == 16 * (6 - int(path.board.degree[x])) + sum(
            expected.values()
        )
    indices = np.array(sorted(reached))
    head_records = []
    load_records = []
    for i in indices:
        x, words = path.basis[i]
        load_records.append(x)
        internal = [a if (words[j] >> bit) & 1 else b for j, bit, a, b in seam_cells]
        head_records.append([*internal, *reached[i]])
    edges = [(a, b) for _, _, a, b in seam_cells] + list(path.bridges)
    return indices, np.array(load_records), np.array(head_records), edges, counts


def enclosing_records(path, sign, load_port, head_end, edges):
    """Derive the all-cut identity from pointwise conserved incoming stock."""
    vertices = np.array([x for cell in path.cells for x in cell])
    initial = np.flatnonzero(load_port == path.cells[0][0])[0]
    signed = sign * (load_port[:, None] == vertices).astype(int)
    incoming = (head_end[:, :, None] == vertices).sum(axis=1)
    np.testing.assert_array_equal(
        incoming - incoming[initial], signed - signed[initial]
    )
    np.testing.assert_array_equal(incoming.sum(axis=1), len(edges))

    region_sets = [set(cell) for cell in path.cells]
    suffixes = [set().union(*region_sets[i:]) for i in range(1, len(path.cells))]
    gateway = {path.bridges[-1][1]}
    right_rim = region_sets[-1] - gateway
    cuts = [*region_sets, *suffixes, gateway, right_rim, set(vertices)]
    for cut in cuts:
        active = [i for i, (u, v) in enumerate(edges) if (u in cut) != (v in cut)]
        heads_inside = np.isin(head_end[:, active], list(cut)).sum(axis=1)
        load_inside = np.isin(load_port, list(cut)).astype(int) * sign
        np.testing.assert_array_equal(
            heads_inside - heads_inside[initial], load_inside - load_inside[initial]
        )
    bridge_ends = head_end[:, -len(path.bridges) :]
    baseline = np.array([bridge[0 if sign == 1 else 1] for bridge in path.bridges])
    flipped = bridge_ends != baseline
    assert np.all(flipped[:, :-1] >= flipped[:, 1:])
    cells_read = flipped.sum(axis=1)
    np.testing.assert_array_equal(
        cells_read,
        np.array(
            [
                next(i for i, cell in enumerate(region_sets) if x in cell)
                for x in load_port
            ]
        ),
    )
    return flipped, cells_read


def run_case(path, sign):
    source = path.source()
    if hasattr(source, "toarray"):
        source = source.toarray()
    np.testing.assert_allclose(source.conj().T @ source, np.eye(2), rtol=0, atol=TOL)
    indices, load_port, head_end, edges, counts = joint_support(
        path, source, sign, local_terms(path)
    )
    pattern, cells_read = enclosing_records(path, sign, load_port, head_end, edges)
    response = expm_multiply(-1j * EXPOSURE * path.h, source)
    weights = abs(response[indices]) ** 2
    np.testing.assert_allclose(weights.sum(axis=0), 1, rtol=0, atol=20 * TOL)
    np.testing.assert_array_equal(path.owner[indices], cells_read)
    patterns = []
    for cell in range(len(path.cells)):
        mask = cells_read == cell
        expected = np.arange(len(path.bridges)) < cell
        np.testing.assert_array_equal(
            pattern[mask], np.broadcast_to(expected, pattern[mask].shape)
        )
        patterns.append(
            {
                "bridge_head_flips": expected.astype(int).tolist(),
                "signed_suffix_head_changes": (sign * expected).astype(int).tolist(),
                "occupied_cell": cell,
                "probability_by_input": weights[mask].sum(axis=0).tolist(),
            }
        )
    return {
        "physical_fans_hub_first": [list(map(int, fan)) for fan in path.cells],
        "physical_bridges_source_to_receiver": [
            list(map(int, edge)) for edge in path.bridges
        ],
        "sign": sign,
        "active_head_stock": len(edges),
        "neutral_background_summed_heads": len(edges) - sign,
        "reached_joint_records": len(indices),
        "compact_chart_records": len(path.basis),
        "directed_steps_checked": dict(counts),
        "once_read_exposure": EXPOSURE,
        "all_start_joint_pattern_probabilities": patterns,
    }


def run():
    from ..fan_path.fixture import FanPath, archived_fixture

    with archived_fixture() as archive:
        cases = [
            run_case(FanPath(archive, size, sign), sign)
            for size in (2, 3, 4)
            for sign in (1, -1)
        ]
    for i in range(0, len(cases), 2):
        for positive, negative in zip(
            cases[i]["all_start_joint_pattern_probabilities"],
            cases[i + 1]["all_start_joint_pattern_probabilities"],
            strict=True,
        ):
            np.testing.assert_allclose(
                positive["probability_by_input"],
                negative["probability_by_input"],
                atol=TOL,
                rtol=0,
            )
    archive_receipt = json.loads((HERE / "results.json").read_text())
    for case, original in zip(cases[:2], archive_receipt["sectors"], strict=True):
        np.testing.assert_allclose(
            case["all_start_joint_pattern_probabilities"][-1]["probability_by_input"],
            original["release"][-1][
                "right_cell_load_and_signed_flux_probability_by_input"
            ],
            atol=TOL,
            rtol=0,
        )
    return {
        "design_sha256": sha256((HERE / "PATH_DESIGN.md").read_bytes()).hexdigest(),
        "two_cell_receipt_sha256": sha256(
            (HERE / "results.json").read_bytes()
        ).hexdigest(),
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(HERE.parent / "replay.py"))["compare"]
        compare(
            result,
            json.loads((HERE / "path_results.json").read_text()),
            "path_charge_records",
        )
        for case in result["cases"]:
            print(
                f"{len(case['physical_fans_hub_first'])} fans {case['sign']:+}, "
                f"{case['reached_joint_records']} written records"
            )


if __name__ == "__main__":
    main()
