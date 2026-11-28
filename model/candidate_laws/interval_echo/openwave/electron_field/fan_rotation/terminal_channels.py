"""Certify terminal selection rules on the original three- and four-fan guides."""

import argparse
from itertools import product
import json
from pathlib import Path
from runpy import run_path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
FIELD = HERE.parent
RIM_REVERSE = (0, 5, 4, 3, 2, 1)


def terminal_reflection(path, cell):
    """Move physical ports and every directed head of just one terminal fan."""
    ports = path.cells[cell]
    assert cell in (0, path.size - 1) and path.entry[cell] == path.exit[cell] == 0
    assert tuple(path.board.degree[list(ports)]) == (6, 5, 5, 5, 5, 5)
    reflection = RIM_REVERSE
    assert tuple(path.heads[x] for x in reflection) == path.heads
    remap = dict(zip(ports, (ports[x] for x in reflection), strict=True))
    word_reflection = {}
    for word in set().union(*path.words_at):
        reflected = 0
        for (a, b), bit in path.field.bits.items():
            target = tuple(sorted((reflection[a], reflection[b])))
            old_head = a if word & bit else b
            if reflection[old_head] == target[0]:
                reflected |= path.field.bits[target]
        word_reflection[word] = reflected

    ordering = np.array(
        [
            path.ids[
                remap.get(port, port),
                words[:cell] + (word_reflection[words[cell]],) + words[cell + 1 :],
            ]
            for port, words in path.basis
        ]
    )
    np.testing.assert_array_equal(ordering[ordering], np.arange(len(ordering)))
    for term in (path.parked, path.bridge_term):
        assert (term[ordering, :][:, ordering] - term).nnz == 0

    portal_words = path.words_at[0]
    word_order = [portal_words.index(word_reflection[word]) for word in portal_words]
    action = path.code.T @ path.code[word_order]
    sign = int(action[0, 1])
    assert abs(sign) == 1
    np.testing.assert_array_equal(action, [[0, sign], [sign, 0]])
    np.testing.assert_array_equal(path.code[word_order], path.code @ action)
    np.testing.assert_array_equal(path.vacuum[word_order], path.vacuum)
    return ordering, word_order, sign


def integer_decoder(path, columns):
    """Twice the published decoder, with integer coefficients."""
    c = path.size - 1
    block = path.blocks[c, 0]
    amplitudes = columns[block.rows].reshape(*block.shape, columns.shape[1])
    output = np.tensordot((2 * path.code.T).astype(int), amplitudes, axes=(1, c))
    return np.moveaxis(output, 0, -2).reshape(-1, 2, columns.shape[1])


def preparations(path, right_order, right_sign, left):
    source = path.source()
    block = path.blocks[0, 0]
    rows = int(np.prod(block.shape))
    scale = rows**0.5
    integral = np.rint(scale * source).astype(int)
    np.testing.assert_allclose(integral, scale * source, atol=5e-15, rtol=0)
    np.testing.assert_array_equal(integral[right_order], integral)
    preload = integral.copy()
    odd_right = (2 * (path.code[:, 0] - right_sign * path.code[:, 1])).astype(int)
    nest = preload[block.rows].reshape(*block.shape, 2)
    nest *= odd_right[:, None]
    np.testing.assert_array_equal(preload[right_order], -preload)
    columns = np.column_stack((integral, preload))
    if left:
        left_order, _, left_sign = left
        convert = np.array([[1, 1], [left_sign, -left_sign]])
        columns = columns @ np.kron(np.eye(2, dtype=int), convert)
        np.testing.assert_array_equal(columns[left_order], columns * [1, -1, 1, -1])
        np.testing.assert_array_equal(
            columns.T @ columns, 2 * rows * np.diag([1, 1, 2, 2])
        )
        names = ("source_reflection_even", "source_reflection_odd")
    else:
        np.testing.assert_array_equal(columns.T @ columns, rows * np.diag([1, 1, 2, 2]))
        names = ("source_first_outer_ring_odd", "source_last_outer_ring_odd")
    return columns, names


def collect_source_span(witnesses, rows, power):
    """Add at most two independent exact rows of the source effect amplitudes."""
    for i, (first, second) in enumerate(rows):
        if len(witnesses) == 2:
            break
        if first == second == 0:
            continue
        if witnesses:
            a, b = witnesses[0]["source_amplitude_integer_row"]
            if first * b == second * a:
                continue
        witnesses.append(
            {
                "power": power,
                "row": i,
                "source_amplitude_integer_row": [int(first), int(second)],
            }
        )


def moments(path, cols, names, right_sign):
    """Use exact integer Taylor numerators, including arbitrary-size integers."""
    h = (2 * path.h).tocsr()
    values = h.data.astype(np.int64)
    np.testing.assert_array_equal(values, h.data)
    assert np.all(np.diff(h.indptr) > 0)
    weights = values[:, None].astype(object)

    sample = (np.arange(2 * len(path.basis)).reshape(-1, 2) % 11) - 5
    np.testing.assert_array_equal(
        integer_decoder(path, sample), 2 * path.decode(sample)
    )
    queries = {}
    for j, name in enumerate(names):
        for preload, parity in enumerate((1, -1)):
            for vacuum in (False, True):
                if vacuum and name == "source_reflection_odd":
                    continue
                queries[preload * 2 + j, parity, vacuum] = {
                    "preparation": name,
                    "initial_receiver_reflection_parity": parity,
                    "required_receiver_reflection_parity": parity,
                    "requires_uniform_environment_vacuum": vacuum,
                }

    records = []
    spans = {(parity, vacuum): [] for parity, vacuum in product((1, -1), (False, True))}
    cols = cols.astype(object)
    for power in range(33):
        decoded = integer_decoder(path, cols)
        for j, parity in enumerate((1, 1, -1, -1)):
            np.testing.assert_array_equal(
                decoded[:, 0, j] - parity * right_sign * decoded[:, 1, j], 0
            )
        if names[1] == "source_reflection_odd":
            np.testing.assert_array_equal(np.sum(decoded[:, :, [1, 3]], axis=0), 0)
        for (parity, vacuum), witnesses in spans.items():
            start = 0 if parity == 1 else 2
            rows = (
                decoded[:, 0, start : start + 2]
                + parity * right_sign * decoded[:, 1, start : start + 2]
            )
            collect_source_span(
                witnesses, rows.sum(axis=0, keepdims=True) if vacuum else rows, power
            )
        for (j, parity, vacuum), record in list(queries.items()):
            channel = decoded[:, 0, j] + parity * right_sign * decoded[:, 1, j]
            if vacuum:
                channel = [sum(channel)]
            nonzero = next(
                ((index, amp) for index, amp in enumerate(channel) if amp), None
            )
            if nonzero:
                index, amp = nonzero
                records.append(
                    record
                    | {
                        "first_nonzero_power": power,
                        "example_environment_row": None if vacuum else index,
                        "nonzero_integer_numerator": int(amp),
                    }
                )
                del queries[j, parity, vacuum]
        if not queries and all(
            len(witnesses) == (1 if vacuum and path.size == 4 else 2)
            for (_, vacuum), witnesses in spans.items()
        ):
            break
        cols = np.add.reduceat(weights * cols[h.indices], h.indptr[:-1], axis=0)
    span_reads = [
        {
            "initial_receiver_reflection_parity": parity,
            "requires_uniform_environment_vacuum": vacuum,
            "stacked_moment_source_rank_lower_bound": len(witnesses),
            "independent_integer_row_witnesses": witnesses,
        }
        for (parity, vacuum), witnesses in spans.items()
    ]
    return records, list(queries.values()), span_reads


def run_sector(path):
    right = terminal_reflection(path, path.size - 1)
    left = terminal_reflection(path, 0) if path.size == 4 else None
    order, right_words, sign = right
    columns, names = preparations(path, order, sign, left)
    moments_checked, unresolved, span_reads = moments(path, columns, names, sign)

    # Verify the exact right-output action on all arriving environmental rows,
    # independently of any dynamical sample. The other words are fixed by R.
    c = path.size - 1
    block = path.blocks[c, 0]
    rows = np.arange(block.rows.start, block.rows.stop).reshape(block.shape)
    np.testing.assert_array_equal(order[block.rows], rows[..., right_words].ravel())
    if left:
        left_order, left_words, _ = left
        np.testing.assert_array_equal(left_order[block.rows], rows[left_words].ravel())

    return {
        "cell_count": path.size,
        "charge_sign": path.sign,
        "gateways": [c[0] for c in path.cells],
        "bridges": [list(e) for e in path.bridges],
        "enumerated_configurations": len(path.basis),
        "terminal_reflections": [
            {
                "hub": path.cells[cell][0],
                "local_neutral_head_permutation": read[1],
                "ring_code_action": [[0, read[2]], [read[2], 0]],
                "commutes_with_all_local_and_bridge_terms_in_aggregate": True,
            }
            for cell, read in zip((0, c), (left, right), strict=True)
            if read
        ],
        "archived_receiver_vacuum_all_time_output_odd_probability_exact": 0,
        "odd_receiver_preload_all_time_output_even_probability_exact": 0,
        "conditioned_receiver_code_output_dimension_upper_bound": 1,
        "vacuum_environment_odd_source_probability_exact": 0 if left else None,
        "nonzero_witnesses": moments_checked,
        "permitted_branches_without_witness_through_power_32": unresolved,
        "source_span_across_all_times_witnesses": span_reads,
    }


def run():
    sys.path.insert(0, str(FIELD))
    from fan_path.fixture import FanPath, archived_fixture

    with archived_fixture() as archived:
        return {
            "sectors": [
                run_sector(FanPath(archived, n, sign))
                for n, sign in product((3, 4), (1, -1))
            ]
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
            record,
            json.loads((HERE / "terminal_channels.json").read_text()),
            "terminal_channels",
        )
        print("Verified the exact longer-guide terminal channel restrictions.")


if __name__ == "__main__":
    main()
