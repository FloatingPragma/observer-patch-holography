"""Audit which measurements can erase the physical five-bridge record."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import pairwise, product
import json
from pathlib import Path
from runpy import run_path
import sys

import networkx as nx

HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[1]


def physical_loop(board, topology):
    fans = topology["physical_fans"]
    bridges = [tuple(edge) for edge in topology["physical_bridge_order"]]
    linkset = set(board.seams)
    assert all(tuple(sorted(edge)) in board.external for edge in bridges)
    assert len(set(x for edge in bridges for x in edge)) == 2 * len(bridges)
    circuit = []
    for j, fan in enumerate(fans):
        center, *rim = fan
        edges = [(center, x) for x in rim] + list(pairwise(rim))
        assert all(tuple(sorted(e)) in linkset for e in edges)
        incoming = bridges[j - 1][1]
        outgoing = bridges[j][0]
        circuit.extend(pairwise(nx.shortest_path(nx.Graph(edges), incoming, outgoing)))
        circuit.append(bridges[j])
    assert all(
        a[1] == b[0] for a, b in zip(circuit, circuit[1:] + circuit[:1], strict=True)
    )
    assert len({a for a, _ in circuit}) == len(circuit)
    return fans, bridges, circuit


def head_delta(oriented_edges):
    """Reversing an edge originally headed toward v gains one head at u."""
    record = Counter()
    for u, v in oriented_edges:
        record[u] += 1
        record[v] -= 1
    return {site: delta for site, delta in sorted(record.items()) if delta}


def permitted_masks(edges, owner):
    """Count the subsets of marked bonds whose reversal preserves a record."""
    passed = []
    for mask in range(1 << len(edges)):
        flipped = [edge for i, edge in enumerate(edges) if mask & (1 << i)]
        if not head_delta((owner(a), owner(b)) for a, b in flipped):
            passed.append(mask)
    return passed


def bare_read(bridges_count):
    """All parts of the hypothetical field-only X read, including sector refusal."""
    outcomes = []
    n = 1 << bridges_count
    # There is at least one unmeasured and mutually orthogonal internal
    # head word. It eliminates port interference even before the Gauss check.
    retained = Fraction(1, 2 * n**2)
    original = Fraction(1, 2 * n)
    for bits in product((1, -1), repeat=bridges_count):
        for phase_port in ("bright", "dark"):
            outcomes.append(
                {
                    "bare_X_outcomes": list(bits),
                    "phase_port": phase_port,
                    "original_gauss_record_probability": str(retained),
                    "gauss_record_refused_probability": str(original - retained),
                }
            )
    assert (
        sum(
            Fraction(r["original_gauss_record_probability"])
            + Fraction(r["gauss_record_refused_probability"])
            for r in outcomes
        )
        == 1
    )
    return {
        "bare_bridge_bits_measured": bridges_count,
        "identical_for_plus_and_minus_exchange_phase": True,
        "original_gauss_record_total_probability": str(Fraction(1, n)),
        "gauss_record_refused_total_probability": str(1 - Fraction(1, n)),
        "outcomes": outcomes,
    }


def whole_loop_read():
    """Read (I±S)/2 on the two directed words after the balanced phase check."""
    comparisons = []
    for sign in (1, -1):
        outcomes = []
        for loop_read, phase_port in product((1, -1), repeat=2):
            # Unknown exchange sign is an input to this sensitivity control;
            # the whole-loop projector is identical in both runs.
            weight = Fraction((1 + sign * loop_read * phase_port) ** 2, 8)
            outcomes.append(
                {
                    "whole_loop_read": loop_read,
                    "interferometer_port": "bright" if phase_port == 1 else "dark",
                    "probability": str(weight),
                }
            )
        assert sum(Fraction(r["probability"]) for r in outcomes) == 1
        comparisons.append(
            {
                "supplied_test_arm_sign": sign,
                "joint_outcomes": outcomes,
                "original_gauss_record_refused_probability": 0,
                "other_loop_pattern_refused_probability": 0,
            }
        )
    return comparisons


def run():
    sys.path.insert(0, str(HERE.parents[6]))
    from model.executable_model import federation

    topology = json.loads((HERE / "results.json").read_text())[
        "physical_cycle_bridge_which_way_control"
    ]
    fans, bridges, loop = physical_loop(federation(2), topology)
    all_bridges = (1 << len(bridges)) - 1
    # Cell-wide charge conservation alone permits only the two windings.
    assert permitted_masks(bridges, lambda p: p // 12) == [0, all_bridges]
    # At physical-port resolution the bare bridge flips need internal corrections.
    assert permitted_masks(bridges, lambda p: p) == [0]
    full_mask = (1 << len(loop)) - 1
    assert permitted_masks(loop, lambda p: p) == [0, full_mask]
    assert not head_delta(loop)
    affected_fan_ports = [
        {p: d for p, d in head_delta(bridges).items() if p in fan} for fan in fans
    ]
    assert all(sorted(changes.values()) == [-1, 1] for changes in affected_fan_ports)
    unchanged_read = [
        {"bridge_word": word, "interferometer_port": port, "probability": "1/4"}
        for word, port in product(("initial", "complement"), ("bright", "dark"))
    ]
    return {
        "physical_fan_hubs": [fan[0] for fan in fans],
        "physical_directed_loop": [list(edge) for edge in loop],
        "physical_directed_bridges": [list(edge) for edge in bridges],
        "per_fan_changed_port_head_counts_if_only_bare_bridges_flipped": [
            [[p, d] for p, d in changes.items()] for changes in affected_fan_ports
        ],
        "bare_bridge_masks_preserving_summed_cell_gauss_record": [0, all_bridges],
        "bare_bridge_masks_preserving_every_port_gauss_record": [0],
        "full_loop_masks_preserving_every_port_gauss_record": [0, full_mask],
        "same_sharp_original_local_head_count_sector_available_on_both_arms": False,
        "local_bridge_diagonal_read_both_signs": unchanged_read,
        "raw_single_edge_and_all_bridge_eraser_records": [
            bare_read(n) for n in (1, len(bridges))
        ],
        "supplied_whole_loop_comparison": whole_loop_read(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(
            result,
            json.loads((HERE / "bridge_eraser.json").read_text()),
            "bridge_eraser",
        )
        print("Verified the local erasure obstruction and whole-loop comparison.")


if __name__ == "__main__":
    main()
