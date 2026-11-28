"""Connect the restored 48-cell route to the archived 1047 exterior meter."""

import argparse
from collections import Counter, deque
from itertools import pairwise
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from ...apron import instrument
from ...wheel import (
    A,
    CENTER,
    EDGES,
    FACES,
    FEED,
    LAND,
    ROUTE,
    edge,
    face_gate,
    flip,
    head,
    in_sector,
    preparation,
    project_word,
    raw_word,
    wheel_basis,
)


HERE = Path(__file__).resolve().parent
PROPAGATION = (
    HERE.parents[5] / "fan_recoil" / "fresh_kick" / "propagation" / "results.json"
)
SOURCE = tuple(map(raw_word, (334, 365)))
LOW_MASK = (1 << 9) - 1


def act(state, pulse, signed=1, edges=EDGES):
    x, word = state
    if len(pulse) == 3:
        return x, flip(word, face_gate(edges, pulse))
    p, q = pulse
    seam = edge(p, q)
    bit = edges.index(seam)
    if x in seam and (head(seam, word, bit) == x) == (signed == 1):
        return p + q - x, word ^ (1 << bit)
    return state


def search(seeds, menu):
    paths = {seeds: ()}
    queued = deque([seeds])
    while queued:
        states = queued.popleft()
        for pulse in menu:
            target = tuple(act(z, pulse) for z in states)
            if target not in paths:
                paths[target] = (*paths[states], pulse)
                queued.append(target)
    return paths


def has_shared_route_orientation(states):
    edges = (FEED, edge(LAND, CENTER))
    return all(
        head(seam, w, EDGES.index(seam)) == source
        for _, w in states
        for seam, source in zip(edges, (A, LAND), strict=True)
    )


def follow(start, pulses, signed=1, edges=EDGES):
    history = [start]
    for pulse in pulses:
        history.append(tuple(act(z, pulse, signed, edges) for z in history[-1]))
    return history


def primary_layout():
    """The exact collector-compatible load, using the original two wheel words."""
    archive = preparation()
    before = tuple(map(tuple, archive["ordered_two_word_trace"][0]))
    paths = search(before, FACES)
    eligible = [
        (states, route)
        for states, route in paths.items()
        if tuple(w & LOW_MASK for _, w in states) == SOURCE
        and len({w >> 9 for _, w in states}) == 1
        and has_shared_route_orientation(states)
    ]
    assert len(paths) == 6 and len(eligible) == 1
    start, backwards = eligible[0]
    prealign = tuple(reversed(backwards))
    load = tuple(map(tuple, archive["native_address_schedule"]))
    history = follow(start, (*prealign, *load))
    assert history[len(prealign) :] == [
        tuple(map(tuple, row)) for row in archive["ordered_two_word_trace"]
    ]
    assert all(x == A for pair in history[: len(prealign) + 1] for x, _ in pair)
    assert all(
        tuple(w & LOW_MASK for _, w in pair) == SOURCE
        for pair in history[: len(prealign) + 1]
    )
    return {
        "commissioned_background": start[0][1] >> 9,
        "paired_wheel_face_orbit_size_at_1045": len(paths),
        "compatible_incoming_pair_count_in_that_orbit": len(eligible),
        "native_address_schedule": [*prealign, *load],
        "prealignment_rings": prealign,
        "archived_load_schedule": load,
        "ordered_two_word_trace": history,
        "wheel_input_words": [project_word(w) for _, w in history[-1]],
    }


def alternate_layout():
    # Same original four-address loading problem, without opening the additional
    # four faces used to align the archive's exact environment.
    face_mask, _ = face_gate(EDGES, ROUTE[2])
    matches = {}
    for background in range(1 << (len(EDGES) - 9)):
        start = tuple((A, w | background << 9) for w in SOURCE)
        possibilities = search(start, ROUTE)
        for states, route in possibilities.items():
            (pos, u), (other_pos, v) = states
            if (
                pos == other_pos == LAND
                and u ^ v == face_mask
                and all(in_sector(pos, project_word(w)) for w in (u, v))
            ):
                matches[background] = (start, route)
                break
    assert list(matches) == [186, 218]
    compatible = [
        (start, route)
        for start, route in matches.values()
        if has_shared_route_orientation(start)
    ]
    assert len(compatible) == 1
    start, route = compatible[0]
    history = follow(start, route)
    return {
        "backgrounds_admitted_by_original_loading_check": list(matches),
        "commissioned_background": start[0][1] >> 9,
        "native_address_schedule": route,
        "ordered_two_word_trace": history,
        "wheel_input_words": [project_word(w) for _, w in history[-1]],
    }


def gauss(states, edges, signed):
    result = []
    for x, w in states:
        counts = Counter(head(seam, w, i) for i, seam in enumerate(edges))
        counts[x] -= signed
        result.append({v: n for v, n in counts.items() if n})
    return result


def collector_replay(layout, ports, source_words):
    route = tuple(edge(a, b) for a, b in pairwise(ports))
    tail = [seam for seam in route if seam not in EDGES]
    assert set(route) & set(EDGES) == {FEED, edge(LAND, CENTER)}
    joints = [*EDGES, *tail]
    assert len(set(joints)) == len(joints)
    # Every vacancy route arrow starts pointing at its next address.
    commissioned = sum(
        int(seam[0] == dst) << joints.index(seam)
        for seam, dst in zip(route, ports[1:], strict=True)
        if seam in tail
    )
    complement = (1 << len(EDGES)) - 1
    local = layout["ordered_two_word_trace"]
    start = tuple((x, (w ^ complement) | commissioned) for x, w in local[0])
    assert [w & LOW_MASK for _, w in start] == list(map(raw_word, source_words))
    for x, w in start:
        assert x == ports[0] and all(
            head(seam, w, joints.index(seam)) == dst
            for seam, dst in zip(route, ports[1:], strict=True)
        )
    out = follow(start, route, -1, joints)
    returned = follow(out[-1], tuple(reversed(route)), -1, joints)
    assert returned == list(reversed(out))
    # This verifies all position/prefix encodings, as consumed by the
    # collector acting within the guide. It does not model the collector clock.
    sectors = gauss(start, joints, -1)
    assert all(
        all(x == port for x, _ in states) and gauss(states, joints, -1) == sectors
        for port, states in zip(ports, out, strict=True)
    )
    joined = follow(returned[-1], layout["native_address_schedule"], -1, joints)
    expected = [
        tuple((x, (w ^ complement) | commissioned) for x, w in states)
        for states in local
    ]
    assert joined == expected and all(gauss(s, joints, -1) == sectors for s in joined)
    return {
        "physical_route_ports_certified": len(ports),
        "guide_edges_shared_with_wheel_loading_fixture": [
            seam for seam in route if seam in EDGES
        ],
        "vacancy_head_destinations_at_shared_start": list(ports[1:3]),
        "returned_collector_word_matches_before_handoff": True,
        "vacancy_local_loading_trace": [
            [(x, w & complement) for x, w in s] for s in joined
        ],
    }


def click(layout, basis):
    shut, launch, flags, _ = instrument(basis, layout, 0)
    live, launch_live, live_flags, _ = instrument(basis, layout, 0.5)
    assert launch == launch_live
    np.testing.assert_array_equal(flags, live_flags)
    face = np.eye(len(shut)) - 2 * (live - shut)
    np.testing.assert_array_equal(face @ face, np.eye(len(face)))
    observed = face[flags][:, launch]
    return observed.T @ observed


def run():
    archive = preparation()
    primary, alternate = primary_layout(), alternate_layout()
    assert primary["wheel_input_words"] == alternate["wheel_input_words"] == [218, 698]
    board = federation(2)
    faces = set(board.faces())
    ports_doc = json.loads(PROPAGATION.read_text())
    ports = ports_doc["ports"]
    assert ports[:3] == [A, LAND, CENTER] and len(ports) == len(set(ports)) == 144
    assert all(edge(a, b) in board.seams for a, b in pairwise(ports))
    for layout in (primary, alternate):
        for pulse in layout["native_address_schedule"]:
            assert (
                edge(*pulse) in board.seams
                if len(pulse) == 2
                else tuple(sorted(pulse)) in faces
            )
        layout["collector_compatibility"] = collector_replay(
            layout, ports, ports_doc["gate_words"]
        )
    assert primary["ordered_two_word_trace"][-1] == tuple(
        map(tuple, archive["ordered_two_word_trace"][-1])
    )
    exterior = ((1045, 1052), (1045, 1053))
    changed = sum(1 << EDGES.index(seam) for seam in exterior)
    assert all(
        z ^ w == changed
        for (_, z), (_, w) in zip(
            primary["ordered_two_word_trace"][-1],
            alternate["ordered_two_word_trace"][-1],
            strict=True,
        )
    )
    basis = wheel_basis()
    saved, first, second = (
        click(route, basis) for route in (archive, primary, alternate)
    )
    np.testing.assert_array_equal(first, saved)
    np.testing.assert_array_equal(second, np.eye(len(basis)) - saved)
    return {
        "primary_exact_archived_outside_heads": primary,
        "alternate_complementary_outside_heads": alternate,
        "physical_outside_seams": exterior,
        "unchanged_primary_exterior_head_effect": True,
        "alternate_changed_head_effect_is_complement": True,
        "basis_size": len(basis),
        "source_map_in_ordered_two_word_code": [[1, 0], [0, 1]],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    body = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as handle:
            handle.write(body)
    else:
        assert json.loads((HERE / "results.json").read_text()) == json.loads(body)
        print(
            "Verified both collector-to-wheel transfers and the two binary exterior effects."
        )


if __name__ == "__main__":
    main()
