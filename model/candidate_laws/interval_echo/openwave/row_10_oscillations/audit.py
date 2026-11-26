"""Finite checks behind the conjugate-source ambiguity on the neutral release.

Run from the repository root using python -m on this module. The checks use
the native contact matrices and actual source pulses. No oscillation frequency,
transition fit, propagation time, remote counts or mixing matrix is an input.
"""

import json
from functools import cache
from itertools import permutations

import numpy as np

from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    LAW,
    blocks,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.apparatus import (
    fixture,
    native_graph,
    permutations_of_chart,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.instruments import (
    material,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release.geometry import (
    layout,
)


def permute_word(word, permutation):
    return sum(((word >> i) & 1) << j for i, j in enumerate(permutation))


def source_reflections():
    """For each sharp three-bit tip word, exhibit a reflection fixing its effect."""
    chart = permutations_of_chart()
    source = material("hub")["character"]
    target = np.exp(2j * np.pi * np.arange(3) / 3) / 3
    for character, column in ((1, target), (2, target.conj())):
        ideal = np.zeros(8, complex)
        ideal[[1, 2, 4]] = column
        phase = np.vdot(ideal, source[:, character]) * 3
        np.testing.assert_allclose(abs(phase), 1, rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            source[:, character], phase * ideal, rtol=0, atol=1e-12
        )

    charged = native_graph(-1)
    assert (charged != native_graph(1)).nnz == 0
    assert np.array_equal(8 * charged.data, np.rint(8 * charged.data))
    for perm, index in chart.items():
        assert (charged[index][:, index] != charged).nnz == 0
        # Verify the complete hub-head sum at the physical receiving tip.
        tip = fixture()["ports"].index(fixture()["tip"])
        for word in range(8):
            old = [64 * tip + 8 * word + h for h in range(8)]
            new = {64 * tip + 8 * permute_word(word, perm) + h for h in range(8)}
            assert set(index[old]) == new
    stabilizers = {}
    for word in range(8):
        reflection = next(
            perm
            for perm in chart
            if sum(i != j for i, j in enumerate(perm)) == 2
            and permute_word(word, perm) == word
        )
        # Write the relative winding phases as exponents modulo three. Every
        # reflection sends j to c-j, hence conjugates the incoming character.
        assert len({(j + reflection[j]) % 3 for j in range(3)}) == 1
        stabilizers[str(word)] = list(reflection)
    return dict(
        charged_chart_dimension=charged.shape[0],
        exact_commuting_arm_permutations=len(chart),
        sharp_tip_word_stabilizers=stabilizers,
    )


@cache
def local_moves(marked, reversed_mark, loads, head):
    """Read the exact off-diagonal possibilities directly from the integer 8H."""
    matrix = blocks()[(1, 1) if marked else "blank"]
    native_head = head ^ reversed_mark if marked else 0
    column = 9 * native_head + LAW.SRC.index(loads)
    moves = []
    for row in np.flatnonzero(matrix[:, column]):
        if row == column:
            continue
        new_loads = LAW.SRC[row % 9]
        new_head = int(row // 9) ^ reversed_mark if marked else 0
        moves.append((new_loads, new_head))
    return moves


def release_graph():
    """Build all ordered vacancy/excess words and both sharp head bits."""
    ring, marks, leaves, _, _ = layout()
    sites = tuple(ring) + tuple(leaves)
    contacts = [(u, v, None) for u, v in zip(ring, ring[1:] + ring[:1], strict=True)]
    contacts.extend((u, v, j) for j, (u, v) in enumerate(marks))
    graph = {}
    for vacancy, excess in permutations(sites, 2):
        for mask in range(4):
            adjacent = set()
            for u, v, j in contacts:
                a, b = sorted((u, v))
                before = tuple(
                    0 if p == vacancy else 2 if p == excess else 1 for p in (a, b)
                )
                head = (mask >> j) & 1 if j is not None else 0
                for after, new_head in local_moves(
                    j is not None, int(u > v), before, head
                ):
                    next_vacancy = (a, b)[after.index(0)] if 0 in after else vacancy
                    next_excess = (a, b)[after.index(2)] if 2 in after else excess
                    next_mask = (
                        mask ^ ((new_head ^ head) << j) if j is not None else mask
                    )
                    adjacent.add((next_vacancy, next_excess, next_mask))
            graph[vacancy, excess, mask] = adjacent
    for word, neighbors in graph.items():
        assert all(word in graph[n] for n in neighbors)
    return graph


def component(graph, start):
    reached = {start}
    queue = [start]
    for current in queue:
        for word in graph[current] - reached:
            reached.add(word)
            queue.append(word)
    return reached


def environment(word, leaves):
    vacancy, excess, mask = word
    loads = tuple(0 if p == vacancy else 2 if p == excess else 1 for p in leaves)
    return loads, mask


def sector_from_leaf_read(word, leaves):
    loads, mask = environment(word, leaves)
    changed = sum(int(load != 1) << i for i, load in enumerate(loads))
    return mask ^ changed


def compare_contact_support():
    """A contrast change affects exchange amplitudes, while retaining the same support."""
    nonzero = []
    for kind in ("blank", (1, 1)):
        support = []
        for strength in (4, 8):
            off_diagonal = blocks(strength)[kind] != 0
            np.fill_diagonal(off_diagonal, False)
            support.append(off_diagonal)
        np.testing.assert_array_equal(*support)
        nonzero.append(int(np.count_nonzero(support[0])))
    return nonzero


def release_sectors():
    graph = release_graph()
    leaves = layout()[2]
    for word, neighbors in graph.items():
        assert all(
            sector_from_leaf_read(word, leaves) == sector_from_leaf_read(other, leaves)
            for other in neighbors
        )
    result = {}
    for sign in (-1, 1):
        sites = (160, 167) if sign == -1 else (167, 160)
        keys = {}
        sizes = []
        for launch_mask in range(4):
            reachable = component(graph, (*sites, launch_mask))
            sizes.append(len(reachable))
            for word in reachable:
                key = environment(word, leaves)
                assert sector_from_leaf_read(word, leaves) == launch_mask
                if key in keys:
                    assert keys[key] == launch_mask
                keys[key] = launch_mask
        result[str(sign)] = dict(
            reachable_words_by_mask=sizes, environment_keys=len(keys)
        )
    return dict(
        full_release_words=len(graph),
        compared_contrasts=[4, 8],
        off_diagonal_contact_entries=compare_contact_support(),
        directed_nonzero_exchanges=sum(map(len, graph.values())),
        active_leaf_ports=list(leaves),
        launch_sectors=result,
    )


def audit():
    return dict(source=source_reflections(), release=release_sectors())


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
