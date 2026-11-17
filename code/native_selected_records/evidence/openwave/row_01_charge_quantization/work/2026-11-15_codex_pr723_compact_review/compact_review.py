"""Standalone verification of known finite row-1 results, not a physical award.

Only standard-library containers, enumeration and JSON formatting are used.
The program does not read files, use the network, launch processes, or import
the producer package. Source states below are explicitly full-state audits.
"""
from collections import deque
from itertools import combinations
import json

N = 12
SEAMS = ((0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5),
         (1, 7), (1, 8), (1, 9), (2, 3), (2, 4), (2, 6),
         (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9),
         (4, 5), (4, 9), (4, 11), (5, 9), (5, 11), (6, 7),
         (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11))


def step(state, edge, coin, law="native"):
    """Disclosed endpoint rule; current is oriented from edge[0] to edge[1]."""
    a, b = edge
    out = list(state)
    if law == "native":
        total = state[a] + state[b]
        low, high = total // 2, total - total // 2
        out[a], out[b] = (low, high) if coin == 0 else (high, low)
    else:
        assert law == "swap_hold"
        if coin == 1:
            out[a], out[b] = state[b], state[a]
    return tuple(out), state[a] - out[a]


def census():
    counts = dict(states=0, native_branches=0, binary_states=0,
                  binary_exit_branches=0, unreachable_states=0,
                  max_shortest_entry_events=0)
    for total in range(5):
        # Stars and bars enumerate every nonnegative twelve-load state once.
        states = []
        for bars in combinations(range(total + N - 1), N - 1):
            cuts = (-1,) + bars + (total + N - 1,)
            states.append(tuple(cuts[i + 1] - cuts[i] - 1 for i in range(N)))
        index = {state: i for i, state in enumerate(states)}
        assert len(index) == len(states)
        predecessors = [[] for _ in states]
        balanced = [max(state) <= 1 for state in states]
        for i, state in enumerate(states):
            for edge in SEAMS:
                for coin in (0, 1):
                    out, current = step(state, edge, coin)
                    assert sum(out) == total and min(out) >= 0
                    assert current == out[edge[1]] - state[edge[1]]
                    target = index[out]
                    predecessors[target].append(i)
                    counts["native_branches"] += 1
                    counts["binary_exit_branches"] += balanced[i] and max(out) > 1
        distance = [0 if yes else -1 for yes in balanced]
        queue = deque(i for i, yes in enumerate(balanced) if yes)
        while queue:
            target = queue.popleft()
            for previous in predecessors[target]:
                if distance[previous] == -1:
                    distance[previous] = distance[target] + 1
                    queue.append(previous)
        counts["states"] += len(states)
        counts["binary_states"] += sum(balanced)
        counts["unreachable_states"] += distance.count(-1)
        counts["max_shortest_entry_events"] = max(
            counts["max_shortest_entry_events"], max(distance))
    assert counts == dict(states=1820, native_branches=109200, binary_states=794,
                          binary_exit_branches=0, unreachable_states=0,
                          max_shortest_entry_events=3)
    return counts


def controls():
    entries = {"native": 0, "swap_hold": 0}
    samples = {}
    for port in range(N):
        initial = tuple(2 if i == port else 0 for i in range(N))
        for edge in SEAMS:
            for coin in (0, 1):
                for law in entries:
                    out, current = step(initial, edge, coin, law)
                    entries[law] += max(out) <= 1
                    if law == "swap_hold":
                        assert sorted(out) == sorted(initial)
                    if port == 0 and edge == (0, 1) and coin == 1:
                        samples[law] = dict(initial=initial, edge=edge, coin=coin,
                                            current=current, final=out)
    assert entries == {"native": 120, "swap_hold": 0}
    uniform_two = (2,) * N
    fixed = 0
    for edge in SEAMS:
        for coin in (0, 1):
            out, current = step(uniform_two, edge, coin)
            assert out == uniform_two and current == 0
            fixed += 1
    return dict(double_branches_per_law=720, double_binary_entries=entries,
                double_example_full_state_audit=samples,
                uniform_two_fixed_branches=fixed, uniform_two_is_binary=False)


if __name__ == "__main__":
    assert len(SEAMS) == len(set(SEAMS)) == 30
    assert all(sum(port in edge for edge in SEAMS) == 5 for port in range(N))
    result = dict(verification_only=True, same_group=True,
                  census=census(), controls=controls(),
                  original_long_traces_verified=False,
                  current_only_decoder_verified=False,
                  physical_charge_quantization_established=False,
                  G6_awarded=False, G9_awarded=False, passed=True)
    print(json.dumps(result, indent=2, sort_keys=True))
