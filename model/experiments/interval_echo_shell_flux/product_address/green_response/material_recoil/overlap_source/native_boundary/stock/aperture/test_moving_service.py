from collections import deque
from itertools import permutations

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    moving_service as lab,
)


def test_stationary_support_equals_the_actual_monitored_orbit():
    # The first three edges form a blank connected bypass. The closing edge
    # carries one physical head. Explore every native outcome explicitly.
    half = lab.gates()[1]
    for start in ((0, 0, 1, 2), (1, 1, 1, 2), (0, 2, 0, 2)):
        pop = tuple(start.count(i) for i in range(3))
        components = lab.stationary_heads(pop)
        for h in range(11):
            seen = {(start, h)}
            pending = deque(seen)
            while pending:
                loads, mark = pending.popleft()
                for u, v in ((0, 1), (1, 2), (2, 3), (0, 3)):
                    marked = (u, v) == (0, 3)
                    col = 9 * (mark if marked else 0) + 3 * loads[u] + loads[v]
                    for out in np.flatnonzero(half[:, col]):
                        new = list(loads)
                        new[u], new[v] = divmod(int(out) % 9, 3)
                        key = (tuple(new), int(out) // 9 if marked else mark)
                        if key not in seen:
                            seen.add(key)
                            pending.append(key)
            assert seen == {
                (order, label)
                for order in permutations(start)
                for label in components[h]
            }


def test_entropy_of_dependent_gate_history_includes_the_initial_record():
    # The second gate shares a port with the first and sees its changed load.
    full, half = lab.gates()
    paths = {}
    for word in set(permutations((0, 0, 2))):
        for x in np.flatnonzero(half[:, 3 * word[0] + word[1]]):
            a, b = divmod(int(x), 3)
            p = half[x, 3 * word[0] + word[1]] / 24
            for y in np.flatnonzero(half[:, 3 * b + word[2]]):
                paths[word, int(x), int(y)] = p * half[y, 3 * b + word[2]] / 8
    probs = np.array(list(paths.values()))
    np.testing.assert_allclose(probs.sum(), 1)
    np.testing.assert_allclose(-probs @ np.log2(probs), np.log2(3) + 2 * 2 / 3)
    # A returning material word and a perfectly decoded transmitted symbol
    # are different readings even for a single pristine blank gate.
    assert full[2, 2] == 0 and np.all(lab.readers()[2][:2, :9] == 16)
