"""Read an outside head through the original adjoining face at 1045."""

from collections import deque
import numpy as np

from .wheel import EDGES, FACES, RIM, WHEEL, edge, face_gate, flip, head, sides


OUT = (1045, 1052, 1053)


def instrument(basis, route, rate):
    extra = sorted({edge(a, b) for a, b in sides(OUT)} - set(WHEEL))
    es = [*WHEEL, *extra]
    end = route["ordered_two_word_trace"][-1]

    def outside(w):
        return sum(
            int(head(e, w, EDGES.index(e)) == e[0]) << (10 + j)
            for j, e in enumerate(extra)
        )

    ext = outside(end[0][1])
    assert ext == outside(end[1][1])
    gates = [face_gate(es, f) for f in (*FACES, OUT)]

    def transitions(x, w):
        for i, (u, v) in enumerate(sides(RIM)):
            e = edge(u, v)
            j = es.index(e)
            if x in e and head(e, w, j) == x:
                yield (v if x == u else u, w ^ (1 << j)), 4
        for i, gate in enumerate(gates):
            z = flip(w, gate)
            if z != w:
                yield (x, z), 0.5 if i < 5 else rate

    embed_states = [(x, w | ext) for x, w in basis]
    seen = set(embed_states)
    todo = deque(embed_states)
    while todo:
        for row, _ in transitions(*todo.popleft()):
            if row not in seen:
                seen.add(row)
                todo.append(row)
    states = sorted(seen)
    ids = {row: i for i, row in enumerate(states)}
    h = np.diag(np.full(len(states), 8.0))
    for i, row in enumerate(states):
        for new, r in transitions(*row):
            h[i, i] += r
            h[ids[new], i] -= r
    np.testing.assert_array_equal(h, h.T)
    embed = [ids[z] for z in embed_states]
    changed = np.array([w >> 10 != ext >> 10 for _, w in states])
    assert set(w >> 10 for _, w in states) == {ext >> 10, 3 ^ (ext >> 10)}
    assert len({(x, w & 1023) for x, w in states}) == len(
        states
    )  # outside words mark distinct inside sectors
    return h, embed, changed, extra


def response(basis, route, p, prepared, extrema):
    cases = []
    for rate in (0, 0.25, 0.5):
        h, embed, changed, extra = instrument(basis, route, rate)
        en, v = np.linalg.eigh(h)
        assert np.linalg.norm(h @ v - v * en, ord=2) < 1e-12
        dwell = []
        for t in (0, 0.25, 1, 4, 8):
            u = (v * np.exp(-1j * t * en)) @ v.T[:, embed]
            probe = u[changed].conj().T @ u[changed]
            same = p @ u[embed]
            intact = same.conj().T @ same
            # Requiring the original doublet and a changed outside word is
            # impossible: the changed outside word fixes a different wheel Gauss sector.
            rows = {}
            for key, rho in {**prepared, **extrema}.items():
                passed = np.trace(rho).real
                detected = np.trace(probe @ rho).real
                keep = np.trace(intact @ rho).real
                rows[key] = {
                    "source_admission": passed,
                    "outside_reversed": detected,
                    "outside_reversed_given_source": detected / passed,
                    "original_energy_and_outside_word_given_source": keep / passed,
                    "other_source_energy_with_outside_unchanged": passed
                    - keep
                    - detected,
                }
                if rate == 0:
                    assert abs(detected) < 1e-13 and abs(keep - passed) < 1e-12
            dwell.append({"wait": t, "records": rows})
        cases.append(
            {
                "exterior_face_rate": rate,
                "outside_edges": [list(e) for e in extra],
                "joint_configuration_count": len(h),
                "dwells": dwell,
            }
        )
    return cases
