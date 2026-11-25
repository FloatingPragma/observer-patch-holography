import json

import numpy as np

from .measure import HERE, SHOTS, geometry, law, local_reader, relocate, score


def test_every_sharp_input_against_local_spectrum():
    for kappa in (0, 8):
        table = local_reader(kappa)
        for a, b in law.SRC:
            expected = np.zeros(64)
            base = kappa * (a - b) ** 2
            if a == b:
                expected[0] = 1
            else:
                expected[base : base + 2] = 0.5
            np.testing.assert_allclose(table[0, a, b], expected, atol=2e-13)
            np.testing.assert_allclose(table[1, a, b], np.eye(64)[0], atol=2e-13)


def test_complete_work_and_controls_from_independent_sharp_formula():
    for label in ("records", "private"):
        for row in json.loads((HERE / (label + ".json")).read_bytes())["rows"]:
            cells = row["cells"]
            level = (20, 80, 320).index(cells)
            car, edges, external, _, gates, old, targets = geometry(level)
            new, k = row["new"], row["kappa"]
            rng = np.random.default_rng(221037 + level)
            loads = rng.integers(0, 3, car.n)
            collars = rng.integers(0, 2, len(edges))
            for owner, load in zip((old, *targets), (2, 0, 0), strict=True):
                loads[12 * owner : 12 * (owner + 1)] = load
            changed = np.flatnonzero(gates[old] != gates[new])
            if label == "private":
                for owner in (old, new):
                    loads[12 * owner + changed] = 2 * gates[owner, changed]
            second = relocate(edges, loads, collars, old, new)
            for arr, restored in zip(
                (loads, collars), relocate(edges, *second, old, new), strict=True
            ):
                np.testing.assert_array_equal(arr, restored)

            def energy(state, memory):
                a, b = state[edges].T
                delta = (k * (a - b) ** 2 + 0.5 * (a != b)) * (1 - memory)
                ref = k * (state[car.degree == 5] ** 2).sum()
                return np.array([delta.sum() + ref, delta[external].sum(), ref])

            works = energy(*second) - energy(loads, collars)
            for name, value in zip(
                ("whole", "exterior", "private"), works, strict=True
            ):
                work, err = row["work"][name]
                assert abs(work - value) <= err
            assert works[0] == works[1] + works[2]
            assert works[2] == (4 * k * len(changed) if label == "private" else 0)
            assert abs(works[0]) <= row["relocation_bound"] <= 390
            sc = score(row["banks"], k)
            assert abs(sc["boundary_residual"]) <= sc["joint_allowance"]
            for bank in row["banks"]:
                for instrument in bank.values():
                    assert sum(instrument["counts"]) == (
                        SHOTS if instrument["terms"] else 0
                    )


def test_soft_source_bound_with_noncommuting_fields():
    def ground(a):
        return np.linalg.eigvalsh(a)[0]

    def width(a):
        return np.ptp(np.linalg.eigvalsh(a))

    rng = np.random.default_rng(3874)
    for _ in range(30):
        mats = rng.normal(size=(3, 7, 7)) + 1j * rng.normal(size=(3, 7, 7))
        h, a, b = mats + mats.conj().transpose(0, 2, 1)
        connected = ground(h + a + b) - ground(h + a) - ground(h + b) + ground(h)
        assert abs(connected) <= min(width(a), width(b)) + 1e-11
    a = np.diag([0.0, 33.0])
    b = a[::-1, ::-1]
    assert ground(a + b) - ground(a) - ground(b) == min(width(a), width(b))
