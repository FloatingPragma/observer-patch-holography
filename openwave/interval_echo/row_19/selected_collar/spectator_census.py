"""Check the material source law with mobile, densely occupied and coherent spectators."""

from itertools import product

import numpy as np

from . import contact_census as experiment

SEED = 28122912


def column(key, edges, degrees):
    diag = 16 * sum((6 - int(d)) * q * q for d, q in zip(degrees, key, strict=True))
    out = {}
    for x, y in edges:
        diff = key[x] - key[y]
        if diff:
            diag += 1 + 16 * diff * diff
            k = list(key)
            k[x], k[y] = k[y], k[x]
            k = tuple(k)
            out[k] = out.get(k, 0) - 1
    out[key] = diag
    return out


def projected(baths, q, edges, degrees, ab, inside):
    left, right = [], []
    for bath in baths:
        for g, store in zip(ab, (left, right), strict=True):
            prep = list(bath)
            prep[g] = q
            k = tuple(prep)
            store.append((k, column(k, edges, degrees)))

    def occupied(k):
        return sum(k[j] != 0 for j in inside)

    n = len(baths)
    acceleration = np.zeros((n, n), int)
    energy = np.zeros((n, n), int)
    for i, j in np.ndindex(n, n):
        u, hu = left[i]
        v, hv = right[j]
        common = hu.keys() & hv.keys()
        acceleration[i, j] = sum(
            hu[w] * hv[w] * (2 * occupied(w) - occupied(u) - occupied(v))
            for w in common
        )
        energy[i, j] = left[i][1].get(left[j][0], 0) - right[i][1].get(right[j][0], 0)
    return acceleration, energy


def grades(bath, ab, edges, deg):
    vals = []
    for g in ab:
        near = [a + b - g for a, b in edges if g in (a, b)]
        vals.append(
            np.array(
                [
                    int(deg[g]) - len(near),
                    sum(bath[j] == 1 for j in near),
                    sum(bath[j] == 2 for j in near),
                ]
            )
        )
    return vals[0] - vals[1]


def run():
    setup, *_ = experiment.prepare()
    _, _, ports, ab, _, degree, edges, second = setup
    inside = [i for i, g in enumerate(ports) if g // 12 == ports[ab[0]] // 12]
    rng = np.random.default_rng(SEED)
    cases = elements = offdiagonal = 0
    examples = []
    for flags in product(range(2), repeat=2):
        active = [
            e
            for e in edges
            if not any(z and e in w for z, w in zip(flags, second, strict=True))
        ]
        for s in range(16):
            bath = rng.integers(3, size=len(ports))
            if s < 3:
                bath[:] = s
            bath[list(ab)] = 0
            sharp = [tuple(map(int, bath))]
            for i, j in active:
                if i in ab or j in ab or bath[i] == bath[j]:
                    continue
                other = bath.copy()
                other[i], other[j] = bath[j], bath[i]
                sharp.append(tuple(map(int, other)))
                if len(sharp) >= 6:
                    break
            rows = np.array([grades(b, ab, active, degree) for b in sharp])
            energies = []
            for q in (1, 2):
                A, E = projected(sharp, q, active, degree, ab, inside)
                want = [
                    int(
                        4
                        * sum(
                            c * int(x)
                            for c, x in zip(
                                experiment.COEFFICIENTS[q - 1], d, strict=True
                            )
                        )
                    )
                    for d in rows
                ]
                np.testing.assert_array_equal(A, np.diag(want))
                np.testing.assert_array_equal(E, -A)
                energies.append(E)
                # Apply complex relative phases to a superposition of different
                # interacting spectator configurations with the same populations.
                for phase in (1, 1j, -1):
                    amp = np.array([phase**j for j in range(len(sharp))]) / np.sqrt(
                        len(sharp)
                    )
                    assert abs(np.vdot(amp, A @ amp) - np.mean(want)) < 2e-12
                elements += A.size
                offdiagonal += A.size - len(sharp)
            want = [
                int(
                    2
                    * sum(
                        c * int(x)
                        for c, x in zip(experiment.COEFFICIENTS[2], d, strict=True)
                    )
                )
                for d in rows
            ]
            np.testing.assert_array_equal(energies[1] - energies[0], np.diag(want))
            cases += 1
            if s == 4:
                examples.append(
                    {
                        "marks": list(flags),
                        "sharp_spectators": [list(v) for v in sharp],
                        "measurable_gradients": rows.tolist(),
                        "twice_clock_diagonal": want,
                    }
                )
    return {
        "seed": SEED,
        "ensembles": cases,
        "crossing_matrix_elements_checked": elements,
        "offdiagonal_elements": offdiagonal,
        "examples": examples,
        "inputs": experiment.save_inputs("CONTACT_CENSUS.md", "routed_contacts.json"),
    }


if __name__ == "__main__":
    result = run()
    experiment.emit("spectator_census", result)
    print(
        result["ensembles"],
        result["crossing_matrix_elements_checked"],
        result["offdiagonal_elements"],
    )
