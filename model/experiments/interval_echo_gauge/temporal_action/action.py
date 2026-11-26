"""Material E/B updates and independently compiled potential-window action."""

from collections import defaultdict

import numpy as np

from .material import geometry, source

STEPS = 4
L = 2 * 4**STEPS
CASES = ("idle", "loop", "reverse", "open_quiet", "open_swap")


def dot_terms(row, operand, weight=1):
    return [(int(weight * c), int(operand[j])) for j, c in enumerate(row) if c]


def update(meter, case, source_sign=1, magnetic=1):
    d, c, seams, _, (u, _, _), hole = geometry()
    rho, j, material = source(case)
    e = np.zeros(30, dtype=int)
    start = seams.index(tuple(sorted((u, hole))))
    e[start] = 1 if seams[start] == (hole, u) else -1
    a = np.zeros(30, dtype=int)
    b = np.zeros(20, dtype=int)
    potential = []
    electric = []
    flux = []
    for n in range(STEPS + 1):
        scale = 4**n
        potential.append(a * (L // scale))
        electric.append(e * (L // scale))
        flux.append(b * (L // scale))
        # Times are in exchange-independent action units. These A, E, B
        # integers use a finer exact unit at each subsequent appointment.
        if n == STEPS:
            potential.append(
                np.array([meter.read([(2, x), (-1, y)]) for x, y in zip(a, e)])
                * (L // (2 * scale))
            )
            break
        a = np.array([meter.read([(4, x), (-2, y)]) for x, y in zip(a, e)])
        mid = np.array(
            [meter.read([(2, v), *dot_terms(row, e, -1)]) for v, row in zip(b, c)]
        )
        e = np.array(
            [
                meter.read(
                    [
                        (4, v),
                        *dot_terms(row, mid, magnetic),
                        (-2 * scale * source_sign, int(q)),
                    ]
                )
                for v, row, q in zip(e, c.T, j[n])
            ]
        )
        b = np.array([meter.read([(2, v)]) for v in mid])
    return dict(
        A=np.array(potential),
        E=np.array(electric),
        B=np.array(flux),
        phi=np.zeros((STEPS + 1, 12), dtype=int),
        J=np.vstack((j, np.zeros(30, int))) * L,
        rho=rho * L,
        record=material,
    )


def gauge(meter, hist):
    d = geometry()[0]
    # End gauge vanishes. The compact gauge reaches all 12 ports.
    chi = np.zeros((STEPS + 2, 12), dtype=int)
    for n in range(1, STEPS + 1):
        chi[n] = L * ((np.arange(12) * (n + 1) + n) % 5 - 2)
    aa = np.array(
        [
            [meter.read([(v,), *dot_terms(row, chi[n])]) for v, row in zip(h, d)]
            for n, h in enumerate(hist["A"])
        ]
    )
    phi = np.array(
        [
            [
                meter.read([(v,), (-2, nxt), (2, old)])
                for v, nxt, old in zip(row, chi[n + 1], chi[n])
            ]
            for n, row in enumerate(hist["phi"])
        ]
    )
    return {**hist, "A": aa, "phi": phi}


class Action:
    """Each local term is a signed square or an explicit source coupling.

    Varying one coordinate changes only its incident terms. Unaffected gates
    cancel between the positive and negative paths of the variation meter.
    """

    def __init__(self, meter, hist, source_sign=1, magnetic=1, kinetic=1, port_sign=1):
        self.meter = meter
        self.terms = []
        self.incident = defaultdict(list)
        d, c = geometry()[:2]
        a, phi, rho, j = (hist[k] for k in ("A", "phi", "rho", "J"))
        for n in range(STEPS + 1):
            for edge in range(30):
                operands = [
                    (("A", n + 1, edge), -2, a[n + 1, edge]),
                    (("A", n, edge), 2, a[n, edge]),
                ]
                operands += [
                    (("phi", n, u), -int(d[edge, u]), phi[n, u])
                    for u in np.flatnonzero(d[edge])
                ]
                self.square(kinetic, operands)
                self.linear(
                    2 * source_sign, j[n, edge], a[n + 1, edge], ("A", n + 1, edge)
                )
            for row in c:
                self.square(
                    -magnetic,
                    [
                        (("A", n + 1, k), int(row[k]), a[n + 1, k])
                        for k in np.flatnonzero(row)
                    ],
                )
            for u in range(12):
                self.linear(2 * port_sign, rho[n, u], phi[n, u], ("phi", n, u))

    def square(self, weight, operands):
        value = self.meter.read(
            [(coef, int(input_word)) for _, coef, input_word in operands]
        )
        self.terms.append((weight, value, value))
        for key, coef, _ in operands:
            self.incident[key].append((weight, value, coef, True))

    def linear(self, weight, source, operand, key):
        self.terms.append((weight, int(source), int(operand)))
        self.incident[key].append((weight, int(source), 1, False))

    def value(self):
        return self.meter.read(self.terms)

    def directions(self):
        out = {}
        for tag, steps, dim in (
            ("A", range(1, STEPS + 1), 30),
            ("phi", range(STEPS + 1), 12),
        ):
            result = []
            for n in steps:
                row = []
                for i in range(dim):
                    terms = []
                    for weight, value, delta, squared in self.incident[tag, n, i]:
                        if squared:
                            plus = self.meter.read([(value,), (delta,)])
                            minus = self.meter.read([(value,), (-delta,)])
                            terms += [(weight, plus, plus), (-weight, minus, minus)]
                        else:
                            terms.append((2 * weight, value, delta))
                    row.append(self.meter.read(terms))
                result.append(row)
            out[tag] = result
        return out


def independent(hist, source_sign=1, magnetic=1, kinetic=1, port_sign=1):
    """From potential and source coordinates, evaluate the full integer action."""
    d, c = geometry()[:2]
    a, phi, rho, j = (hist[k] for k in ("A", "phi", "rho", "J"))
    electric = -2 * np.diff(a, axis=0) - phi @ d.T
    flux = a[1:] @ c.T
    return int(
        kinetic * np.sum(electric**2)
        - magnetic * np.sum(flux**2)
        + 2 * source_sign * np.sum(j * a[1:])
        + 2 * port_sign * np.sum(rho * phi)
    )


def independently_vary(hist, key, **weights):
    outputs = []
    steps = range(1, STEPS + 1) if key == "A" else range(STEPS + 1)
    for n in steps:
        row = []
        for i in range(hist[key].shape[1]):
            a = hist[key].copy()
            a[n, i] += 1
            plus = independent({**hist, key: a}, **weights)
            a[n, i] -= 2
            minus = independent({**hist, key: a}, **weights)
            row.append(plus - minus)
        outputs.append(row)
    return outputs
