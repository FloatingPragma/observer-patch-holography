"""Calibrate the applied hopping phase by opposite physical rim circuits."""

from collections import deque
from math import asin, ceil, floor

import numpy as np
from scipy.linalg import expm

from .....fan_recoil.self_start.measure import bound, draw
from .wheel import FACES, RIM, WHEEL, edge, face_gate, flip, head, sides


def reported_interval(lo, hi):
    """Outward display precision for count intervals from large binomial banks."""
    return [floor(lo * 10**8) / 10**8, ceil(hi * 10**8) / 10**8]


def tour(basis, start):
    gates = [face_gate(WHEEL, f) for f in FACES]
    queue = deque([((start, 0), ())])
    seen = {(start, 0)}
    while queue:
        (state, n), schedule = queue.popleft()
        if n == 5 and state == start:
            return schedule
        p, w = basis[state]
        neighbors = []
        if n < 5:
            i = RIM.index(p)
            dest = RIM[(i + 1) % 5]
            e = edge(p, dest)
            b = WHEEL.index(e)
            if head(e, w, b) == p:
                neighbors.append((basis.index((dest, w ^ (1 << b))), n + 1, i))
        for i, gate in enumerate(gates):
            v = flip(w, gate)
            if v != w:
                neighbors.append((basis.index((p, v)), n, 5 + i))
        for dest, count, pulse in neighbors:
            if (dest, count) not in seen:
                seen.add((dest, count))
                queue.append(((dest, count), (*schedule, pulse)))
    raise ValueError("A sharp incoming word has no restorable positive rim circuit")


def pulse_matrix(basis, i, phase):
    term = np.zeros((len(basis), len(basis)), complex)
    for j, (p, w) in enumerate(basis):
        if i < 5:
            u, v = sides(RIM)[i]
            e = edge(u, v)
            b = WHEEL.index(e)
            if p not in e or head(e, w, b) != p:
                continue
            dest = v if p == u else u
            new = (dest, w ^ (1 << b))
            z = np.exp(1j * (1 if p == u else -1) * phase / 5)
        else:
            word = flip(w, face_gate(WHEEL, FACES[i - 5]))
            if word == w:
                continue
            new = (p, word)
            z = 1
        term[basis.index(new), j] = -z
        term[j, j] = 1
    np.testing.assert_allclose(term, term.conj().T, rtol=0, atol=2e-15)
    return expm(-0.5j * np.pi * term)


def ruler(basis, start, prob, face, shots, pad, rng):
    schedule = tour(basis, start)
    runs = []
    for phi in (-1 / 32, -1 / 64, 0, 1 / 64, 1 / 32):
        u = [pulse_matrix(basis, i, phi) for i in range(10)]
        arms = []
        for order in (schedule, schedule[::-1]):
            a = np.eye(len(basis), dtype=complex)[:, start]
            for i in order:
                a = u[i] @ a
            assert np.linalg.norm(np.delete(a, start)) < 2e-13
            arms.append(a[start])
        fringe = arms[0] * arms[1].conjugate()
        np.testing.assert_allclose(fringe, np.exp(2j * phi), atol=2e-13, rtol=0)
        banks = []
        for v in (fringe.real, fringe.imag):
            bright = prob * (1 + v) / 2
            c = draw([bright, prob - bright, face - prob, 1 - face], shots, rng)
            bright_interval = bound(c[0], shots, pad)
            pass_interval = bound(sum(c[:2]), shots, pad)
            banks.append(
                {
                    "bright_dark_otherword_face_refused": c,
                    "quadrature_bounds": reported_interval(
                        max(-1, 2 * bright_interval[0] / pass_interval[1] - 1),
                        min(1, 2 * bright_interval[1] / pass_interval[0] - 1),
                    ),
                }
            )
        assert banks[0]["quadrature_bounds"][0] > 0
        iv = [asin(v) / 2 for v in banks[1]["quadrature_bounds"]]
        runs.append(
            {
                "given_phase": phi,
                "measured_phase_interval": reported_interval(*iv),
                "counts": banks,
            }
        )
    return {
        "calibration_word_probability_per_original": prob,
        "forward_pulse_addresses": [
            list(sides(RIM)[i] if i < 5 else FACES[i - 5]) for i in schedule
        ],
        "uncertainty_uses_small_flux_branch": True,
        "measurements": runs,
    }


def rotor_check(modes, flux, ruler_data):
    e0, e1, e2 = [modes[k][0][0] for k in (0, 1, 2)]
    gap = e1 - e0
    held = e2 - e0
    split = abs(np.diff(flux["phase_probes"][-1]["energies_plus_one_minus_one"])[0])
    phase = next(
        r["measured_phase_interval"]
        for r in ruler_data["measurements"]
        if r["given_phase"] == 1 / 64
    )
    # For E_k(f) = E0 + gap*(k+f)^2, f=q*flux/(2π),
    # the first pair splits by (2*gap/π)*flux.
    coef = [
        np.pi * (split + s * 2e-11) / (2 * (gap - s * 2e-11) * v)
        for s, v in ((-1, phase[1]), (1, phase[0]))
    ]
    return {
        "first_zero_field_gap": gap,
        "second_zero_field_gap": held,
        "rigid_rotor_predicted_second_gap": 4 * gap,
        "first_gap_flux_coefficient_interval": reported_interval(*coef),
    }
