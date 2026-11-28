"""Compare one incident face read after stationary energy selection with source 98."""

from itertools import product
from math import log, sqrt

import numpy as np

from .....fan_recoil.self_start.measure import DRAW_PAD, PROBABILITY_TICKS
from ...measure import PAULI
from ..measure import moments, partial_transpose
from .wheel import FACES, RIM, WHEEL, face_gate


def reconstructed(t):
    return (
        sum(
            t[i, j] * np.kron(a, b)
            for i, a in enumerate(PAULI)
            for j, b in enumerate(PAULI)
        )
        / 4
    )


def witness_counts(t, coeff, den, incoming, at_port, pad, rng):
    shots = 2**42
    x, y = np.array(list(product((1, -1), repeat=2))).T
    prob = t[0, 0]
    estimate = coeff[0, 0]
    total_range = 0.0
    banks = []
    for i, j in product((1, 2, 3), repeat=2):
        sharp = (prob + x * t[i, 0] + y * t[0, j] + x * y * t[i, j]) / 4
        probabilities = np.array(
            [*sharp, at_port - prob, den - at_port, incoming - den, 1 - incoming]
        )
        ticks = np.rint(probabilities * PROBABILITY_TICKS).astype(np.int64)
        ticks[ticks.argmax()] += PROBABILITY_TICKS - ticks.sum()
        np.testing.assert_allclose(
            ticks / PROBABILITY_TICKS, probabilities, atol=DRAW_PAD, rtol=0
        )
        counts = rng.multinomial(shots, ticks / PROBABILITY_TICKS).tolist()
        weights = coeff[i, j] * x * y
        if i == j:
            weights = weights + coeff[i, 0] * x + coeff[0, j] * y
        accepted = sum(counts[:4])
        estimate += counts[:4] @ weights / accepted
        total_range += np.ptp(weights) ** 2 / accepted
        banks.append(
            {"word_reference_axes": "IXYZ"[i] + "IXYZ"[j], "original_counts": counts}
        )
    # Trace(W rho) for partial transpose of a unit pure projector is contractive.
    sys = (2 * pad + 16 * DRAW_PAD * (1 + abs(coeff).sum())) / (
        prob - pad - 4 * DRAW_PAD
    )
    radius = sqrt(log(2 * 2 / 0.01) * total_range / 2) + sys
    return {
        "original_starts_per_setting": shots,
        "pauli_coefficients": coeff.tolist(),
        "predicted_conditional_score": (coeff * t).sum() / prob,
        "count_score": estimate,
        "simultaneous_99pct_interval": [estimate - radius, estimate + radius],
        "banks": banks,
    }


def read(basis, embed, project, t_input, pad, rng, counted):
    rho = reconstructed(t_input)
    transfer = project @ embed
    passed = float(
        np.trace(np.kron(transfer.conj().T @ transfer, np.eye(2)) @ rho).real
    )
    local = []
    for j, x in enumerate(RIM):
        here = transfer[[i for i, (p, _) in enumerate(basis) if p == x]]
        at_port = float(np.trace(np.kron(here.conj().T @ here, np.eye(2)) @ rho).real)
        for f in (FACES[(j - 1) % 5], FACES[j]):
            mask, cw = face_gate(WHEEL, f)
            pairs = {}
            for i, (p, w) in enumerate(basis):
                if p == x and w & mask in (cw, cw ^ mask):
                    pairs.setdefault(w & ~mask, {})[int((w & mask) != cw)] = i
            assert len(pairs) == 1 and all(set(a) == {0, 1} for a in pairs.values())
            order = next(iter(pairs.values()))
            k = transfer[[order[0], order[1]]]
            joint = np.kron(k, np.eye(2))
            result = joint @ rho @ joint.conj().T
            en = np.linalg.eigvalsh(partial_transpose(result))
            smin = np.linalg.svd(k, compute_uv=False)[-1]
            rec = {
                "charge_port": x,
                "incident_face": list(f),
                "charge_port_per_original": at_port,
                "local_cyclic_per_original": result.trace().real,
                "word_map_least_singular_value": smin,
                "minimum_partial_transpose_interval": [en[0] - pad, en[0] + pad],
                "local_word_and_source_moments": moments(result).tolist(),
            }
            if counted and x == 1046 and f == FACES[0]:
                assert smin > 1e-5
                initial = np.array([2, 1, 1, -1 - 1j])
                fitted = np.linalg.solve(joint.conj().T, initial)
                fitted /= np.linalg.norm(fitted)
                coeff = np.array(
                    [
                        [
                            np.vdot(fitted, np.kron(a, b.T) @ fitted).real / 4
                            for b in PAULI
                        ]
                        for a in PAULI
                    ]
                )
                rec["local_witness"] = witness_counts(
                    moments(result), coeff, passed, t_input[0, 0], at_port, pad, rng
                )
            local.append(rec)
    return local
