"""Forecast the two quadratures of a supplied source-side eraser and a remote pair click."""

import json
from itertools import permutations

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.apparatus import (
    native_graph,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.instruments import (
    material,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release.flight import (
    chart,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release.geometry import (
    layout,
)

TIMES = (256, 512, 2048)
OWNERS = (12, 17)


def launch():
    """Reduced matrices on the two kept direction records, per original start."""
    generator = -32j * native_graph()
    source = np.zeros((generator.shape[0], 3), dtype=complex)
    source[:8] = material("hub")["character"]
    tip = expm_multiply(generator, source, traceA=generator.diagonal().sum())[64:128]
    # Axes: third (inactive) tip bit, two active bits, untouched hub bits, recipe.
    tip = tip.reshape(2, 4, 8, 3)
    return np.einsum("rmhj,rnhj->jmn", tip, tip.conj())


def gram_by_owner():
    """Overlap at the far click for each pair of source sectors, with order unread."""
    ports = {owner: [p for p in layout()[0] if p // 12 == owner] for owner in OWNERS}
    coords = {owner: list(permutations(p, 2)) for owner, p in ports.items()}
    amplitudes = {owner: [] for owner in OWNERS}
    for mask in range(4):
        words, generator = chart((160, 167, mask))
        energies, basis = eigh(generator, driver="evd")
        initial_row = words.index((160, 167, mask))
        for owner in OWNERS:
            index = [words.index((u, v, mask)) for u, v in coords[owner]]
            phases = np.exp(-1j * np.outer(TIMES, energies))
            amplitudes[owner].append(
                (phases * basis[initial_row].conj()) @ basis[index].T
            )
    return {
        owner: np.einsum("mtp,ntp->tmn", np.array(rows), np.array(rows).conj())
        for owner, rows in amplitudes.items()
    }


def read(rho):
    """All four probabilities of each separately supplied local quadrature."""
    result = {}
    for axis, phase in (("X", 1), ("Y", 1j)):
        vectors = np.eye(4, dtype=complex)
        vectors[:, 1] = np.array([0, 1, phase, 0]) / np.sqrt(2)
        vectors[:, 2] = np.array([0, 1, -phase, 0]) / np.sqrt(2)
        np.testing.assert_allclose(vectors.conj().T @ vectors, np.eye(4), atol=1e-14)
        probabilities = np.diag(vectors.conj().T @ rho @ vectors).real
        assert min(probabilities) >= -1e-20
        np.testing.assert_allclose(probabilities.sum(), np.trace(rho).real, atol=1e-20)
        result[axis] = probabilities.tolist()
    return result


def forecast():
    density = launch()
    launch_probability = np.trace(density, axis1=1, axis2=2).real
    held = []
    for owner, gram in gram_by_owner().items():
        for time, overlap in zip(TIMES, gram, strict=True):
            outcomes = []
            for source, rho in enumerate(density):
                joint = rho * overlap
                p = np.trace(joint).real
                assert 0 < p <= launch_probability[source] < 1 / 3
                outcomes.append(
                    dict(
                        recipe=source,
                        original_refused=2 / 3,
                        initial_arrival_failed=1 / 3 - launch_probability[source],
                        far_pair_missed=launch_probability[source] - p,
                        granted_joint=read(joint),
                    )
                )
            plus, minus = (density[j] * overlap for j in (1, 2))
            p_plus, p_minus = (np.trace(a).real for a in (plus, minus))
            np.testing.assert_allclose(p_plus, p_minus, rtol=1e-9, atol=0)
            delta = plus / p_plus - minus / p_minus
            phase_contrast = [2 * delta[1, 2].real, -2 * delta[1, 2].imag]
            held.append(
                dict(
                    owner=owner,
                    time=time,
                    original_pair_probability=p_plus,
                    difference_in_conditional_exit_contrasts=phase_contrast,
                    contrast_magnitude=float(np.linalg.norm(phase_contrast)),
                    outcomes=outcomes,
                )
            )
    return dict(
        original_source_arrivals=launch_probability.tolist(),
        preflight_imaginary_coherence=[
            float((r / p)[1, 2].imag)
            for r, p in zip(density, launch_probability, strict=True)
        ],
        pointer_order=["00", "positive", "negative", "11"],
        held=held,
    )


if __name__ == "__main__":
    print(json.dumps(forecast(), indent=2))
