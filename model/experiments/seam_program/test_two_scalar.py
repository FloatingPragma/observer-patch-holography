"""Independent checks of full binary pulses, logical two-packet evolution and spatial-band dilution."""

from itertools import product
import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from .two_scalar import ARMS, HERE, run


def test_native_resonant_binary_gate_and_contact_control():
    indices = [law.SRC.index(pair) for pair in ((0, 0), (0, 1), (1, 0), (1, 1))]
    inputs = np.eye(law.DIM)[:, indices]
    exchange, contrast = law.local_term_parts()
    swap = np.eye(4)[:, [0, 2, 1, 3]]
    for q in (1, 4):
        t = q * np.pi / 4
        expected = inputs @ expm(-1j * t * (np.eye(4) - swap) / 2)
        for kappa in (0, 8):
            np.testing.assert_allclose(
                expm(-1j * t * (exchange + kappa * contrast)) @ inputs,
                expected,
                atol=3e-14,
            )
    extra = np.diag([0, 0, 0, 1])
    changed = expm(-1j * np.pi / 4 * ((np.eye(4) - swap) / 2 + extra))
    native = expm(-1j * np.pi / 4 * (np.eye(4) - swap) / 2)
    np.testing.assert_allclose(changed[:, :3], native[:, :3], atol=1e-14)
    np.testing.assert_allclose(
        changed[3, 3], np.exp(-1j * np.pi / 4) * native[3, 3], atol=1e-14
    )


def logical_joint(launch, contact_phase, steps):
    # Ordered, symmetric first-quantized matrix. The excluded diagonal always stays zero.
    slots = np.arange(192).reshape(2, 2, 2, 3, 8)
    stations = slots.reshape(8, 24)
    a, b = [stations[i] for i in launch]
    state = np.zeros((192, 192), dtype=complex)
    state[np.ix_(a, b)] = state[np.ix_(b, a)] = 1 / (24 * np.sqrt(2))
    permutation = np.arange(192)
    for address in np.ndindex(2, 2, 2, 3):
        for axis, length in enumerate(slots.shape[:-1]):
            target = list(address)
            target[axis] = (target[axis] + 1) % length
            i, j = slots[address][2 * axis], slots[tuple(target)][2 * axis + 1]
            permutation[i], permutation[j] = j, i
    z = np.exp(-1j * np.pi / 4)
    for _ in range(steps):
        for block in slots.reshape(-1, 8):
            for mask in range(1, 8):
                for i in range(8):
                    if i >= (i ^ mask):
                        continue
                    a, b = block[i], block[i ^ mask]
                    both = state[a, b]
                    va, vb = state[[a, b]]
                    state[a] = state[:, a] = ((1 + z) * va + (1 - z) * vb) / 2
                    state[b] = state[:, b] = ((1 - z) * va + (1 + z) * vb) / 2
                    state[a, a] = state[b, b] = 0
                    state[a, b] = state[b, a] = contact_phase * both
        state = state[np.ix_(permutation, permutation)]
    mass = (abs(state) ** 2).reshape(8, 24, 8, 24).sum(axis=(1, 3))
    np.testing.assert_allclose(mass.sum(), 1, atol=1e-12)
    return [
        mass[i, j] if i == j else mass[i, j] + mass[j, i]
        for i in range(8)
        for j in range(i, 8)
    ]


def test_complete_two_mobile_records_and_independent_logical_contact():
    saved = json.loads((HERE / "two_scalar.json").read_text())
    fresh = run()
    assert fresh["protocol_sha256"] == saved["protocol_sha256"]
    assert tuple(saved["arms"]) == ARMS
    assert saved["labels"] == [
        list(x) if isinstance(x, tuple) else x for x in fresh["labels"]
    ]
    periods = np.array(saved["period_probabilities"])
    np.testing.assert_allclose(
        fresh["period_probabilities"], periods, atol=2e-12, rtol=0
    )
    np.testing.assert_allclose(periods.sum(axis=-1), 1, atol=2e-11, rtol=0)
    assert np.all(periods[..., -1] == 0)
    for j, launch in enumerate(saved["launches"]):
        for arm, phase in enumerate((1, np.exp(-1j * np.pi / 4))):
            np.testing.assert_allclose(
                logical_joint(launch, phase, 3),
                periods[3, j, arm, :-1],
                rtol=0,
                atol=2e-12,
            )
    for new, old in zip(fresh["stopped"], saved["stopped"], strict=True):
        assert new["clock"] == old["clock"] and list(new["launch"]) == old["launch"]
        np.testing.assert_allclose(new["laws"], old["laws"], atol=2e-12, rtol=0)
        distributions = np.array(old["laws"])
        np.testing.assert_allclose(distributions.sum(axis=1), 1, atol=2e-11, rtol=0)
        assert np.all(distributions[:, -1] == distributions[0, -1])
        for arm, key in ((1, "tv_contact"), (2, "tv_free")):
            np.testing.assert_allclose(
                abs(distributions[0] - distributions[arm]).sum() / 2,
                old[key],
                atol=2e-11,
                rtol=0,
            )


def test_collision_bound_attained_by_correlated_bandlimited_pair():
    modes = np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0)])
    for n in (4, 8):
        sites = np.array(list(product(range(n), repeat=3)))
        waves = np.exp(2j * np.pi * sites @ modes.T / n) / np.sqrt(n**3)
        state = waves @ waves.conj().T / np.sqrt(len(modes))
        np.testing.assert_allclose(np.linalg.norm(state), 1, atol=1e-14)
        np.testing.assert_allclose(
            np.sum(abs(np.diag(state)) ** 2), len(modes) / n**3, rtol=1e-14
        )
