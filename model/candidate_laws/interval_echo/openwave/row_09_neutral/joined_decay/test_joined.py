import json

import numpy as np
from scipy.linalg import eigh, expm

from ..passive_decay.measure import SHOTS, hamiltonians
from .measure import HERE, PERIOD, core_pulse, delay_kernel, run


def test_fourier_reader_matches_direct_time_order_and_keeps_all_bins():
    source, _, split, join = hamiltonians(4)
    energy, frame = eigh(split + join)
    z = np.diag(1 - 2 * (np.arange(16) % 2))
    start = np.kron(core_pulse(eigh(source)[1]), np.eye(4)) @ frame[:, 0]
    probe = frame.T @ z @ frame
    count = 16
    interval = 2 * np.pi / PERIOD
    paths = []
    for n in range(count):
        before = expm(-1j * (split + join) * n * interval)
        after = expm(-1j * (split + join) * (count - 1 - n) * interval)
        paths.append(after @ z @ before @ start)
    observed = (abs(np.fft.ifft(paths, axis=0)) ** 2).sum(axis=1)
    kernel = delay_kernel(energy[None, :] - energy[:, None], np.arange(count), count)
    filtered = (kernel * probe) @ (frame.T @ start)
    np.testing.assert_allclose((abs(filtered) ** 2).sum(axis=1), observed, atol=2e-13)
    np.testing.assert_allclose(observed.sum(), 1)


def test_positive_bin_sign_on_a_known_lowering_transition():
    # One native calibration gap equal to three bins in this short diagnostic clock.
    n = 32
    delta = np.array([[0.0, 3 * PERIOD / n], [-3 * PERIOD / n, 0.0]])
    density = abs(delay_kernel(delta, np.arange(n), n)) ** 2
    assert np.argmax(density[:, 0, 1]) == 3
    assert np.argmax(density[:, 1, 0]) == n - 3


def test_held_emission_and_cut_controls():
    record = run()
    assert record == json.loads((HERE / "records.json").read_text())
    # Hoeffding union bound over 192 instantaneous and 18 time-bin probabilities.
    error = np.sqrt(np.log(2 * 210 / 0.01) / (2 * SHOTS))
    intact, cold, cut = record["panels"]
    for a, b, c in zip(
        intact["absorption"], cold["absorption"], cut["absorption"], strict=True
    ):
        loss = a["loss_gain_other"][0] / SHOTS
        assert loss - error > intact["far_probe_light_free_bound"]
        assert (
            loss
            > max(b["loss_gain_other"][0], c["loss_gain_other"][0]) / SHOTS + 2 * error
        )
    for panel in record["panels"]:
        for read in panel["flight"]:
            assert np.sum(read["core_by_far_x"]) == SHOTS

    def far_rate(panel, i):
        return np.sum(np.asarray(panel["flight"][i]["core_by_far_x"])[:, 1]) / SHOTS

    assert abs(far_rate(intact, 0) - far_rate(cold, 0)) < 2 * error
    assert far_rate(intact, 5) > max(far_rate(cold, 5), far_rate(cut, 5)) + 2 * error
