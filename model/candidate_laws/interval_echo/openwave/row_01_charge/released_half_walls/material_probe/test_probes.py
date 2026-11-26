import json

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from .. import ring
from . import flight, live, measure, readback


def archived(name):
    return json.loads((measure.HERE / name).read_text())


def test_phase_and_physical_mount():
    flight.verify_plates()
    for lvl, n in flight.BOARDS:
        ports = ring.geometry(lvl, n)
        stations, remote = flight.mount(lvl, n)
        word = np.arange(n) % 3
        for m, far in zip(stations, remote, strict=True):
            j, a, b = m
            for s in (-1, 1):
                assert flight.sensed_sum(lvl, n, (a, b), word, s, 8) == -16 * s * (
                    word[j] + word[(j + 1) % n]
                )
                assert flight.sensed_sum(lvl, n, far, word, s, 8) == 0
            assert all(p not in ports for p in (a, b))


def test_complete_records_and_discriminants():
    raw = archived("counts.json")
    stats = readback.analyze(raw)
    assert stats["calibration_failures"] == 0
    assert stats["counts"]["native"]["full_match"] == 256
    assert stats["counts"]["native"]["half_detuning_excluded"] == 161
    for name in ("distant", "uncoupled"):
        assert stats["counts"][name]["zero_match"] == 256
    raw["measurements"].append(raw["measurements"][0])
    with pytest.raises(ValueError, match="duplicate"):
        readback.analyze(raw)


def test_stopped_live_generator_reproduces_two_separate_material_probabilities():
    k, src, probe, t = 8, -1, -1, 1 / 16
    mounts, _ = flight.mount(0, 14)
    bonds = tuple(m[0] for m in mounts)
    words, gate, weights = measure.source(14, 1, src, k, bonds)
    h, mark = live.combined(k, src, probe, "stopped")
    np.testing.assert_array_equal(mark, gate)
    a = -1j * t * h
    after = expm_multiply(
        a, np.kron(live.resident(k), np.ones(4) / 2), traceA=a.diagonal().sum()
    ).reshape(-1, 4)
    for setting in flight.ANALYZERS:
        pl, pr = [
            flight.probability(
                flight.sensed_sum(0, 14, m[1:], words, probe, k), t, setting
            )
            for m in mounts
        ]
        reference = [
            weights[(~gate, gate)[i]] @ (left * right)[(~gate, gate)[i]]
            for i in (0, 1)
            for left in (1 - pl, pl)
            for right in (1 - pr, pr)
        ]
        np.testing.assert_allclose(
            live.bins(after, gate, setting), reference, rtol=0, atol=5e-12
        )


def test_live_change_has_observer_counts_and_every_bank():
    raw = archived("live_counts.json")
    d = live.decode(raw)[8]
    assert d["direct_wall_rates_resolving_backaction"] == [4, 8]
    assert d["held_single_probe_rates_resolving_live_from_stopped"] == [17, 32]
    raw["banks"].pop()
    with pytest.raises(ValueError, match="Missing"):
        live.decode(raw)


def test_literal_cut_removes_all_four_complete_contact_terms():
    from . import cut

    k, src, probe = 8, 1, -1
    base, _ = live.combined(k, src, probe, "live")
    missing, _ = live.combined(k, src, probe, "cut")
    delta = base - missing
    mounts, _ = flight.mount(0, 14)
    bonds = tuple(m[0] for m in mounts)
    words, _, _ = measure.source(14, 1, src, k, bonds)
    energies = []
    for a in (1, 1 - probe):
        for b in (1, 1 - probe):
            energy = sum(
                k * (words[:, i] - tip) ** 2
                for j, tip in zip(bonds, (a, b), strict=True)
                for i in (j, (j + 1) % 14)
            )
            energies.append(energy)
    np.testing.assert_allclose(
        delta.diagonal(), np.stack(energies, axis=1).ravel(), atol=2e-14, rtol=0
    )
    assert delta.count_nonzero() == np.count_nonzero(delta.diagonal())
    assert cut.decode(archived("cut_counts.json"))[8]["comparisons"]["live_detached"][
        "wall_rate"
    ] == [4, 8]
