"""Test joint residual/recoil spectra, physical controls and conditional work."""

from hashlib import sha256
import json
from itertools import combinations
import numpy as np
from . import strand_recoil as recoil
from .measure import HERE
from scipy.linalg import expm
from scipy import sparse


def test_full_state_flight_against_dense_exponential_on_multiple_clocks():
    rng = np.random.default_rng(4590)
    raw = rng.normal(size=(24, 24))
    h = raw + raw.T
    stock = rng.normal(size=(24, 4))
    times = [0, 1.1, 7, 12]
    observed = recoil.open_flights(sparse.csr_matrix(h), stock, times)
    for t, p in zip(times, observed, strict=True):
        np.testing.assert_allclose(p, expm(-1j * t * h) @ stock, atol=4e-12)


def test_residual_detector_on_two_entangled_ports_and_both_residents():
    rng = np.random.default_rng(89364)
    terminal = [41, 42]
    cell = [51, 52, 53, 54]
    outside = [55]
    cw = list(combinations(cell, 2))
    words = list(combinations((*terminal, *cell, *outside), 3))
    rows = recoil.embedded_rows(words, terminal, cw)
    sf, _ = np.linalg.qr(rng.normal(size=(4, 4)))
    cf, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    wave = rng.normal(size=(2 * len(words), 2)) + 1j * rng.normal(
        size=(2 * len(words), 2)
    )
    wave, _ = np.linalg.qr(wave)
    sgroups = [np.array([0]), np.arange(1, 4)]
    cgroups = [np.arange(3), np.arange(3, 6)]
    results = recoil.final_joint(wave, rows, sf, cf, sgroups, cgroups)
    rows_by_word = {frozenset(w): i for i, w in enumerate(words)}
    for i, sg in enumerate(sgroups):
        for j, cg in enumerate(cgroups):
            expected = np.zeros((2, 2), complex)
            for s in sg:
                for c in cg:
                    detector = np.zeros(len(wave))
                    for a in (0, 1):
                        for e in (0, 1):
                            for p, w in enumerate(cw):
                                k = (
                                    a * len(words)
                                    + rows_by_word[frozenset((terminal[1 - e], *w))]
                                )
                                detector[k] = sf[2 * a + e, s] * cf[p, c]
                    amp = detector @ wave
                    expected += np.outer(amp.conj(), amp)
            np.testing.assert_allclose(expected, results[i, j], atol=3e-13)


def test_complete_two_resident_neutral_orbit():
    b, d, h, _ = recoil.ws.ring_layout(8, True, False)
    _, cell, edges = recoil.exterior.geometry(b, d, h)
    assert recoil.neutral_orbit(cell, edges, list(combinations(cell, 2))) == 66
    calibration = recoil.calibrate_receiver()
    known = calibration["site_phase"][str(calibration["entrance"])]
    assert 0.84 < known["true"] < 0.85
    assert abs(known["true"] - known["read"]) < 3e-6
    assert calibration["cold_false"] == 0
    assert calibration["true_clock_min"] > 0.999999
    assert calibration["three_read_false_max"] < 3.5e-10
    assert calibration["three_read_true_min"] > 0.999998


def test_closed_exit_factorizes_and_source_probabilities_use_all_originals():
    data = json.loads((HERE / "strand_recoil.json").read_text())
    assert (
        data["protocol_sha256"]
        == sha256((HERE / recoil.DESIGN).read_bytes()).hexdigest()
    )
    assert data["tested"]["6"]["standards"]["first_higher_cell"] > 6
    rec = data["tested"]["8"]
    for rows in rec["trials"].values():
        for row in rows:
            cond = np.array(row["conditional_by_incoming_level"])
            assert np.max(cond.sum(axis=(0, 1))) <= 1 + 2e-10
            assert (
                abs(row["open_work"][2] + row["close_work"][2] - row["net_work"][2])
                < 2e-9
            )
            for arm, s in row["source"].items():
                p = np.array(rec["original_mode_probabilities"][arm])
                np.testing.assert_allclose(
                    s["incoming_second_residual_first_cell_doublet"],
                    p[1] * cond[1, 1, 2],
                    atol=1e-13,
                )
                np.testing.assert_allclose(
                    p.sum(), np.sum(s["keep"]) + s["inventory_refused_keep"], atol=1e-11
                )
    for row in rec["trials"]["link_cut"]:
        raw = np.array(row["conditional_by_incoming_level"])
        assert raw[1, 1, 2] < 1e-17
    q = rec["bare_board_separate_standard"]
    assert abs(q[-1]["mismatch"]) < 0.00009595 and q[-1]["residual"] < 1e-11
    assert abs(q[0]["mismatch"]) > 0.001


def test_resonance_energy_cost_and_spectral_refusal():
    data = json.loads((HERE / "strand_recoil.json").read_text())["tested"]["8"]
    ref = data["standards"]
    strand = ref["strand_low"]
    cell = [s[0] for s in ref["cell_low"]]
    for line in strand[1:]:
        mismatch = min(
            abs(line - residual - remote) for residual in strand for remote in cell[1:]
        )
        assert mismatch > 0.0038
    assert 0.0038 < strand[2] - strand[1] - cell[1] < 0.0039
    assert 0.09 < data["far_port_transition"] < 0.10
