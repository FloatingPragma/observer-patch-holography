import json

import numpy as np
import pytest

from ...full_geometry.measure import leg_kernels, pairing, ranging
from ..measure import blank_clock, contacts, density_center
from . import measure as m


def test_first_arrivals_and_native_clock_gaps():
    saved = json.loads((m.HERE / "results.json").read_text())
    board = pairing.formation.federation(1)
    marks = pairing.formation.region(board)
    states = density_center()[0]
    for charge in (-1, 1):
        weights = np.zeros((2, 11))
        for parity, block in enumerate(states):
            h = ranging.generator(board, marks, charge, 8, bool(parity))
            panels = leg_kernels(h, (marks[1], marks[2]))
            for sender in (0, 1):
                for n, a in enumerate(ranging.first_hits(panels, sender, 1 - sender)[1:12]):
                    weights[sender, n] += np.trace(a @ block @ a.conj().T).real
        np.testing.assert_allclose(weights, saved["first_hit"][str(charge)], rtol=0, atol=2e-13)
        gaps = np.linalg.eigvalsh(blank_clock(charge, 8))
        np.testing.assert_allclose(gaps, [8, 9], atol=1e-13, rtol=0)
    np.testing.assert_allclose(np.linalg.eigvalsh(contacts(0, 8)), [8, 8, 9, 9], atol=1e-13, rtol=0)


def test_counts_against_local_readers_and_bilateral_inference():
    record = json.loads((m.HERE / "results.json").read_text())
    assert len(record["panels"]) == 8
    for row in record["panels"]:
        beta = m.infer(row)["radial_beta"]
        target = 0.8 if row["late_wait_multiple"] else 0.0
        assert beta[0] < target < beta[1]
        assert max(abs(x - target) for x in beta) < (0.00008 if target else 0.00054)
        for direction, arm in enumerate(row["directions"]):
            p = np.array(record["first_hit"][str(row["charge"])][direction][: row["ticks"]])
            tick = 0.5 * np.arange(1, row["ticks"] + 1)
            for i, t in enumerate((0, row["period"])):
                for tag, z in (("local", np.exp(1j * t)),
                               ("remote", p @ np.exp(1j * (tick + (1 + row["late_wait_multiple"]) * t)))):
                    if tag == "remote":
                        for bank in arm[tag][i]:
                            observed_bins = np.sum(bank["by_hit"], axis=1) / m.SHOTS
                            np.testing.assert_allclose(observed_bins, p, rtol=0, atol=m.SCORE_ERROR)
                    observed, _ = m.phase_and_error(arm[tag][i])
                    assert max(abs((z - observed).real), abs((z - observed).imag)) < m.SCORE_ERROR


def test_no_arrivals_and_incomplete_bank():
    with pytest.raises(ValueError, match="no resolved"):
        m.phase_and_error([dict(by_hit=[[0, 0]], missed=m.SHOTS)] * 2)
    with pytest.raises(ValueError, match="incomplete"):
        m.phase_and_error([dict(by_hit=[[0, 1]], missed=m.SHOTS)] * 2)
