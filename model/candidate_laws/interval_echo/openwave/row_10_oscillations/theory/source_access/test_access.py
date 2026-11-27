"""Contact support, recording of the spectator and complete-trial discrimination."""

import json
from pathlib import Path

import numpy as np
import pytest

from model.candidate_laws.interval_echo.openwave.row_10_oscillations.theory.source_access.access import (
    calculations,
    partitions,
    record,
    sharp_read,
    tv,
)


@pytest.fixture(scope="module")
def replay():
    return calculations()


def test_each_detour_opens_its_own_word_difference():
    for strength in (4, 8):
        for setting in range(4):
            graph = partitions(setting, strength)
            assert graph["shortest_from_01_to_10"] == (10 if setting == 3 else None)
            assert len(graph["hub_components"]) == 2 ** (2 - setting.bit_count())
    for bit in (0, 1):
        np.testing.assert_array_equal(sum(sharp_read(bit)), np.eye(4))


def test_single_loop_controls_and_all_archived_probabilities(replay):
    archive = json.loads(Path(__file__).with_name("probabilities.json").read_text())
    assert len(replay["panels"]) == len(archive["panels"]) == 6
    for row, expected in zip(replay["panels"], archive["panels"], strict=True):
        assert (row["owner"], row["time"]) == (expected["owner"], expected["time"])
        np.testing.assert_allclose(
            row["original_arrival"], expected["original_arrival"], rtol=3e-7, atol=0
        )
        np.testing.assert_allclose(
            row["conditional_tv_ceiling"], expected["conditional_tv_ceiling"], atol=1e-7
        )
        assert len(row["records"]) == 22 and row["conditional_tv_ceiling"][0] < 1e-9
        for bank, held in zip(row["records"], expected["records"], strict=True):
            for key in ("arm", "panel", "analyzer_strength"):
                assert bank[key] == held[key]
            for key in ("plus", "minus"):
                np.testing.assert_allclose(bank[key], held[key], atol=1e-7)
            if bank["panel"] in ("park", "cut", "dephased"):
                assert bank["tv_joint"] < 1e-9
            if bank["panel"] == "open" and bank["analyzer_strength"] == 8:
                assert bank["tv_vacancy"] > 9e-5


def test_two_rotated_native_reads_survive_all_stable_word_phases(replay):
    # The exact bound is proved algebraically. These phases exercise two independent
    # spectator sectors and detect wrong signs or a misplaced phase gate in the pulses.
    for row in replay["panels"]:
        delta = np.zeros((4, 4), complex)
        for i, j, x, y in row["conditional_source_coherence_difference"]:
            delta[i, j] = complex(x, y)
            delta[j, i] = delta[i, j].conjugate()
        for arm in (0, 1):
            unit = row["conditional_tv_ceiling"][1 << arm]
            for rotate in (0, 0.31, 1.45, 2.38):
                phases = np.exp(1j * np.array([0.8, -2.3, 0.3, 1.7]) * rotate)
                d = delta * phases[:, None] * phases.conj()
                a, b = np.eye(4) / 4 + d / 4, np.eye(4) / 4 - d / 4
                # Test states differ by d/2. Sum, rather than mean, compensates
                # for that factor when comparing with the original trace distance.
                two_reads = sum(
                    tv(record(a, arm, 8, p), record(b, arm, 8, p))
                    for p in ("open", "rotated")
                )
                assert two_reads + 1e-10 >= np.sqrt(3) * unit / 4
                assert two_reads <= unit + 1e-10
