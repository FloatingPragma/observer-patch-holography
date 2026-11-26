import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo.openwave.row_03_clock.polar_island import (
    measure,
)


def test_block_filter_against_every_delay_and_two_local_projectors(monkeypatch):
    # Independent dense word Hamiltonian, address-resolved read and explicit
    # enumeration of every pair of delay draws on a small ring.
    n, g, length, rest = 4, 3, 5, 15.4
    dim = 2**n
    words = np.arange(dim)
    bits = (words[:, None] >> np.arange(n)) % 2
    walls = (bits != np.roll(bits, 1, axis=1)).sum(axis=1)
    h = np.diag(g * (2 * walls + bits.sum(axis=1))).astype(float)
    for j in range(n):
        h[words, words ^ (1 << j)] = -0.5
    energies, modes = eigh(h)
    starts = measure.source_words(n)
    expected = np.zeros((2, 3, 4))
    for r in range(length):
        for s in range(length):
            t = abs(r - s) / 128
            evolved = modes @ (np.exp(-1j * t * energies)[:, None] * (modes.T @ starts))
            for j in range(n):
                upper = bits[:, j] == 1
                other = evolved[words ^ (1 << ((j + 1) % n))]
                table = []
                for positive in [upper, ~upper]:
                    for sign in [1, -1]:
                        table.append(
                            np.sum(abs((evolved + sign * other)[positive]) ** 2, axis=0)
                            / 4
                        )
                p = np.stack(table, axis=1) / (n * length**2)
                expected[0] += p
                expected[1] += np.cos(rest * t) * p
    monkeypatch.setattr(measure, "WINDOWS", (length,))
    all_reads, *_ = measure.block_measures(n, g, rest)
    for s in range(3):
        ordinary = measure.joint(1, all_reads[0, 0, 0, s])
        signed = measure.joint(measure.kernel(rest, length), all_reads[0, 0, 2, s])
        np.testing.assert_allclose(
            [ordinary, signed], expected[:, s], atol=3e-13, rtol=0
        )


def test_all_banks_and_grade():
    archived = json.loads(Path(__file__).with_name("records.json").read_text())
    assert measure.run() == archived
    for panel in archived["panels"]:
        assert panel["band_width"] < 0.2
        assert panel["rest_hand"] > 19
        assert sum(b["length"] == measure.WINDOWS[-1] for b in panel["banks"]) == 27
        assert panel["n"] * 3 < 10 * panel["kappa"]
        for b in panel["banks"]:
            counts = np.array(b["counts"])
            assert counts.shape == (4, 2) and np.all(counts >= 0)
            assert counts.sum() == measure.SHOTS
            assert measure.tagged_contrast(counts) == b["tagged_contrast"]

        def triplet(source, control):
            return [
                b["tagged_contrast"]
                for b in panel["banks"]
                if b["length"] == measure.WINDOWS[-1]
                and b["source"] == source
                and b["control"] == control
            ]

        sham = triplet("all_lower", "free")
        dephased = triplet("all_lower", "ground_dephased")
        assert abs(sham[1]) - max(abs(sham[0]), abs(sham[2])) > 8 / np.sqrt(
            measure.SHOTS
        )
        assert max(map(abs, dephased)) < 2 / np.sqrt(measure.SHOTS)
        for source in measure.SOURCES[:2]:
            a, center, b = map(abs, triplet(source, "ground_dephased"))
            assert center - max(a, b) < 8 / np.sqrt(measure.SHOTS)
            assert max(map(abs, triplet(source, "selected_island"))) < 4 / np.sqrt(
                measure.SHOTS
            )


def test_all_excited_transition_enclosures():
    from model.candidate_laws.interval_echo.openwave.row_03_clock.polar_island import (
        exclusion,
    )

    archived = json.loads(Path(__file__).with_name("enclosures.json").read_text())
    for result, saved in zip(exclusion.run(), archived, strict=True):
        assert result["rest_order_index"] == saved["rest_order_index"] == 1
        assert (result["n"], result["kappa"]) == (saved["n"], saved["kappa"])
        assert 0 < result["ordered_eigenvalue_error"] < 9e-9
        assert (
            1e-6 < result["strict_detuning_lower"] < result["min_excited_gap_detuning"]
        )
        np.testing.assert_allclose(
            result["min_excited_gap_detuning"],
            saved["min_excited_gap_detuning"],
            atol=1e-11,
            rtol=0,
        )
    # Direct dense eigensolve and exhaustive pair enumeration on a different
    # finite ring check the sorted insertion-neighbor search independently.
    n, contrast = 6, 4
    certificate = exclusion.certify(n, contrast)
    dimension = 2**n
    h = np.zeros((dimension, dimension))
    for word in range(dimension):
        digit = [(word >> i) % 2 for i in range(n)]
        h[word, word] = contrast * (
            sum(digit) + 2 * sum(digit[i] != digit[i - 1] for i in range(n))
        )
        for i in range(n):
            h[word, word ^ (1 << i)] = -0.5
    levels = np.linalg.eigvalsh(h)
    expected = levels[1] - levels[0]
    nearest = np.min(abs(levels[1:, None] - levels[None, 1:] - expected))
    assert (
        abs(expected - certificate["ground_to_rest"])
        <= 2 * certificate["ordered_eigenvalue_error"]
    )
    assert nearest >= certificate["strict_detuning_lower"]
    assert abs(nearest - certificate["min_excited_gap_detuning"]) < 1e-11
