"""Independent labelled hard-core lift, raw contact columns, and complete counter checks."""

import gzip
import json
from itertools import permutations
import numpy as np
import pytest
from scipy.linalg import expm
from . import certificate, exchange as x, measure


def labelled(k):
    words = list(permutations(range(8), 2))
    h = np.zeros((len(words), len(words)))
    raw0, raw1 = x.LAW.local_term_parts()
    raw = raw0 + k * raw1
    for col, w in enumerate(words):
        for a in range(8):
            b = (a + 1) % 8
            word = tuple(int(i not in w) for i in (a, b))
            j = x.LAW.SRC.index(word)
            support = np.flatnonzero(abs(raw[:, j]) > 1e-12)
            assert all(i < 9 and x.LAW.SRC[i] in (word, word[::-1]) for i in support)
            h[col, col] += raw[j, j]
            if word[0] != word[1]:
                old, new = (a, b) if a in w else (b, a)
                hopped = tuple(new if i == old else i for i in w)
                h[words.index(hopped), col] += raw[x.LAW.SRC.index(word[::-1]), j]
    return words, h


@pytest.mark.parametrize("kind,sign", [("native", 1), ("exchange_minus", -1)])
def test_two_ordered_vacancies_give_the_complete_guide_and_both_arms(kind, sign):
    words, h = labelled(8)
    cases, holes, actual = x.generator(2, 8, kind)
    frame = np.zeros((len(words), len(cases)))
    for i, (a, b) in enumerate(cases):
        frame[words.index((a, b)), i] = 1 / np.sqrt(2)
        frame[words.index((b, a)), i] = sign / np.sqrt(2)
    np.testing.assert_allclose(h @ frame, frame @ actual, atol=2e-12)
    # Independently evolve tagged positions and then recover symmetric or
    # antisymmetric joint amplitudes. The native read has no private tags.
    dur = 2
    step = 0.25
    waves = []
    original = frame[:, cases.index(x.source(2))]
    for direction in (1, -1):
        psi = original.astype(complex)
        for aa, bb, pp in x.PAIR:
            a, b, p = [(direction * i) % 8 for i in (aa, bb, pp)]
            for t in range(8):
                s = (t + 0.5) / 8
                load = np.array(
                    [
                        int(p in w) + (1 - s) * int(a in w) + s * int(b in w)
                        for w in words
                    ]
                )
                phase = np.exp(4j * step * load)
                psi = phase * (expm(-1j * step * h / 8) @ (phase * psi))
        waves.append(frame.T @ psi)
    _, st, e, _ = x.trace(2, kind, dur, step)
    for i in range(2):
        np.testing.assert_allclose(abs(waves[i]) ** 2, st[-1][2][i], atol=2e-12)
    for i, sig in enumerate((1, -1)):
        np.testing.assert_allclose(
            abs(waves[0] + sig * waves[1]) ** 2 / 4, e[i], atol=2e-12
        )


def test_all_written_words_and_native_sharp_read():
    cert = certificate.local_certificate()
    assert cert["written_integer_columns_checked"] == 180
    assert x.geometry()["ports"] == list(x.RING)
    assert len({p[-1] for p in x.initial_stocks()}) == 2
    n = certificate.reflection_sectors()
    assert [r["sharp_target_arm_sign"] for r in n] == [1, 1, -1, 1, -1, -1]
    for count in (1, 2, 3, 4, 5, 6, 7):
        # Reordering across the cut changes sign only for an even particle number.
        h = x.generator(count, 8, "exchange_minus")[2]
        expect = x.generator(count, 8, "flux_pi" if count % 2 == 0 else "native")[2]
        np.testing.assert_array_equal(h, expect)


def test_every_outcome_and_winding_calibration():
    deposit = json.loads(gzip.decompress(measure.RECEIPT.read_bytes()))
    for p in deposit["panels"]:
        n, law = p["holes"], p["law"]
        words = deposit["word_orders"][str(n)]
        ix = words.index(list(x.target(n)))
        negative = law == "flux_pi" or (n == 2 and law == "exchange_minus")
        raw = np.asarray(p["fringe"]["counts"])
        assert raw.sum() == measure.SHOTS and raw[1].sum() == p["dark_counts"]
        assert p["fringe"]["probability"][int(not negative)][ix] < 1e-15
        assert raw[int(not negative), ix] == 0
        for sample in p["stages"]:
            for arm in sample["arms"]:
                assert sum(arm["counts"]) == measure.SHOTS
                assert (
                    sum(arm["counts"][words.index(w)] for w in arm["tracked_words"])
                    == arm["localized_counts"]
                )
        assert (
            abs(sum(p["incoherent"]["counts"][0]) / measure.SHOTS - 0.5) < measure.EPS
        )
        if p["step_duration"] == 4096 and law in (
            "native",
            "exchange_minus",
            "flux_pi",
        ):
            assert (
                min(
                    a["localized_counts"] / measure.SHOTS
                    for s in p["stages"]
                    for a in s["arms"]
                )
                > 0.95
            )
            assert abs(p["dark_counts"] / measure.SHOTS - int(negative)) < 0.05
        if law in ("parked", "zero_contrast"):
            assert p["dark_counts"] == 0
    for c in deposit["clock_calibration"]:
        assert c["interval"][0] < c["setting"] < c["interval"][1]
