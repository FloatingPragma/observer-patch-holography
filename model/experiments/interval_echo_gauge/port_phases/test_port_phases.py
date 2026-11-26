"""Regressions for signed gain, missing enrollments, Zeno branch and recorded charge."""

import json
from math import pi

import numpy as np
import pytest

from .apparatus import LAW, writer
from .certificate import local_certificate
from .code_frame import N, evaluate, membership, unitary_word
from .code_readout import instruments, read
from .measure import HERE, decode, unpack
from .exact_bands import certify


def test_all_ternary_writes_return_guard_with_linear_relative_phase():
    cols = [LAW.SRC.index((n, 1)) for n in range(3)]
    for angle in (0, pi / 11, pi, 2 * pi + pi / 7):
        gate = writer(angle, "native")
        phases = gate[cols[0], cols[0]] * np.exp(-1j * angle * np.arange(3))
        np.testing.assert_allclose(
            gate[:, cols], np.eye(9)[:, cols] * phases, atol=8e-14, rtol=0
        )
    assert local_certificate()["exact_dressed_charge_commutators_zero"]


def test_missing_records_do_not_count_as_measured_failures():
    raw = json.loads((HERE / "records.json").read_text())
    assert decode(raw)["maximum_held_rate_residual"]["native"] < 0.006
    with pytest.raises(ValueError):
        unpack(
            raw["flights"][:-1], iter([tuple(row[:-1]) for row in raw["flights"]]), 3
        )


def test_signed_compiler_preserves_unaddressed_sign_and_pointer_reads_spare():
    target = np.eye(N, dtype=complex)
    target[:2, :2] = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
    actual = evaluate(unitary_word(target, 1))
    np.testing.assert_allclose(actual[0], np.eye(N), atol=3e-13)
    np.testing.assert_allclose(actual[1], target, atol=3e-13)
    for mask in membership()[0]:
        assert list(np.flatnonzero(mask.diagonal())) == [1, 3, 5]


def test_all_held_code_frequencies_and_controls():
    _, ideal, fine = instruments()
    raw = json.loads((HERE / "code_records.json").read_text())
    summary = read(raw, ideal, fine)
    errs = summary["maximum_raw_frequency_errors"]
    assert (
        errs["declared"] < 0.003
        and errs["fine"] > 0.03
        and errs["reversed_order"] > 0.015
    )
    raw[0][-1][-1] += 1
    with pytest.raises(ValueError):
        read(raw, ideal, fine)


def test_exact_spectral_projectors_and_current_rank():
    result = certify()
    assert result["compressed_current_gram_determinant_mod_prime"] == 1
    assert [row["leading_rank"] for row in result["exact_leading_bands"]] == [3, 2]
