import json
from pathlib import Path

import pytest
from sympy import Matrix, Rational, diag

from . import certificate


def test_coherent_neutral_input_and_daughters():
    q, g = certificate.exterior()
    minus, plus, dark = certificate.signed_effects(q, g)
    # A coherence across neutral labels leaves each signed count unchanged.
    basis = [0, 15, 16, 31]
    v = Matrix([1, 2, -2, 1])
    assert [
        (v.T * op.extract(basis, basis) * v)[0] / v.dot(v) for op in (minus, plus, dark)
    ] == [Rational(2, 5), Rational(2, 5), Rational(1, 5)]
    assert g["T-"][23, 15] == g["T+"][8, 16] == 1


def test_resolved_address_and_wrong_doublet():
    q, g = certificate.exterior()
    # Adding the two electric readings gives zero; subtracting them
    # measures the pair with a squared contrast of four.
    for a, b in [(7, 8), (24, 23), (7, 24)]:
        assert (q[a, a] + q[b, b]) ** 2 == 0
        assert (q[a, a] - q[b, b]) ** 2 == 4
    with pytest.raises(AssertionError):
        certificate.matches(g, [7], [8, 16], [8, 16])
    # A resolved spectator operation cannot be intertwined with an empty site.
    local_pair = diag(q[7, 7] ** 2, q[7, 7] ** 2)
    assert local_pair.det() == 1 and local_pair != diag(0, 0)


def test_receipt():
    receipt = Path(__file__).with_name("results.json")
    assert certificate.run() == json.loads(receipt.read_text())
