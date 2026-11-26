from fractions import Fraction

import numpy as np

from model.executable_model import federation
from .finite_spectral_certificate import certify, positive_leading_minors, sharp_pair_bound


def test_exact_gap_and_coupling_certificates_on_all_enrolled_boards():
    for level, cells, seams in ((0, 20, 630), (1, 80, 2520), (2, 320, 10080)):
        witness = certify(federation(level))
        assert (witness["cells"], witness["seams"]) == (cells, seams)
        assert witness["minimum_integer_slack"] > 0
        assert witness["old_rayleigh_excess"] > 0
    assert sharp_pair_bound(Fraction(5701, 2000)) == 1
    assert sharp_pair_bound(3) == 1
    assert sharp_pair_bound(8) == Fraction(17350000, 106069401)
    assert sharp_pair_bound(128) == Fraction(17350000, 62649589401)


def test_integer_inertia_check_rejects_a_singular_or_indefinite_block():
    assert positive_leading_minors(np.array([[2, 1], [1, 2]]))
    assert not positive_leading_minors(np.array([[1, 1], [1, 1]]))
    assert not positive_leading_minors(np.array([[1, 2], [2, 1]]))
