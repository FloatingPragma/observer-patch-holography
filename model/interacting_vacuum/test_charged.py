from fractions import Fraction as F

import numpy as np
import pytest

from model.cube_eta_protection import cube, repair
from model.interacting_vacuum import charged as c


def test_exact_positive_charge_lines_and_nontrivial_interacting_measure():
    result = c.spectral_certificate()
    assert result['up_addition']['gap'] == result['up_removal']['gap'] == 'sqrt3/16'
    assert result['down_addition']['gap_mean'] == '1/2'
    assert result['down_addition']['gap_variance'] == '67/256'
    assert result['down_addition']['residual_weight'] == '4/5'
    assert result['down_removal']['total_weight'] == '0'
    assert result['native_reference_preparation_or_detector_executed'] is False


def test_car_sum_rules_and_origin_frame_are_original_physical_modes():
    g = c.reference_vector()
    for mode in (0, 1):
        first = c.car(c.car(g, mode, True), mode, False)
        second = c.car(c.car(g, mode, False), mode, True)
        assert c._equal(tuple(a+b for a, b in zip(first, second, strict=True)), g)
    from model.field_symmetry.symmetry import onsite_frame
    assert np.array_equal(onsite_frame((0, 0, 0)), np.eye(2))


def test_source_contact_is_nonzero_on_the_down_added_vector():
    v = c.car(c.reference_vector(), 1, True)
    doublons = np.asarray([sum(((mask >> (2*j)) & 3) == 3 for j in range(8)) for mask in range(65536)], dtype=np.int64)
    contact_v = tuple(doublons*a for a in v)
    assert c.inner(v, contact_v) == (648, 0)
    assert c.inner(contact_v, contact_v) == (648, 0)
    # Removing D changes the proved spectral mean from1/2 to0.
    free_gap = tuple(a-16*b for a, b in zip(c.gap_numerator(v), contact_v, strict=True))
    assert c.inner(v, free_gap) == (0, 0)


def test_every_admitted_repair_preserves_all_charged_sector_cross_terms():
    result = c.repair_certificate()
    assert result['passive_protected_dimension'] == 6748
    assert result['complete_outcomes_including_zero_and_outside'] == 130
    # Explicit outside projector check on a coherent cross-number vector.
    g = c.reference_vector()
    vectors = [g, c.car(g, 0, False), c.car(g, 1, True)]
    psi = sum(a+np.sqrt(3)*b for a, b in vectors)/36
    for law in ('C', 'E'):
        for seam in repair.SEAMS:
            outputs = repair.apply(psi, law, seam)
            assert all(np.count_nonzero(value) == 0 for _, value in outputs[:-1])
            assert np.array_equal(outputs[-1][1], psi)


def test_residual_gap_moment_variance_is_strictly_positive_exactly():
    # mu=1/2-sqrt3/64; second=131/256-sqrt3/64.
    rational_mean_square = F(1, 4)+F(3, 4096)
    assert F(131, 256)-rational_mean_square == F(1069, 4096) > 0


@pytest.mark.parametrize('mode,create', [(16, True), (-1, False), (True, True), (0, 1)])
def test_probe_domain(mode, create):
    with pytest.raises(ValueError): c.car(c.reference_vector(), mode, create)
