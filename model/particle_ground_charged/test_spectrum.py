"""Full ground binding, exact frame admission, all channels and offset controls."""
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
import sympy as s

from model.particle_ground_charged import spectrum as q


def test_full_source_reference_support_and_complete_car_probe_weights():
    state = q.ground_input()
    assert (len(state['u']), len(state['Au'])) == (53, 160)
    assert {mask.bit_count() for scalar, mask in state['u']} == {6}
    assert {mask.bit_count() for scalar, mask in state['Au']} == {8}
    for addition, empty_number, occupied_number in ((False, 6, 8), (True, 12, 10)):
        first = sum(sum(v*v for v in q.probe(state['u'], i, addition).values()) for i in range(18))
        second = sum(sum(v*v for v in q.probe(state['Au'], i, addition).values()) for i in range(18))
        assert first == empty_number*state['norm']
        assert second == occupied_number*state['image_norm']


def test_every_representative_block_uses_original_source_and_retained_complement():
    dimensions = {(0, 'left'): (165, 80), (0, 'singlet'): (53, 240),
                  (1, 'left'): (60, 212), (1, 'singlet'): (144, 93),
                  (2, 'left'): (60, 212), (2, 'singlet'): (144, 93)}
    for (family, kind), expected in dimensions.items():
        for addition, dimension in zip((False, True), expected, strict=True):
            mode = 6*family+(0 if kind == 'left' else 4)
            block = q.source_block(mode, addition)
            assert block['H'].shape == (dimension, dimension)
            assert dimension+block['orthogonal_complement_dimension'] == 786432
            assert (block['H']-block['H'].T).nnz == 0
            assert {mask.bit_count()-2*(scalar != 0) for scalar, mask in block['basis']} == {7 if addition else 5}
            indices = {state: i for i, state in enumerate(block['basis'])}
            for column, state in enumerate(block['basis']):
                expected_column = dict(q.literal_vertex(*state))
                if state[0]:
                    expected_column[state] = expected_column.get(state, 0)+1
                actual = block['H'][:, column].tocoo()
                assert {block['basis'][int(row)]: int(value) for row, value in zip(actual.row, actual.data)} == expected_column


def test_independent_small_matrix_energy_and_full_projector_enclosures():
    h = sparse.csr_matrix([[0, 2], [2, 1]], dtype=np.int64)
    frame = q.admit_frame(h, *q.exact_proposal(h))
    e = [(1-s.sqrt(17))/2, (1+s.sqrt(17))/2]
    H = s.Matrix([[0, 2], [2, 1]])
    assert len(frame['groups']) == 2
    for index, group in enumerate(frame['groups']):
        low, high = group['lower'], group['upper']
        assert s.Rational(low) < e[index] < s.Rational(high)
        exact = (H-e[1-index]*s.eye(2))/(e[index]-e[1-index])
        v = s.Matrix(frame['q'][:, group['indices']].tolist())/frame['scale']
        error = np.linalg.norm(np.array((exact-v*v.T).evalf(30)).astype(float), 2)
        assert error <= float(group['projector_operator_error'])
    assert sum(len(group['indices']) for group in frame['groups']) == 2


def test_frame_admission_rejects_missing_columns_floats_bad_residual_and_overflow():
    h = sparse.csr_matrix([[0, 2], [2, 1]], dtype=np.int64)
    v, e, scale = q.exact_proposal(h)
    with pytest.raises(ValueError):
        q.admit_frame(h, v[:, :1], e, scale)
    with pytest.raises(ValueError):
        q.admit_frame(h.astype(float), v, e, scale)
    with pytest.raises(ArithmeticError):
        q.admit_frame(h, np.zeros((2, 2), dtype=np.int64), e, scale)
    with pytest.raises(ArithmeticError):
        q.admit_frame(h, v, e+scale, scale)
    with pytest.raises(OverflowError):
        q.admit_frame(h, np.full((2, 2), 1 << 61, dtype=np.int64), e, scale)
    huge = sparse.csr_matrix(np.diag(np.array([np.iinfo(np.int64).min, 1], dtype=np.int64)))
    with pytest.raises(OverflowError):
        q.admit_frame(huge, v, e, scale)
    wrapped_asymmetry = sparse.csr_matrix([[0, np.iinfo(np.int64).min], [0, 0]], dtype=np.int64)
    with pytest.raises(ValueError):
        q.admit_frame(wrapped_asymmetry, v, e, scale)


def _source_action(vector):
    result = {}
    for state, value in vector.items():
        if state[0]:
            result[state] = result.get(state, 0)+value
        for destination, coefficient in q.literal_vertex(*state):
            result[destination] = result.get(destination, 0)+int(coefficient)*value
    return {state: value for state, value in result.items() if value}


def test_full_multiplet_moment_equalities_use_every_original_car_mode():
    source = q.ground_input()
    for family in range(3):
        for kind, modes in (('left', range(6*family, 6*family+4)), ('singlet', range(6*family+4, 6*family+6))):
            for addition in (False, True):
                expected = None
                for mode in modes:
                    x = q.probe(source['u'], mode, addition)
                    y = q.probe(source['Au'], mode, addition)
                    hx, hy = dict(x), dict(y); moments = []
                    for power in range(3):
                        moments.append((sum(value*hx.get(state, 0) for state, value in x.items()),
                                        sum(value*hy.get(state, 0) for state, value in y.items()),
                                        sum(value*hy.get(state, 0) for state, value in x.items())))
                        hx, hy = _source_action(hx), _source_action(hy)
                    if expected is None:
                        expected = moments
                    else:
                        assert moments == expected


def test_all_certified_spectral_groups_retain_rank_and_original_weight():
    totals_lo = totals_hi = F(0)
    for family in range(3):
        for kind in ('left', 'singlet'):
            for addition in (False, True):
                channel = q.certified_channel(family, kind, addition)
                groups = channel['all_spectral_groups']
                assert sum(row['rank_in_reducing_source_block'] for row in groups) == channel['actual_probe_reducing_block_dimension']
                assert all(row['excitation_interval'][0] > 0 for row in groups)
                assert any(row['unscaled_single_mode_weight'][0] == 0 for row in groups)
                assert all(0 <= row['unscaled_single_mode_weight'][0] <= row['unscaled_single_mode_weight'][1] <= 1 for row in groups)
                lo, hi = channel['complete_original_probe_weight_interval']
                assert sum(row['original_complete_instrument_weight'][0] for row in groups) <= hi
                assert sum(row['original_complete_instrument_weight'][1] for row in groups) >= lo
                totals_lo += lo; totals_hi += hi
                assert not channel['dark_or_zero_groups_discarded']
    assert totals_lo <= 1 <= totals_hi


def test_nonzero_dominant_charged_lines_are_actual_ground_readouts():
    removal = q.certified_channel(0, 'left', False)['all_spectral_groups'][0]
    addition = q.certified_channel(0, 'left', True)['all_spectral_groups'][0]
    assert F(70, 1000) < removal['excitation_interval'][0] < removal['excitation_interval'][1] < F(72, 1000)
    assert removal['unscaled_single_mode_weight'][0] > F(66, 100)
    assert F(48, 100) < addition['excitation_interval'][0] < addition['excitation_interval'][1] < F(49, 100)
    assert addition['unscaled_single_mode_weight'][0] > F(24, 100)


def test_same_reference_offset_cannot_symmetrize_full_or_conditional_laws():
    report = q.offset_obstruction()
    assert len(report['channel_controls']) == 6
    assert all(row['equal_full_raw_laws_excluded'] for row in report['channel_controls'])
    assert all(row['equal_translated_conditional_shapes_excluded'] for row in report['channel_controls'])
    first = report['channel_controls'][0]['offset_required_to_equalize_conditional_means']
    second = report['channel_controls'][2]['offset_required_to_equalize_conditional_means']
    assert first.hi < 0 < second.lo
    assert not report['physical_CPT_violation_claimed']
    assert report['outside_that_window_ground_must_be_recomputed']


def test_literal_moment_zero_matches_original_weight_without_conditioning():
    for mode in (0, 4, 6, 10, 12, 16):
        for addition in (False, True):
            block = q.source_block(mode, addition)
            moments = q.matrix_moments(block)
            assert moments[0]['trial_interval'] == q.original_probe_weight(block)
            assert len(moments) == 5
            assert all(row['true_ground_interval'].lo <= row['trial_interval'].lo <= row['trial_interval'].hi <= row['true_ground_interval'].hi
                       for row in moments)


def test_primary_source_paths_and_complete_readout_scope():
    root = Path(__file__).resolve().parents[2]
    assert all((root/path).is_file() for path in q.SOURCE_INPUTS)
    data = q.certificate()
    assert data['physical_particle_assignments'] == data['physical_mass_predictions'] == []
    assert data['measured_particle_inputs'] == []
    assert not data['floating_eigenvalues_accepted_as_proof']
    assert not data['conditional_branch_histograms_substituted']
    assert data['probe_reference_compensates_gauge_representation']
