from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from . import certificate as c
from model.openwave_coupled_haar_source.source import hamiltonian


def test_all_source_spectra_records_and_independent_control():
    result = c.run()
    assert len(result['records']) == 6
    for row in result['records']:
        assert F(row['full_Haar_gap_interval_exact'][0]) > 0
        assert F(row['positive_full_complement_floor_exact']) > 0
        assert row['scope']['cap_is_only_a_proof_coordinate']
        assert len(row['complete_joint_observer_instrument']['four_outcomes']) == 4
        covariance = tuple(map(F, row['connected_Wilson_covariance_interval_exact']))
        if row['shared_edge']:
            assert covariance[0] > 0
            assert row['first_three_levels_identified_by_disjoint_windows_and_fourth_floor']
        else:
            assert covariance[0] <= 0 <= covariance[1]
    primary = result['records'][2]
    gap = tuple(map(F, primary['full_Haar_gap_interval_exact']))
    assert gap[1]-gap[0] < F(1, 10**9)
    # A larger proof coordinate gives an independent numerical sanity check;
    # authority remains the exact residual and full-complement certificate.
    eigenvalues = np.linalg.eigvalsh(np.asarray(hamiltonian(9), dtype=float))
    assert float(gap[0]) <= eigenvalues[1]-eigenvalues[0] <= float(gap[1])


def test_source_and_complete_frame_mutations_are_rejected(monkeypatch):
    original = deepcopy(c.witnesses())
    key = 'shared_g2_2_cap4'
    for mutation in ('source', 'labels', 'frame', 'trial'):
        damaged = deepcopy(original)
        if mutation == 'source':
            damaged[key]['exact_source_matrix_sha256'] = '0'*64
        elif mutation == 'labels':
            damaged[key]['ordered_spin_network_labels'].reverse()
        elif mutation == 'frame':
            damaged[key]['complete_finite_frame']['dyadic_eigenframe'][0][0] = '0'
        else:
            damaged[key]['low_trial_integer_coordinates'][1] = ['0']*55
        monkeypatch.setattr(c, 'witnesses', lambda: damaged)
        c.calculate.cache_clear()
        with pytest.raises((AssertionError, ValueError)):
            c.calculate(key)
    monkeypatch.setattr(c, 'witnesses', lambda: original)
    c.calculate.cache_clear()


def test_algebraic_rounding_encloses_actual_operator():
    matrix = hamiltonian(2)
    rounded, error = c.rounded_operator(matrix, denominator=2**12)
    discrepancy = np.asarray(matrix, float)-np.asarray(rounded, float)/2**12
    assert max(abs(np.linalg.eigvalsh(discrepancy))) <= float(error)
