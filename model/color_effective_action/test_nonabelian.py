"""Conditional nonabelian regulator and actual source-identity controls."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s
from scipy.sparse.linalg import expm_multiply

from model.interacting_extension.higgs_join.interacting_carrier import NativeCarrier, build_join
from . import nonabelian


def symmetric_lift(generator, degree):
    d = generator.rows
    basis = [x for x in product(range(degree+1), repeat=d) if sum(x) == degree]
    index = {x:i for i, x in enumerate(basis)}
    answer = s.zeros(len(basis))
    for j, occupation in enumerate(basis):
        for a in range(d):
            for b in range(d):
                if not occupation[b] or not generator[a, b]:
                    continue
                target = list(occupation)
                target[b] -= 1
                target[a] += 1
                answer[index[tuple(target)], j] += occupation[b]*generator[a, b]
    return answer


def test_native_normalization_and_missing_conjugates():
    row = nonabelian.growth_certificate()
    assert row['minimum_dimension_after_all_conjugates_added'] == 119
    assert row['native_link_dimension'] == 60
    assert row['zero_electric_energy_vacuum_coordinate_unitarity_defect'] == '1'
    assert row['source_normalized_Fierz_identities'][0]['full_tensor_matrix_entries_checked'] == 81
    assert row['conditional_same_control_rotor_bounds']['total_bounded_perturbation'] == '301'
    assert not row['actual_native60_source_selects_this_completion']


@pytest.mark.parametrize('degree', (1, 2, 3, 4))
def test_source_symmetric_tensor_Casimirs_without_dimension_guessing(degree):
    carrier = NativeCarrier()
    color_indices = carrier.data['blocks'][0]['indices']
    weak_indices = carrier.data['blocks'][1]['indices']
    color = tuple(s.I*g.extract(color_indices, color_indices)/2 for g in carrier.anti[:7])
    color += (s.I*carrier.anti[7].extract(color_indices, color_indices)/(2*s.sqrt(3)),)
    weak = tuple(s.I*g.extract(weak_indices, weak_indices)/2 for g in carrier.anti[8:11])
    dims = []
    for generators, expected in ((color, s.Rational(degree*(degree+3), 3)),
                                (weak, s.Rational(degree*(degree+2), 4))):
        lifted = [symmetric_lift(g, degree) for g in generators]
        n = lifted[0].rows
        casimir = sum((g*g for g in lifted), s.zeros(n)).applyfunc(s.simplify)
        assert casimir == expected*s.eye(n)
        dims.append(n)
    row = nonabelian.highest_Q_power(degree)
    assert row['representation_dimension'] == dims[0]*dims[1]
    assert F(row['native_Casimir']) == 2*F(degree*(degree+3), 3)+F(degree*(degree+2), 8)+F(degree*degree, 12)


def test_full_cutoff_coordinate_identity_is_a_low_energy_core_property():
    cutoff = 5
    labels = list(range(-cutoff, cutoff+1))
    shift = np.zeros((len(labels),)*2)
    for i, charge in enumerate(labels):
        if charge+1 in labels:
            shift[labels.index(charge+1), i] = 1
    core = [i for i, charge in enumerate(labels) if abs(charge) <= cutoff-1]
    assert np.array_equal((shift.T@shift)[:, core], np.eye(len(labels))[:, core])
    assert np.array_equal((shift@shift.T)[:, core], np.eye(len(labels))[:, core])
    assert np.linalg.norm(shift.T@shift-np.eye(len(labels)), 2) == 1
    rng = np.random.default_rng(99)
    vector = rng.normal(size=len(labels))+1j*rng.normal(size=len(labels))
    vector /= np.linalg.norm(vector)
    energy = float(np.sum(np.array(labels)**2*abs(vector)**2))
    defect = 1-np.linalg.norm(shift@vector)**2
    assert defect <= energy/(cutoff-1)**2


def test_conditional_nonabelian_path_bounds_retain_growing_norms():
    for radius in (8, 16, 32, 64):
        row = nonabelian.cutoff_bounds(F(radius))
        count = radius//2+1
        assert row['minimum_exit_insertions'] == count
        assert F(row['conditional_full_to_cutoff_resolvent_bound']) == F(2, 301)*F(1, 2)**count
        assert F(row['conditional_local_action_resolvent_bound']) == F(2, 301)*F(1, 2)**(2*count)
        assert row['complete_coordinate_family_premise_required']
        assert not row['source_capacity_selection_inferred']


def test_completed_source_experiment_is_distinguishable_at_finite_time():
    receipt = nonabelian.completion_response_witness()
    native, coordinate = receipt['finite_time_experiments']
    a, b = map(F, native['complete_flux_over_Higgs_probability_ratio'])
    c, d = map(F, coordinate['complete_flux_over_Higgs_probability_ratio'])
    assert a < F(9, 2) < b < c < 9 < d
    source = build_join()['sectors'][0]['H']
    from model.color_loop_readout.source import Sector
    sector = Sector(0)
    seed = np.zeros(438)
    seed[sector.seed] = 1
    evolved = expm_multiply(-1j*float(F(receipt['same_source_time']))*source, seed)
    flux = [i for i, row in enumerate(sector.labels) if row['link0_block'] != -1 or row['link1_block'] != -1]
    scalar = [i for i, row in enumerate(sector.labels) if row['N_H'] == 1]
    ratio = np.sum(abs(evolved[flux])**2)/np.sum(abs(evolved[scalar])**2)
    assert float(a) < ratio < float(b)
    assert not receipt['identical_complete_source_records']
    assert not receipt['different_low_energy_universality_class_proved']


def test_invalid_continuum_premises_cannot_generate_a_bound():
    with pytest.raises(ValueError, match='margins'):
        nonabelian.cutoff_bounds(F(8), mu=F(301))
    with pytest.raises(ValueError, match='margin'):
        nonabelian.coordinate_core_error(2, 2, 1)
    with pytest.raises(TypeError, match='exact'):
        nonabelian.coordinate_core_error(8., 2, 1)
    with pytest.raises(ValueError, match='integer'):
        nonabelian.highest_Q_power(-1)
    with pytest.raises(ValueError, match='comparison margin'):
        nonabelian.completion_response_witness(F(1))
