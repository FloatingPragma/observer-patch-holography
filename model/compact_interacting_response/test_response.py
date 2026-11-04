from fractions import Fraction as F

import numpy as np
import pytest

from . import algebra, finite


def test_exact_spectrum_and_actual_charge_algebra():
    result = algebra.exact_identities()
    assert result['localized_annihilators'] == 180
    assert result['charge_3B_minus_3L_Fock_range'] == [-144, 108]


@pytest.mark.parametrize('values,regime', [((1., .5, .1, .2), 'paired'),
                                         ((.001, .01, .004, .005), 'blocked'),
                                         ((.6, .1, .01, 2.), 'single_odd')])
def test_full_car_spectrum_and_complete_readout(values, regime):
    row = finite.block_readout(*values)
    assert row['slice_ground_regime'] == regime
    assert len(row['complete_Fock_energies']) == 16
    assert abs(sum(float(x[1]) for x in row['addition_plus_removal_poles_and_weights'])-1) < 1e-10


def test_spatial_kinetic_term_changes_operator_and_ground():
    complete = finite.numerical_block(.001, .01, .004, .005)[0]
    dropped = finite.numerical_block(.001, .01, 0., 0.)[0]
    assert np.linalg.norm(complete-dropped) > .01
    assert finite.block_readout(.001, .01, 0., 0.)['slice_ground_regime'] == 'paired'
    assert finite.block_readout(.001, .01, .004, .005)['slice_ground_regime'] == 'blocked'


def test_mass_interval_encloses_retained_source_calculation():
    from model.mcpr_native_memory_bridge.construction import design
    for box, value in zip(algebra.mass_boxes(), design('forward', 0)['masses_over_v'], strict=True):
        assert float(box.lo)-1e-17 <= value <= float(box.hi)+1e-17
        assert box.hi-box.lo < F(1, 10**18)


def test_memory_norm_includes_nonzero_kinetics_and_180_modes():
    value = algebra.compact_memory_bound(1, F(1, 128), F(1, 12), 1, 2**30)
    expected = 90*F(7, 4)/60+180+F(8100, 128)
    assert F(value['core_Fock_norm_upper']) == expected
    assert 0 < F(value['complete_history_arbitrary_CQ_error_squared_upper']) < 1
    with pytest.raises(ValueError):
        algebra.compact_memory_bound(1., F(1, 128), F(1, 12), 1, 1024)


@pytest.mark.parametrize('values', [(.001, .01, .004, .005), (1., .5, .1, .2),
                                   (.6, .1, .01, 2.)])
def test_variational_vacuum_bound_against_complete_Fock_spectrum(values):
    m, u, e1, e2 = values
    quadratic_ground = -max(m, e1)-max(m, e2)
    even = int(m > e1)+int(m > e2)
    upper = min(quadratic_ground+u*even, -e1-e2)
    full = np.linalg.eigvalsh(finite.numerical_block(*values)[0])[0]
    assert quadratic_ground-1e-12 <= full <= upper+1e-12


@pytest.mark.parametrize('values', [(0, .1, .2, .3), (1, -.1, .2, .3),
                                    (1, .1, -.2, .3), (1, .1, np.nan, .3)])
def test_bad_physical_control_domains_rejected(values):
    with pytest.raises(ValueError):
        finite.numerical_block(*values)
