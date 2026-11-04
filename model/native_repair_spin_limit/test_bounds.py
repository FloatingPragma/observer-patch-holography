"""Rate, time and domain controls for the analytic local theorem."""
from fractions import Fraction as F

import pytest

from . import bounds as b


def test_fixed_positive_rate_and_growing_archive_have_rational_small_error():
    fixed = b.repair_bounds(2**256, 2**512, support=8, observed_edges=12)
    growing = b.repair_bounds(2**768, 2**1536, time=2**64, support=8, observed_edges=12)
    for value in (fixed, growing):
        assert value['repair_rate_per_edge'] == value['exchange'] == 1
        assert 0 < value['faithful_input_and_bare_readout_joint_error'] < F(1, 10**6)
        assert not value['depends_on_total_volume']
        assert not value['rate_tends_to_zero']
    assert b.growing_time_exponents(F(1, 12))['quantum_error_power'] == F(-1, 2)


def test_fixed_exchange_error_decays_with_c_while_penalty_ratio_diverges():
    first = b.repair_bounds(2**128, 2**256, positive_energy_density=0)
    second = b.repair_bounds(2**129, 2**258, positive_energy_density=0)
    ratio = second['bare_input_and_readout_joint_error']/first['bare_input_and_readout_joint_error']
    assert F(49, 100) < ratio <= F(1, 2)
    assert second['penalty_over_exchange'] == 4*first['penalty_over_exchange']


def test_local_record_readout_cost_depends_on_observed_edges_but_not_total_volume():
    no_record = b.repair_bounds(2**128, 2**256, observed_edges=0)
    one = b.repair_bounds(2**128, 2**256, observed_edges=1)
    seven = b.repair_bounds(2**128, 2**256, observed_edges=7)
    assert no_record['complete_record_error_uncapped'] == 0
    assert seven['complete_record_error_uncapped'] == 7*one['complete_record_error_uncapped']
    assert seven['dressed_quantum_error_uncapped'] == no_record['dressed_quantum_error_uncapped']


def test_energy_domain_covers_bare_ground_and_thermal_inputs_without_stationarity_assumption():
    ground = b.repair_bounds(2**128, 2**256, positive_energy_density=0)
    thermal = b.repair_bounds(2**128, 2**256, positive_energy_density=2)
    assert ground['bare_defect_density_upper'] < thermal['bare_defect_density_upper']
    assert thermal['bare_readout_joint_error'] < F(1, 10**6)
    assert b.sqrt_upper(F(1, 9)) == F(1, 3)
    assert b.sqrt_upper(F(2, 7))**2 >= F(2, 7)


@pytest.mark.parametrize('kwargs', [
    {'c': 1, 'U': 1}, {'c': True}, {'c': 1.}, {'rate': 0},
    {'time': -1}, {'support': True}, {'observed_edges': -1},
    {'observed_edges': True}, {'positive_energy_density': -1},
])
def test_inputs_outside_the_declared_exact_domain_are_rejected(kwargs):
    args = {'c': 2**128, 'U': 2**256, **kwargs}
    with pytest.raises(ValueError):
        b.repair_bounds(**args)


def test_expanding_time_claim_does_not_include_the_uncontrolled_endpoint():
    with pytest.raises(ValueError):
        b.growing_time_exponents(F(1, 6))
