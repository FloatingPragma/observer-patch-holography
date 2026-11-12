import math
from fractions import Fraction as F

import pytest

from .thermodynamics import (species,massless_limit,integrals,gas,bessel_control,
                             finite_source_gibbs)


def test_complete_calibrated_inventory_and_exact_relativistic_counts():
    rows=species();assert len(rows)==17
    assert sum(x['degrees'] for x in rows if x['fermion'])==90
    assert sum(x['degrees'] for x in rows if not x['fermion'])==28
    result=massless_limit()
    assert result['energy_and_entropy_degrees']=='427/4'
    assert result['quark_number_fraction']=='108/191'
    assert result['lepton_number_fraction']=='27/191'
    assert F(result['quark_energy_fraction'])==F(252,427)
    assert result['extra_sterile_neutrino_degrees']==0


def test_entropy_identity_against_independent_temperature_derivative():
    for fermion in (True,False):
        for ratio in (.1,1.,10.):
            t=1/ratio;h=t*1e-4
            pressure=lambda temp:temp**4*integrals(1/temp,fermion)['pressure']
            derivative=(-pressure(t+2*h)+8*pressure(t+h)-8*pressure(t-h)+pressure(t-2*h))/(12*h)
            expected=t**3*integrals(ratio,fermion)['entropy']
            assert abs(derivative-expected)<max(1e-12,abs(expected)*2e-8)
            assert bessel_control(ratio,fermion)['numerical_quadrature_disagreement']<2e-10


def test_actual_finite_global_cap_prevents_independent_bose_replacement():
    result=finite_source_gibbs(10.)
    records=result['all_boson_occupation_records']
    assert len(records)==1653 and abs(sum(r['probability'] for r in records)-1)<2e-14
    assert result['boson_mean_number']<2<result['uncapped_same_mode_mean_boson_number']
    assert len(result['complete_factorized_fermion_log_records'])==180
    assert all(math.isfinite(r['log_probability']) for r in records)


def test_conditional_thermodynamic_consistency_and_domain_controls():
    for t in (.001,1.,10.,1000.):
        result=gas(t)
        assert 0<result['sound_speed_squared']<=1/3+1e-12
        assert not result['physical_temperature_evolution_predicted']
        assert abs(sum(r['number_fraction'] for r in result['species'])-1)<1e-12
        assert result['thermodynamic_residuals']['trace_anomaly_identity']<1e-8
    suppressed=integrals(10000.,False)
    assert suppressed['number']==0 and suppressed['each_dimensionless_integral_log_error_upper']<0
    assert suppressed['underflowed_upper_retained_as_log']
    with pytest.raises(ValueError):gas(0)
    with pytest.raises(ValueError):gas(float('nan'))
    with pytest.raises(ValueError):integrals(-1.,False)
    with pytest.raises(ValueError):finite_source_gibbs(float('nan'))
