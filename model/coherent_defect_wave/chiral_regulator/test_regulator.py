import sympy as s
import pytest
from .model import (symbol,spectral_momenta,momentum_matrix,projectors,spectral_certificate,
                    position_certificate,filled_reference_certificate,instrument_certificate,tail_bound)


def test_literal_spectral_derivative_and_nonlocality():
    result=position_certificate()
    assert result['literal_side_two_cubic_columns_checked']==16
    assert result['finite_nonlocality_witness']['momentum_matrix_element']=='-I/2'
    assert result['finite_strict_light_cone_proved'] is False


def test_chiral_particle_and_hole_are_not_same_charge_mirrors():
    result=spectral_certificate()
    assert result['particle_helicity']=='-1/2' and result['antiparticle_helicity']=='1/2'
    assert result['hole_is_independent_same_charge_right_Weyl_mirror'] is False
    assert s.sympify(result['Tr_Y_power5'])==-s.sympify(result['Tr_conjugate_Y_power5'])


def test_unique_filled_reference_has_only_a_volume_gap():
    result=filled_reference_certificate()
    assert result['gap_vanishes_in_infinite_volume']
    assert result['infinite_volume_rest_mass']=='0'
    assert [row['occupied_negative_modes'] for row in result['rows']]==[120,960,3240]


def test_complete_free_records_retain_failure_created_pair():
    result=instrument_certificate()
    assert result['complete_records']==16 and result['nonzero_records']==6
    populated=sorted(s.sympify(row['probability']) for row in result['records'] if row['probability']!='0')
    assert populated==sorted([(3-s.sqrt(3))/12,(2+s.sqrt(3))/12,s.Rational(1,12),
                              (3+s.sqrt(3))/12,(2-s.sqrt(3))/12,s.Rational(1,12)])
    assert sum(s.sympify(row['probability']) for row in result['records'] if row['energy']=='sqrt(3)')==s.Rational(1,6)


def test_exact_domains_and_cutoff_tail():
    assert tail_bound(3,2,10)==1
    assert symbol((1,2,2))**2==9*s.eye(2)
    assert set(spectral_momenta(4))=={-s.Rational(3,2),-s.Rational(1,2),s.Rational(1,2),s.Rational(3,2)}
    for side in (0,1,3,4.0):
        with pytest.raises(ValueError):spectral_momenta(side)
    with pytest.raises(ValueError):symbol((1.0,2,2))
    with pytest.raises(ValueError):projectors((0,0,0))
    with pytest.raises(ValueError):tail_bound(0,2,1)
