import pytest
import sympy as s

from .spectrum import symbol,symbol_lattice,energy,projectors,continuum_bound,source_mass
from .source import certificate as source_certificate
from .certificate import algebra_certificate,continuum_certificate
from .position import certificate as position_certificate
from .fock import rest_instrument,act


def test_literal_gauss_source_and_all_occupations():
    result=source_certificate()
    assert result['full_local_graded_occupation_dimension']==1048576
    assert result['literal_slot_charge_equals_occupation']
    assert not result['bosonic_link_tensor_product_alone_proves_local_CAR']


def test_spectrum_car_and_charge_conjugation():
    result=algebra_certificate()
    assert result['twenty_mode_CAR_column_checks']==4800
    assert result['all_four_mode_exterior_lift_dimensions']==[1,4,6,4,1]
    assert len(result['Floquet_even_corner_small_phase_branches'])==4


def test_complete_filled_sea_particle_hole_instrument():
    result=rest_instrument()
    rows=result['complete_preparation_and_particle_hole_records']
    assert len(rows)==36 and sum(s.Rational(row['probability']) for row in rows)==1
    assert result['retained_zero_probability_records']==30
    successes=[row for row in rows if row['probability']!='0' and row['preparation_record'] in ('add_from_reservoir','remove_to_reservoir')]
    assert len(successes)==2
    assert {(row['relative_pair_charge'],row['excitation_energy'],row['probability']) for row in successes}=={(1,'2','1/4'),(-1,'2','1/4')}


def test_literal_spatial_source_and_continuum_bounds():
    assert position_certificate()['literal_source_Bloch_columns_checked']==20
    result=continuum_certificate()
    assert result['one_finite_energy_momentum_valley_per_face_channel']
    assert [row['positive_energy'] for row in result['Wilson_corners_at_spacing_one_half']].count('2')==1


def test_clock_and_momentum_input_contract():
    a=s.Rational(1,2);p=(s.pi,0,0)
    assert symbol(a,p)==symbol_lattice(a,tuple(a*x for x in p))
    assert energy(a,(0,0,0))==source_mass()==2
    plus,minus=projectors(a,(0,0,0))
    assert plus+minus==s.eye(4)
    for invalid in (0,-1,0.5):
        with pytest.raises(ValueError):symbol(invalid,p)
    with pytest.raises(ValueError):symbol(a,(0,0))
    with pytest.raises(ValueError):continuum_bound(a,-1)
    with pytest.raises(ValueError):act(1<<20,0,True,20)
