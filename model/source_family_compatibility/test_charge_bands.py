from fractions import Fraction as F
import pytest
import sympy as s
from . import charge_bands as b


def test_readable_bulk_band_preserves_original_preparation_weight():
    row=b.broad_band(F(1,8),1,6)
    assert row['each_original_preparation_weight']==F(1,2)
    assert row['each_original_second_moment']==F(35,256)
    assert row['positive_energy_window_squared']==[F(35,256),F(2744,5)]
    assert row['each_original_band_probability_lower']==F(25,100352)
    assert not row['band_is_a_spectral_gap_or_particle_mass']


def test_neighbour_degrees_affect_fourth_moment_not_local_second_moment():
    local=b.broad_band(F(1,8),1,1)
    boundary=b.broad_band(F(1,8),1,1,6)
    assert local['each_original_second_moment']==boundary['each_original_second_moment']
    assert local['each_original_band_probability_lower']>boundary['each_original_band_probability_lower']


@pytest.mark.parametrize('args',[(1,0,6),(F(1,8),1,6),(F(3,2),4,1)])
def test_band_probability_is_positive_without_strong_coupling_projection(args):
    row=b.broad_band(*args)
    low,high=row['positive_energy_window_squared']
    assert 0<low<high
    assert 0<row['each_original_band_probability_lower']<F(1,2)


@pytest.mark.parametrize('args',[(0,1,6),(1,-1,6),(True,1,6),(0.125,1,6),
                                 (1,1,0),(1,1,True),(1,1,6,5)])
def test_invalid_or_inexact_source_parameters_rejected(args):
    with pytest.raises(ValueError):b.broad_band(*args)


def test_dimer_full_source_moments_and_composite_clock():
    row=b.dimer()
    moments=[s.Matrix([[s.sympify(x) for x in r] for r in m])
             for m in row['source_retarded_moments_0_through_4']]
    h=moments[1]
    assert moments[2]==h*h+s.eye(4)/4
    assert moments[3]==h**3-h/4
    assert moments[3]!=h**3+h/4
    assert row['charged_second_moment_sum_rule_on_whole_Fock_space']
    assert row['actual_composite_clock']=='B=-3h'
    assert row['second_response_closes_on_this_dimer']
    assert not row['Dirac_composite_clock_B_equals_minus_h']


def test_dimer_complete_charge_energy_measurement_keeps_zero_records():
    row=b.dimer();records=row['complete_records']
    assert len(records)==48
    assert sum(s.sympify(r['original_probability']) for r in records)==1
    assert sum(r['original_probability']=='0' for r in records)==32
    for mode in range(4):
        for charge in ('addition','removal'):
            weight=sum(s.sympify(r['original_probability']) for r in records
                       if r['mode']==mode and r['charge']==charge)
            assert weight==s.Rational(1,8)
    assert not row['spatial_continuum_or_measured_mass_assigned']


def test_second_parameter_pair_is_recomputed_from_source():
    row=b.dimer(F(3,2),4)
    assert list(map(s.sympify,row['source_positive_charge_energies']))==[s.Rational(7,4),s.Rational(13,4)]
    assert list(map(s.sympify,row['each_mode_each_charge_original_weights']))==[s.Rational(2,5),s.Rational(1,10)]


def test_certificate_distinguishes_original_bands_from_physical_particles():
    row=b.certificate()
    assert row['finite_and_infinite_ground_same_original_weight_bound']
    assert not row['charged_pole_mass_or_Lorentz_covariance_derived']
