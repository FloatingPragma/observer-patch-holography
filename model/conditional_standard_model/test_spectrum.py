"""Independent complete-kernel, generic flavor and gauge normalization checks."""
from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from model.common_source import calibrated_matter as calibrated
from model.common_source.chiral_candidate import action, kinetic
from . import spectrum as module


def numeric(value):
    return np.asarray(value.evalf(30), complex)


@pytest.fixture(scope='module')
def catalogue():
    return module.build_catalogue()


def test_default_uses_actual_existing_K45_and_every_family_color(catalogue):
    assert catalogue['exact_same_calibrated_K45'] is True
    actual = module.mass_kernel()
    assert module._zero(actual-calibrated.mass_kernel())
    singular = np.sort(np.linalg.svd(numeric(actual), compute_uv=False))
    expected = [0.]*3
    for channel in module.CHANNELS:
        for mass in calibrated.mass_coordinates(channel):
            expected.extend([float(mass)]*(2 if channel == 'charged_lepton' else 6))
    np.testing.assert_allclose(singular, np.sort(expected), atol=1e-12, rtol=1e-10)
    assert catalogue['kernel_checks']['K45_rank_exact'] == 42
    assert catalogue['kernel_checks']['K45_nullity_exact'] == 3
    for channel, record in catalogue['families'].items():
        for found, expected_value in zip(record['singular_values_exact'], calibrated.yukawa_coordinates(channel), strict=True):
            assert s.simplify(s.sympify(found)-expected_value) == 0
    assert catalogue['counts']['internal_Weyl_components'] == 45
    assert catalogue['counts']['continuum_fermion_on_shell_states_including_antiparticles'] == 90
    assert catalogue['counts']['static_CAR_mode_count_is_identified_with_on_shell_degrees'] is False


def test_electroweak_mass_normalization_from_covariant_derivative_and_potential():
    data = replace(module.calibrated_inputs(), vEW=s.Rational(7), gW=s.Rational(3, 5),
                   gY=s.Rational(2, 7), lambda_SM=s.Rational(1, 11))
    result = module.bosonic_spectrum(data)
    w1, w2, w3, b = s.symbols('W1 W2 W3 B', real=True)
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    higgs = s.Matrix([0, data.vEW/s.sqrt(2)])
    covariant = (data.gW*sum((w*g/2 for w, g in zip((w1,w2,w3),pauli,strict=True)),s.zeros(2))
                 +data.gY*b*s.eye(2)/2)*higgs
    kinetic_mass_form = s.expand((covariant.H*covariant)[0])
    real_mass_squared = s.hessian(kinetic_mass_form,(w1,w2,w3,b))
    assert real_mass_squared[:2,:2] == s.eye(2)*data.gW**2*data.vEW**2/4
    assert real_mass_squared[2:,2:] == s.Matrix([[data.gW**2,-data.gW*data.gY],
                                                [-data.gW*data.gY,data.gY**2]])*data.vEW**2/4
    rows = {r['name']:r for r in result['rows']}
    assert s.sympify(rows['Z']['tree_mass_squared_exact']) == (data.gW**2+data.gY**2)*data.vEW**2/4
    assert s.sympify(rows['Higgs_radial']['tree_mass_squared_exact']) == 2*data.lambda_SM*data.vEW**2
    assert sum(r['multiplicity']*r['free_polarizations_per_species'] for r in result['rows']) == 28
    assert result['tree_m2_is_inferred_from_SMDR_loop_minimization'] is False


def test_complex_yukawas_full_kernel_and_independent_singular_values():
    inputs = replace(module.calibrated_inputs(), vEW=s.Rational(3),
        Yu=s.Matrix([[1,s.I,2],[0,2-s.I,1],[1,0,3]]),
        Yd=s.Matrix([[2,1+s.I,0],[s.I,3,1],[1,0,1-s.I]]),
        Ye=s.Matrix([[1,2,0],[0,s.I,1],[s.I,0,1]]))
    result = module.build_catalogue(inputs)
    full = np.linalg.svd(numeric(module.mass_kernel(inputs)),compute_uv=False)
    expected = [0.]*3
    for channel, y in inputs.yukawas().items():
        independent = np.sort(np.linalg.svd(numeric(y),compute_uv=False))
        reported = [float(s.N(s.sympify(v),30)) for v in result['families'][channel]['singular_values_exact']]
        np.testing.assert_allclose(reported,independent,rtol=3e-13,atol=2e-14)
        for singular in independent:
            expected.extend([float(inputs.vEW/s.sqrt(2))*singular]*(2 if channel=='charged_lepton' else 6))
    np.testing.assert_allclose(np.sort(full),np.sort(expected),rtol=5e-13,atol=5e-14)
    assert result['scope']['default_uses_existing_empirical_calibration'] is False


def test_exact_flavor_rotations_leave_masses_and_ratios_invariant():
    left=s.Matrix([[s.Rational(3,5),s.Rational(4,5),0],[-s.Rational(4,5),s.Rational(3,5),0],[0,0,s.I]])
    right=s.Matrix([[1,0,0],[0,s.Rational(5,13),s.I*s.Rational(12,13)],[0,s.I*s.Rational(12,13),s.Rational(5,13)]])
    assert left.H*left == right.H*right == s.eye(3)
    y=s.Matrix([[1,s.I,2],[0,2,1],[1,0,3]])
    first,_=module.singular_spectrum(y)
    second,_=module.singular_spectrum(left*y*right.H)
    assert first['Gram_characteristic_polynomial'] == second['Gram_characteristic_polynomial']
    assert first['singular_values_exact'] == second['singular_values_exact']


def test_zero_and_repeated_yukawa_spectra_retain_all_neutral_and_flavor_modes():
    inputs=replace(module.calibrated_inputs(),Yu=s.zeros(3),Yd=s.diag(0,2,2),Ye=s.diag(1,1,1))
    result=module.build_catalogue(inputs)
    assert result['kernel_checks']['K45_rank_exact'] == 18
    assert result['kernel_checks']['K45_nullity_exact'] == 27
    assert result['families']['up']['ratios_to_lightest_in_channel_exact'] is None
    assert result['families']['down']['ratios_to_lightest_in_channel_exact'] is None
    assert result['families']['charged_lepton']['ratios_to_lightest_in_channel_exact'] == ['1','1','1']
    assert result['counts']['continuum_fermion_on_shell_states_including_antiparticles'] == 90
    assert result['scope']['three_massless_neutrinos_fit_observed_oscillations'] is False


def test_mass_shell_uses_actual_complex_K45_not_only_a_mass_table():
    inputs=replace(module.calibrated_inputs(),vEW=2,Yu=s.diag(1,2,3),
                   Yd=s.Matrix([[0,s.I,0],[2,0,0],[0,0,3]]),Ye=s.diag(1,2,4))
    k=module.mass_kernel(inputs)
    symbol=kinetic.bdg_symbol((.2,-.3,.4),kernel=k,speed=1.)
    positive=np.repeat(np.sqrt(.29+np.linalg.svd(numeric(k),compute_uv=False)**2),2)
    np.testing.assert_allclose(np.linalg.eigvalsh(symbol),np.sort(np.r_[-positive,positive]),atol=1e-12,rtol=1e-12)
    charge=kinetic.nambu_charge(action.electric_charge())
    np.testing.assert_allclose(symbol@charge,charge@symbol,atol=1e-14,rtol=0)


@pytest.mark.parametrize('name',['vEW','gY','gW','g3','lambda_SM'])
@pytest.mark.parametrize('bad',[0,-1,True,.1,s.Float('.1'),s.Symbol('unknown'),s.oo])
def test_broken_domain_and_exact_inputs_are_required(name,bad):
    with pytest.raises(ValueError):
        replace(module.calibrated_inputs(),**{name:bad})


@pytest.mark.parametrize('bad',[s.zeros(2),[[1,2],[3,4]],s.eye(3)*s.Float('.1'),s.eye(3)*s.Symbol('unknown')])
def test_invalid_yukawa_matrices_rejected(bad):
    with pytest.raises(ValueError):
        replace(module.calibrated_inputs(),Yu=bad)


def test_catalogue_scope_and_default_calibration_do_not_hide_pole_or_source_selection(catalogue):
    assert catalogue['calibration']['calibration_contract']['all_nine_Yukawas_CKM_and_v_are_empirical_inputs']
    assert all(row['physical_pole_mass'] is None for row in catalogue['fermions']+catalogue['bosons']['rows'])
    assert catalogue['scope']['these_fields_are_identified_as_radial_topological_defects'] is False
    assert catalogue['scope']['vacuum_or_defect_sector_selected'] is False
    assert catalogue['finite_representation']['gauge_current_generators'] == 12
    assert sum(row['dimension_per_family'] for row in catalogue['finite_representation']['multiplets_in_all_left_convention']) == 15
