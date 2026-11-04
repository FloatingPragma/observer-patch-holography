"""Full matrix, CAR and arbitrary-Higgs controls for the Majorana bridge."""
from dataclasses import replace
from itertools import product

import numpy as np
import pytest
import sympy as s

from ..common_source import calibrated_matter as calibrated
from ..common_source.chiral_candidate import action, kinetic
from ..conditional_standard_model import spectrum as charged
from . import spectrum as module


def numeric(value):
    return np.asarray(value.evalf(30), complex)


def matter():
    return charged.ActionInputs(2, 1, 1, 1, 1,
                                s.diag(1, 2, 3), s.Matrix([[0, s.I, 0], [2, 0, 0], [0, 0, 3]]),
                                s.diag(1, 3, 4), ('exact finite control',))


def neutrinos(masses=(1, 2, 4)):
    first = s.Matrix([[s.Rational(3, 5), s.Rational(4, 5), 0],
                      [-s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, (3+4*s.I)/5]])
    second = s.Matrix([[1, 0, 0], [0, s.Rational(5, 13), s.I*s.Rational(12, 13)],
                       [0, s.I*s.Rational(12, 13), s.Rational(5, 13)]])
    return module.MajoranaInputs(tuple(map(s.Rational, masses)), (first, second), s.sqrt(2),
                                 ('exact complex Takagi control',))


def test_calibrated_neutral_block_resolves_small_masses_without_large_quark_spectrum():
    result = module.catalogue()
    assert result['same_existing_calibrated_K45_exact']
    assert (result['neutrino_rank'], result['full_K45_rank'], result['full_K45_nullity']) == (2, 44, 1)
    indices = action.ordering()['massless_neutrino_internal_family_indices']
    actual = module.full_kernel().extract(indices, indices)
    assert charged._zero(actual-calibrated.neutrino_family_kernel())
    p = (2e-11, -3e-11, 4e-11)
    # Scaling the isolated exact invariant block prevents the O(100GeV)
    # charged sector from hiding O(1e-11GeV) neutrino frequencies.
    block = kinetic.bdg_symbol(p, kernel=actual, speed=1.)/1e-11
    masses = np.array([float(x)/1e-11 for x in calibrated.neutrino_mass_coordinates()])
    expected = np.repeat(np.sqrt(29+masses*masses), 2)
    np.testing.assert_allclose(np.linalg.eigvalsh(block), np.sort(np.r_[-expected, expected]), atol=2e-13, rtol=2e-13)
    for mass, bounds in zip(calibrated.neutrino_mass_coordinates(), result['masses_eV_exact_intervals'], strict=True):
        lo, hi = map(s.Rational, bounds)
        assert lo*lo <= s.simplify((mass*10**9)**2) <= hi*hi
        assert hi-lo <= s.Rational(1, 10**14)


@pytest.mark.parametrize('masses,rank', [((0, 0, 0), 0), ((0, 0, 4), 1), ((0, 2, 4), 2), ((1, 2, 4), 3), ((2, 2, 2), 3)])
def test_every_neutrino_rank_and_repeated_mass_on_actual_K45(masses, rank):
    data = neutrinos(masses)
    result = module.catalogue(data, matter=matter())
    kernel = module.full_kernel(matter=matter(), neutrinos=data)
    assert kernel.rank() == 42+rank
    assert result['full_K45_rank'] == 42+rank
    assert result['full_K45_nullity'] == 3-rank
    assert sum(row['massless_Weyl'] for row in result['neutrino_rows']) == 3-rank
    assert all(row['physical_pole_mass'] is None for row in result['neutrino_rows'])


def test_actual_full180_momentum_symbol_for_complex_masses_and_every_channel():
    inputs, data = matter(), neutrinos()
    kernel = module.full_kernel(matter=inputs, neutrinos=data)
    expected_masses = []
    for name, y in inputs.yukawas().items():
        for value in np.linalg.svd(numeric(y), compute_uv=False):
            expected_masses.extend([float(inputs.vEW/s.sqrt(2))*value]*(2 if name == 'charged_lepton' else 6))
    expected_masses.extend(map(float, data.masses))
    assert len(expected_masses) == 45
    for p in ((0., 0., 0.), (.2, -.3, .7), (-.5, .25, -.125)):
        symbol = module.nambu_symbol(p, matter=inputs, neutrinos=data)
        positive = np.repeat(np.sqrt(sum(x*x for x in p)+np.array(expected_masses)**2), 2)
        np.testing.assert_allclose(np.linalg.eigvalsh(symbol), np.sort(np.r_[-positive, positive]), atol=2e-12, rtol=2e-12)
        assert symbol.shape == (180, 180)
    assert module.universal_symbol_certificate()['same_identity_for_every_symmetric_internal_dimension']


def test_same_Higgs_scales_Majorana_quadratically_and_charged_masses_linearly():
    inputs, data = matter(), neutrinos()
    h = (s.I, 2-s.I)
    k = module.full_kernel(h, matter=inputs, neutrinos=data)
    neutrino = module.majorana_kernel(h, data)
    scaled = module.full_kernel(tuple(2*x for x in h), matter=inputs, neutrinos=data)
    assert charged._zero(scaled-(2*(k-neutrino)+4*neutrino))
    # A common Higgs phase doubles in the Majorana coefficient.
    assert charged._zero(module.majorana_kernel(tuple(s.I*x for x in h), data)+neutrino)


def test_arbitrary_Higgs_full_action_covariance_under_all_retained_currents():
    inputs, data = matter(), neutrinos()
    higgs = s.Matrix([1+s.I, 2-s.I])
    kernel = module.full_kernel(tuple(higgs), matter=inputs, neutrinos=data)
    for g, primitive in zip(action.generators(), action._representation()[1], strict=True):
        tangent = primitive[3:, 3:]*higgs
        # The complete action is at most quadratic in H and conjugate(H).
        derivative = (module.full_kernel(tuple(higgs+tangent), matter=inputs, neutrinos=data)
                      -module.full_kernel(tuple(higgs-tangent), matter=inputs, neutrinos=data))/2
        assert charged._zero(g.T*kernel+kernel*g+derivative)


def test_independent_complete64_state_CAR_spectrum_with_complex_pairing():
    data = neutrinos()
    k = data.family_kernel()
    delta = numeric(s.kronecker_product(action.SPIN_EPSILON, k.conjugate()))
    annihilators = []
    for i in range(6):
        a = np.zeros((64, 64), complex)
        for mask in range(64):
            if mask >> i & 1:
                a[mask ^ (1 << i), mask] = (-1)**((mask & ((1 << i)-1)).bit_count())
        annihilators.append(a)
    creation = sum((delta[i, j]*annihilators[i].conj().T@annihilators[j].conj().T/2
                    for i in range(6) for j in range(6)), np.zeros((64, 64), complex))
    actual = np.linalg.eigvalsh(creation+creation.conj().T)
    energies = sorted(sum(n*m for n, m in zip(row, (1, 1, 2, 2, 4, 4), strict=True))-7
                      for row in product((0, 1), repeat=6))
    np.testing.assert_allclose(actual, energies, atol=2e-12, rtol=2e-12)
    assert charged._zero(data.mixing().T*k*data.mixing()-s.diag(*data.masses))
    assert not charged._zero(data.mixing().T*k.conjugate()*data.mixing()-s.diag(*data.masses))


def test_actual_PMNS_support_has_exact_positive_phase_independent_bounds():
    record = module.mixing_support_certificate()
    matrix = numeric(calibrated.pmns_matrix())
    lower = np.array([[float(s.Rational(x)) for x in row] for row in record['squared_modulus_lower']])
    assert np.min(lower) > 0
    assert np.all(abs(matrix)**2+1e-14 >= lower)
    assert record['arbitrary_unitary_mixing_has_this_support'] is False


@pytest.mark.parametrize('change', ['negative_mass', 'float_mass', 'wrong_count', 'nonunitary', 'zero_higgs', 'no_provenance'])
def test_invalid_Takagi_input_is_rejected(change):
    values = dict(masses=(s.Rational(1), s.Rational(2), s.Rational(3)),
                  unitary_factors=(s.eye(3),), reference_higgs_norm=s.Rational(1), provenance=('control',))
    if change == 'negative_mass': values['masses'] = (-1, 2, 3)
    elif change == 'float_mass': values['masses'] = (0., 2, 3)
    elif change == 'wrong_count': values['masses'] = (1, 2)
    elif change == 'nonunitary': values['unitary_factors'] = (2*s.eye(3),)
    elif change == 'zero_higgs': values['reference_higgs_norm'] = 0
    else: values['provenance'] = ()
    with pytest.raises(ValueError): module.MajoranaInputs(**values)
