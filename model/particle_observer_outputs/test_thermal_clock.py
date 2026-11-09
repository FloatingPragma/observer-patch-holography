"""Full-sector thermal, modular and complete-record checks."""
from fractions import Fraction as F

import pytest
import sympy as s

from model.particle_observer_outputs import thermal_clock as t
from model.particle_observer_outputs import higgs_records as h


@pytest.fixture(scope='module')
def thermal():
    return t.thermal_fixture()


def test_exact_Gibbs_weights_retain_every_dark_state(thermal):
    assert thermal['source']['H'].rows == 307
    assert thermal['ranks'] == dict(ground=1, excited=1, dark=305)
    assert thermal['energies'] == dict(ground=-7, excited=8, dark=1)
    assert thermal['partition'] == s.Rational(71809, 256)
    assert thermal['per_state_weights'] == dict(
        ground=s.Rational(32768, 71809), excited=s.Rational(1, 71809), dark=s.Rational(128, 71809))
    assert s.trace(thermal['rho']) == 1
    assert min(thermal['per_state_weights'].values()) == s.Rational(1, 71809)
    assert sum(thermal['per_state_weights'].values()) != 1
    assert h._zero(thermal['rho']*thermal['source']['H']-thermal['source']['H']*thermal['rho'])


def test_modular_log_identity_covers_the_actual_full_matrix(thermal):
    difference = (thermal['modular_hamiltonian']-thermal['beta']*thermal['source']['H']
                  -s.log(thermal['partition'])*s.eye(307))
    assert h._zero(difference)
    report = t.modular_clock_certificate()
    assert report['same_operator_clock_map'] == 'tau=-beta*s'
    assert sum(row['matrix_unit_count'] for row in report['spectral_operator_blocks']) == 307**2
    row = next(x for x in report['spectral_operator_blocks']
               if x['left_block'] == 'ground' and x['right_block'] == 'excited')
    assert row['Heisenberg_frequency'] == '-15'
    assert s.simplify(s.sympify(row['modular_frequency'])-15*s.log(2)) == 0
    # A plus sign in tau=+beta*s would reverse this nonzero frequency.
    assert s.sympify(row['modular_frequency']) != thermal['beta']*s.sympify(row['Heisenberg_frequency'])


def test_full_canonical_purification_and_preparation(thermal):
    data = t.purification_data()
    c, blank = data['coefficient'], data['blank_coefficient']
    difference, normalization = data['householder_difference'], data['householder_norm_squared']
    assert data['joint_dimension'] == 94249 and data['Schmidt_rank'] == 307
    assert data['charged_reference_required_for_full_sector_Gauss_route']
    assert not data['field_density_supported_only_on_singlets']
    assert h._zero(c*c.H-thermal['rho'])
    assert h._zero(c.T*c.conjugate()-thermal['rho'].T)
    # Full-matrix marginal identities above fix ||C||; the complete vector
    # norm follows by polarization without expanding nested radicals.
    polarization = s.trace(thermal['rho'])+s.trace(blank.H*blank)-2*s.trace(blank.H*c)
    assert t._clean(polarization-normalization) == 0
    assert h._zero(blank-2*t._clean(s.trace(difference.H*blank)/normalization)*difference-c)
    # For ||u||=1, the all-column reflection squares to I on the complete
    # joint carrier, not only on the prepared input.
    assert t._clean(-4+4*polarization/normalization) == 0


def test_conjugate_reference_covariance_and_admission(thermal):
    # The same-H flow is a concrete member of the general commuting-symmetry
    # domain. Gauge membership is supplied by the source Ward result.
    report = t.purification_covariance(s.I*thermal['source']['H'])
    assert report['joint_generator_annihilates_purification']
    with pytest.raises(ValueError):
        t.purification_covariance(thermal['source']['H'])
    with pytest.raises(ValueError):
        t.purification_covariance(s.I*thermal['source']['vacuum_projector'])


def test_complete_occupation_weights_and_exact_probabilities(thermal):
    weights = thermal['occupation_weights']
    assert weights == dict(zip(('empty', 'bright', 'dark'),
                              (s.Rational(x, 1077135) for x in (262151, 229384, 585600))))
    assert sum(weights.values()) == 1 and weights['dark'] > s.Rational(1, 2)
    one, two = t.thermal_joint_law(s.pi/45), t.thermal_joint_law(2*s.pi/45)
    assert len(one) == len(two) == 4
    assert all(p > 0 for p in (*one.values(), *two.values()))
    assert sum(one.values()) == sum(two.values()) == 1
    assert one['occupied', 'occupied'] >= weights['dark']
    assert two['empty', 'occupied']/(4*one['empty', 'occupied']) == s.Rational(3, 4)


def test_joint_law_against_actual307_state_Kraus_operators(thermal):
    # Integer spectral phases make this independent full-matrix replay exact
    # without numerical exponentiation. This time is not an inversion window.
    data = thermal['source']
    direct = {key: t._clean(s.trace(k*thermal['rho']*k.H))
              for key, k in h.sequential_kraus(data, s.pi).items()}
    expected = t.thermal_joint_law(s.pi)
    assert direct == expected
    assert sum(direct.values()) == 1


def test_gap_readout_needs_reproducible_preparation_not_a_ground_state(thermal):
    data, time = thermal['source'], s.pi/45
    for weights in (dict(empty=1, bright=0, dark=0),
                    dict(empty=s.Rational(1, 3), bright=s.Rational(1, 6), dark=s.Rational(1, 2)),
                    thermal['occupation_weights']):
        first = h.occupation_joint_law(data, time, **weights)
        second = h.occupation_joint_law(data, 2*time, **weights)
        p, q = first['empty', 'occupied'], second['empty', 'occupied']
        assert h.gap_interval_from_records((p, p), (q, q), time)['gap_interval'] == (15, 15)
    assert not t.preparation_independent_gap_identity()['ground_preparation_required']
    dark = h.occupation_joint_law(data, time, empty=0, bright=0, dark=1)
    assert dark == {('empty', 'empty'): 0, ('empty', 'occupied'): 0,
                    ('occupied', 'empty'): 0, ('occupied', 'occupied'): 1}
    with pytest.raises(ValueError):
        h.gap_interval_from_records((0, 0), (0, 0), time)
    one = h.occupation_joint_law(data, time, empty=1, bright=0, dark=0)
    changed = h.occupation_joint_law(data, 2*time, empty=s.Rational(1, 2), bright=s.Rational(1, 2), dark=0)
    assert changed['empty', 'occupied']/(4*one['empty', 'occupied']) == s.Rational(3, 8)


def test_reference_resources_and_complete_error_are_explicit():
    report = t.preparation_plan()
    assert report['spin_bits'] == 18
    assert report['joint_code_dimension'] == 94249
    assert report['unitary_dimension'] == report['padded_unitary_dimension'] == 512**2
    assert report['reference_retained_after_preparation']
    assert not report['coherent_compiled_tape_executed_here']
    assert not report['passive_repair_Gibbs_preparation_claimed']
    receipt = t.record_certificate()
    assert F(receipt['two_tick_two_detector_complete_error']) == F(5, 10**6)
    assert receipt['recovered_gap'] == '15'
    assert t.certificate()['empirical_mass_inputs'] == []
