"""Complete-sector source identities, record inversion and dark controls."""
import pytest
import sympy as s

from model.particle_observer_outputs import higgs_records as h


@pytest.fixture(scope='module')
def data():
    return h.sector_data(s.diag(1, 2, 3), 1)


def test_actual_full_sector_spectral_resolution(data):
    assert data['H'].shape == (307, 307)
    assert data['g_squared'] == 56
    assert data['gap'] == 15
    assert (data['ground_energy'], data['excited_energy'], data['dark_energy']) == (-7, 8, 1)
    projectors = h.spectral_projectors(data)
    assert [s.trace(p) for p in projectors.values()] == [1, 1, 305]
    assert h._zero(sum(projectors.values(), s.zeros(307))-s.eye(307))
    for name, p in projectors.items():
        assert p == p.H
        assert h._zero(p*p-p)
        for other_name, q in projectors.items():
            if name != other_name:
                assert h._zero(p*q)
    assert h._zero(data['H']+7*projectors['ground']-8*projectors['excited']-projectors['dark'])


def test_same_vertex_phase_pulse_prepares_ground_without_discard(data):
    u = h.preparation_unitary(data)
    assert h._zero(u.H*u-s.eye(307))
    assert h._zero(u*data['vacuum_projector']*u.H-data['ground_projector'])
    assert h._zero(u*data['dark_projector']-data['dark_projector'])
    # V(iY) is the actual source phase-shifted vertex, not a fitted 2x2 matrix.
    j, basis = h.source.zero_charge_block(s.I*data['Y'], pair_cost=0)
    assert basis == data['basis']
    v, omega = data['bright_vector'], data['vacuum']
    assert h._zero(j-s.I*(v*omega.H-omega*v.H))
    scalar_phase = data['vacuum_projector']+s.I*data['occupation']
    assert h._zero(scalar_phase*data['vertex']*scalar_phase.H-j)
    angle = h.preparation_angle(data)
    assert s.simplify(s.cos(angle)**2-s.Rational(8, 15)) == 0


def test_propagator_has_one_time_and_preserves_dark_phase(data):
    assert h._zero(h.propagator(data, 0)-s.eye(307))
    u = h.propagator(data, s.pi)
    assert h._zero(u.H*u-s.eye(307))
    assert h._zero(u*data['dark_projector']+data['dark_projector'])
    assert h._zero(u*data['ground_projector']+data['ground_projector'])
    assert h._zero(u*data['excited_projector']-data['excited_projector'])


@pytest.mark.parametrize('refined,expected_count', [(False, 4), (True, 9)])
def test_every_pointer_outcome_is_retained_and_complete(data, refined, expected_count):
    maps = h.sequential_kraus(data, s.pi, refined)
    assert len(maps) == expected_count
    total = sum((k.H*k for k in maps.values()), s.zeros(307))
    assert h._zero(total-s.eye(307))
    if refined:
        for key in (('empty', 'dark'), ('dark', 'empty'), ('bright', 'dark'), ('dark', 'bright')):
            assert key in maps
            assert h._zero(maps[key])
        assert h._zero(maps['dark', 'dark']+data['dark_projector'])


def test_joint_formula_matches_the_complete_quantum_instrument(data):
    maps = h.sequential_kraus(data, s.pi)
    dark = h.gauge_invariant_dark_probe(data)['density']
    states = [data['ground_projector'], data['vacuum_projector'], dark,
              (data['vacuum_projector']+dark)/2]
    for rho in states:
        weights = {name: h._clean(s.trace(p*rho))
                   for name, p in h.measurement_projectors(data, True).items()}
        expected = h.occupation_joint_law(data, s.pi,
            empty=weights['empty'], bright=weights['bright'], dark=weights['dark'])
        direct = {key: h._clean(s.trace(k*rho*k.H)) for key, k in maps.items()}
        assert direct == expected
        assert sum(direct.values()) == 1


def test_same_source_ground_spectral_residue(data):
    pminus, pplus, dark = h.spectral_projectors(data).values()
    n = data['occupation']
    assert s.trace(pminus*n) == s.Rational(7, 15)
    assert s.trace(pminus*n*pplus*n) == s.Rational(56, 225)
    assert h._zero(dark*n*pminus)
    response = h.ground_correlation(data, s.pi/45)
    assert response['mean'] == s.Rational(7, 15)
    assert response['residue'] == s.Rational(56, 225)
    # The ordered quantum correlator is complex; it is not the classical
    # probability of two occupied projective records.
    assert s.simplify(s.im(response['ordered'])) != 0
    assert h.ground_joint_law(data, s.pi/45)['occupied', 'occupied'].is_real


def test_unconditional_records_recover_the_gap_without_mass_input(data):
    time = s.pi/45
    p = h.ground_joint_law(data, time)['empty', 'occupied']
    q = h.ground_joint_law(data, 2*time)['empty', 'occupied']
    assert (p, q) == (s.Rational(448, 3375), s.Rational(1344, 3375))
    inverse = h.gap_interval_from_records((p, p), (q, q), time)
    assert inverse['ratio_interval'] == (s.Rational(3, 4), s.Rational(3, 4))
    assert inverse['gap_interval'] == (15, 15)
    assert not inverse['conditioning_on_first_record']
    error = s.Rational(1, 100000)
    uncertain = h.gap_interval_from_records((p-error, p+error), (q-error, q+error), time)
    low, high = uncertain['ratio_interval']
    assert low < s.Rational(3, 4) < high


def test_bright_only_truncation_fails_on_a_gauge_invariant_input(data):
    probe = h.gauge_invariant_dark_probe(data)
    x, vector, rho = probe['family_probe'], probe['vector'], probe['density']
    assert s.trace(data['Y'].H*x) == 0
    assert (vector.H*vector)[0] == 1
    assert h._zero(data['H']*vector-vector)
    kept = data['vacuum_projector']+data['bright_projector']
    assert s.trace(kept*rho*kept) == 0
    law = h.occupation_joint_law(data, s.pi/45, empty=0, bright=0, dark=1)
    assert law == {('empty', 'empty'): 0, ('empty', 'occupied'): 0,
                   ('occupied', 'empty'): 0, ('occupied', 'occupied'): 1}
    # Postselected bright normalization of this equal mixture moves half
    # its total original weight and removes a physical neutral input.
    true_law = h.occupation_joint_law(data, s.pi, empty=s.Rational(1, 2), bright=0, dark=s.Rational(1, 2))
    false_law = h.occupation_joint_law(data, s.pi, empty=1, bright=0, dark=0)
    assert sum(abs(true_law[k]-false_law[k]) for k in true_law)/2 == s.Rational(1, 2)


def test_entangled_untouched_reference_is_retained(data):
    dark = h.gauge_invariant_dark_probe(data)['vector']
    purification = s.ImmutableSparseMatrix.hstack(data['vacuum'], dark)/s.sqrt(2)
    maps = h.sequential_kraus(data, s.pi, refined=True)
    reduced_reference = sum(((k*purification).H*(k*purification) for k in maps.values()), s.zeros(2))
    assert h._zero(reduced_reference-s.eye(2)/2)


def test_native_coherent_family_uses_the_same_observer_map():
    data = h.coherent_sector()
    assert s.simplify(data['g_squared']-(425+100*s.sqrt(2))/1944) == 0
    assert s.simplify(data['gap']**2-(911+100*s.sqrt(2))/486) == 0
    mean = h.ground_correlation(data, 0)['mean']
    assert h._clean(s.trace(data['ground_projector']*data['occupation'])-mean) == 0
    records = h.ground_joint_law(data, s.Rational(1, 8))
    assert h._clean(sum(records.values())-1) == 0


@pytest.mark.parametrize('y,cost', [(s.zeros(3), 1), (s.eye(2), 1),
                                   (s.eye(3), s.I), (s.eye(3), .5),
                                   (s.diag(1, 2, .1), 1)])
def test_nonadmitted_exact_inputs_fail_closed(y, cost):
    with pytest.raises(ValueError):
        h.sector_data(y, cost)


@pytest.mark.parametrize('p,q,time', [((0, 0), (0, 0), 1),
                                     ((s.Rational(1, 10), s.Rational(1, 5)), (1, 1), 1),
                                     ((s.Rational(1, 4), s.Rational(1, 4)), (0, 0), 1),
                                     ((s.Rational(1, 4), s.Rational(1, 4)), (s.Rational(1, 2), s.Rational(1, 2)), 0)])
def test_unresolved_inverse_domains_fail_closed(p, q, time):
    with pytest.raises(ValueError):
        h.gap_interval_from_records(p, q, time)
