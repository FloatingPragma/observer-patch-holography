"""Exact controls for event-cost selection and complete quantum-link CAR."""
import pytest
import sympy as s

from model.particle_calibration_selection import families_chirality as f


def test_full_carrier_event_cost_is_laplacian_without_dimension_window():
    d = f.event_cost_data()
    assert d['spectrum'] == {s.Integer(0): 2, s.Integer(2): 20,
                            s.Rational(5, 2): 8,
                            (5-s.sqrt(5))/3: 3, (5+s.sqrt(5))/3: 3}
    assert s.trace(d['projector']) == 3
    assert f._zero(d['cost']*d['invariant'])
    assert not f._zero(d['cost']*d['projector'])
    assert len(d['rotations']) == 12
    assert f._zero(d['all_rotation_cost']-2*(s.eye(36)-d['invariant']))
    assert s.trace(d['all_rotation_cost']) == 68


def test_family_cost_rank_requires_selected_controller_and_matter_attachment():
    row = f.family_certificate()
    assert row['smallest_positive_cost_rank'] == 3
    assert row['alternative_covariant_cooling_rank'] == 20
    assert row['invariant_zero_cost_rank'] == 2
    assert not row['dimension_window_used_in_this_derivation']
    assert not row['passive_corrected_C_E_selects_the_cooling_service']
    assert not row['three_physical_generations_selected']


def test_external_copy_factor_is_invisible_to_source_words():
    values = [f.replica_trace_control(n) for n in (1, 2, 3, 4, 7)]
    assert len(set(values)) == 1
    assert values[0] != 0
    for invalid in (0, 9, True, 2.0):
        with pytest.raises(ValueError):
            f.replica_trace_control(invalid)


def test_canonical_dressing_retains_full_link_and_all_occupations():
    d = f.dressed_data()
    assert d['unitary'].shape == (240, 240)
    assert f._zero(d['unitary'].H*d['unitary']-s.eye(240))
    assert f._zero(d['unitary']*d['total_number']-d['total_number']*d['unitary'])
    # Vacuum and double occupation each retain all60 link states.
    for mask in (0, 3):
        indices = list(range(60*mask, 60*(mask+1)))
        assert d['unitary'].extract(indices, indices) == s.eye(60)
    for i in range(2):
        for j in range(2):
            b, c = d['dressed'][i], d['dressed'][j]
            assert f._zero(b*c+c*b)
            assert f._zero(b*c.H+c.H*b-(s.eye(240) if i == j else s.zeros(240)))


def test_all_native_endpoint_currents_and_wrong_orientation_control():
    d = f.dressed_data()
    assert len(d['gauss']) == 24
    for generator in d['gauss']:
        assert f._zero(generator*d['unitary']-d['unitary']*generator)
    internal = d['internal_index']
    g = next(g for g in f.gauge.currents() if g[internal, internal] != 0)
    n = f.creator(2, 0)*f.creator(2, 0).H
    wrong = s.kronecker_product(g[internal, internal]*n, s.eye(60))-s.kronecker_product(
        s.eye(4), f.gauge.link_generator(g, 0))
    assert not f._zero(wrong*d['unitary']-d['unitary']*wrong)


def test_rank_one_sign_projector_requires_explicit_dark_completion():
    d = f.dressed_data()
    assert s.trace(d['old_projector']) == 1
    assert s.trace(d['one_particle_projector']) == 60
    assert s.trace(d['dark_completion']) == 59
    assert f._zero(d['dark_completion']**2-d['dark_completion'])
    assert f._zero(d['old_projector']*d['dark_completion'])
    assert f._zero(d['kernel']**2-s.eye(120))
    # A unitary never changes the rank of a fixed one-mode x full-link
    # projector. Treating the rank-one kernel as that object must fail.
    assert s.trace(d['old_projector']) % 60 != 0


@pytest.mark.parametrize('strength', [2, 4, 8])
def test_complete_penalty_spectrum_and_controlled_complement(strength):
    d = f.dressed_data()
    h = f._clean(d['parent']+strength*d['heavy'])
    report = f.penalty_certificate(strength)
    low = s.sympify(report['mixed_low_energy'])
    high = s.sympify(report['mixed_high_energy'])
    # Check the exact characteristic factors through their matrix minimal
    # polynomial, with trace and trace-square fixing the listed multiplicities.
    identity = s.eye(240)
    product = f._clean((h-identity)*(h-(1+strength)*identity)
                      *(h-low*identity)*(h-high*identity))
    assert f._zero(product)
    assert s.simplify(s.trace(h)-(119+119*(1+strength)+low+high)) == 0
    assert s.simplify(s.trace(h*h)-(119+119*(1+strength)**2+low**2+high**2)) == 0
    active = {60+d['active_link_index'], 120}
    inactive = [i for i in range(240) if i not in active]
    outside = h.extract(inactive, inactive)
    assert outside == s.diag(*outside.diagonal())
    assert sum(value == 1 for value in outside.diagonal()) == 119
    assert sum(value == 1+strength for value in outside.diagonal()) == 119
    assert f._zero(h.extract(inactive, sorted(active)))
    assert f._zero(h.extract(sorted(active), inactive))
    assert s.simplify(s.sympify(report['schur_bound'])-(1-low)).is_positive
    for g in d['gauss']:
        assert f._zero(g*h-h*g)


def test_unitary_dressing_does_not_remove_parent_spectrum():
    d = f.dressed_data()
    assert f._zero(d['unitary']*d['parent']*d['unitary'].H-d['parent'])
    report = f.chirality_certificate()
    assert report['balanced_sign_control_dark_zero_columns'] == 118
    assert not report['balanced_sign_control_is_a_physical_kinetic_regulator']
    assert not report['interacting_mirror_decoupling_or_Lorentz_limit_proved']
    assert not report['native_state_or_physical_clock_selected']
    assert not report['source_pulse_tape_executed_for_this_dressing']


@pytest.mark.parametrize('strength', [0, 1, -1, 1.5, s.sqrt(2)])
def test_penalty_bound_domain_is_explicit(strength):
    with pytest.raises(ValueError):
        f.penalty_certificate(strength)
