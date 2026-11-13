"""Native residue, complete Dyson and retained energy-disturbance controls."""
import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from .motion import (atom_first_hop, detuned_arrival_bounds, ec_flux_attachment,
                     family_certificate, sampling_budget, spin_effects, validate_ec_pole_source)


def physical_matrix(matrix, gram):
    change = np.linalg.cholesky(np.array(gram, dtype=complex)).conj().T
    result = change@np.array(matrix, dtype=complex)@np.linalg.inv(change)
    assert np.linalg.norm(result-result.conj().T) < 1e-11
    return result, change


@pytest.mark.parametrize('family', ['ec', 'L'])
def test_atom_residue_against_complete_native_target_diagonalization(family):
    from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.model import transfer
    from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source

    certified = atom_first_hop(family)
    ground = source(1)
    h, change = physical_matrix(ground['H'], ground['gram'])
    _, vectors = np.linalg.eigh(h)
    source_vector = vectors[:, 0]
    atom = float(sum(certified['atom'])/2)
    components = []
    for spin in range(2):
        block = transfer(family, spin)
        target, target_change = physical_matrix(block['H'], block['gram'])
        energies, basis = np.linalg.eigh(target)
        image = target_change@np.array(block['map'], dtype=complex)@np.linalg.solve(change, source_vector)
        selected = abs(energies-atom) < 1e-8
        components.append(float(np.linalg.norm(basis[:, selected].conj().T@image)**2))
    weight = sum(components)/2
    vector = components[0]-components[1]
    assert abs(weight-float(sum(certified['weight'])/2)) < 1e-11
    assert abs(vector-float(sum(certified['spin_residue'])/2)) < 1e-11
    assert abs(vector/weight-float(certified['spin_to_weight_ratio'])) < 1e-11
    assert certified['moment'][0] > 0


def test_native_three_axis_contractions_and_unfiltered_limit():
    from model.coherent_defect_wave.chiral_joint.interacting_charged_transport.algebra import spin_effects as broad_effects
    (w, c), selected = spin_effects()
    b, broad = broad_effects()
    # Selecting all atomic target energies gives w=1/2 and c=b at both sites.
    for exact, unfiltered in zip(selected, broad, strict=True):
        assert (exact.subs(w, s.Rational(1, 2))-unfiltered.subs(b, c)).applyfunc(s.expand) == s.zeros(4)
    # The L lowest atom is spin zero. Its route effect is scalar even before
    # tracing over the entangled eight-spin cube reference.
    for exact in selected:
        assert (exact.subs(c, 2*w)-w**3*s.eye(4)/4).applyfunc(s.expand) == s.zeros(4)
        assert (exact-exact.H).applyfunc(s.expand) == s.zeros(4)


def test_shared_ec_atom_rejects_a_different_native_weight(monkeypatch):
    from . import motion
    from model.coherent_defect_wave.chiral_joint.interacting_charged_response.local import positive_response
    actual = atom_first_hop('ec')
    local = positive_response('ec')
    assert actual['shared_ec_atomic_source']['complete_native_pole_and_weight_regenerated_and_equal']
    assert not actual['shared_ec_atomic_source']['preparation_or_timeline_probability_bound_transferred']
    changed = dict(motion.energy_resolved.native_pole_residue())
    changed['pump_pole_weight_interval'] = ['1/4', '1/4']
    monkeypatch.setattr(motion.energy_resolved, 'native_pole_residue', lambda: changed)
    with pytest.raises(ArithmeticError, match='shared spectral source'):
        validate_ec_pole_source(actual['selected_factor'], actual['weight'], local['blocks'])


def test_same_native_flux_refines_every_time_and_retains_all_failures(monkeypatch):
    from . import motion
    attached = ec_flux_attachment()
    assert attached['regenerated_route_cut_identities'] == 1536
    assert attached['regenerated_complete_endpoint_Ward_identities'] == 1416
    assert attached['regenerated_additive_factor_basis_Gauss_comparisons'] == 24840
    assert not attached['other_preparation_arrival_probability_transferred']
    assert 'invariant one-particle ec reference sector' in attached['all_state_domain']
    assert not attached['arrival_effect_restricts_the_complete_reference_carrier']
    assert not attached['fixed_ec_flux_atom_asserted_in_other_reference_sectors']
    value = family_certificate('ec')
    for row in value['clock_records']:
        assert row['joint_native_charge_flux_six_route_probability_from_original_input'] == row['six_route_probability_from_original_input']
    for row in value['complete_records']['records']:
        if row['record'] == 'line_prepared_charged_arrival':
            atom = row['joint_native_charge_flux_atom']
            assert atom['arrival_source_matter_plus_scalar_Y_and_Q'] == '-1'
            assert list(map(s.Rational, atom['three_half_cube_fluxes'])).count(-1) == 1
            assert row['joint_atom_probability_equals_actual_arrival_probability']
        else:
            assert 'all central-Y cut-flux and charge outcomes retained' in row['joint_native_flux_refinement']
    full_flux = motion.native_flux.cut_flux_six
    monkeypatch.setattr(motion.native_flux, 'cut_flux_six', lambda mask, links: -full_flux(mask, links))
    with pytest.raises(ArithmeticError, match='full native cut flux'):
        ec_flux_attachment.__wrapped__()


def test_physical_L_arrival_rejects_unrestricted_ec_flux_atom():
    from collections import defaultdict
    from .motion import native_flux
    from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES

    indices = FAMILIES['L']
    dimension = len(indices)
    gauges = native_flux.representation()[2]
    reference_y = gauges[11][indices[0], indices[0]]
    assert reference_y == -s.Rational(1, 2)
    # This conjugate arrival irrep is present in the actual selected L atom,
    # and the complete native first-hop image has strictly positive norm.
    assert atom_first_hop('L')['moment'][0] > 0
    label_index = {label: index for index, label in enumerate(native_flux.labels())}
    def factor_action(vector, factor, operator):
        output = defaultdict(lambda: s.Integer(0))
        for coordinate, amplitude in vector.items():
            for (row, column), value in operator.todok().items():
                if column == coordinate[factor]:
                    target = list(coordinate)
                    target[factor] = row
                    output[tuple(target)] += value*amplitude
        return {key: value for key, value in output.items() if value}

    supported = {}
    for origin_endpoint in (0, 1):
        entries = {}
        for a, native_a in enumerate(indices):
            for b, native_b in enumerate(indices):
                label = ((1, native_a, native_b) if origin_endpoint == 0
                         else (-1, native_b, native_a))
                link_index = label_index[label]
                entries[a, link_index, b] = s.Rational(1, dimension)
        assert sum(value**2 for value in entries.values()) == 1
        supported[origin_endpoint] = {link for _, link, _ in entries}
        # Exact singlets of all twelve native endpoint currents, with the
        # actual L reference and its conjugate arrival-source representation.
        for current, gauge in enumerate(gauges):
            reference = s.ImmutableSparseMatrix(gauge.extract(indices, indices))
            reference_image = factor_action(entries, 0, reference)
            origin_image = factor_action(entries, 1, native_flux.current(current, origin_endpoint))
            assert origin_image == {key: -value for key, value in reference_image.items()}
            arrival_image = factor_action(entries, 1, native_flux.current(current, 1-origin_endpoint))
            assert arrival_image == factor_action(entries, 2, reference.T)
    terms = {term.path[0]: term for term in native_flux.kinetic_terms(native_flux.Configuration())}
    edges = {edge: (left, right) for edge, left, right in native_flux.graph()[0]}
    for route in native_flux.routes():
        edge, arrival = route['edge'], route['arrival']
        origin_endpoint = 0 if edges[edge][0] == 0 else 1
        key = (1 << (30*arrival+indices[0]), (0,)*24, (0,)*8)
        # Regenerate the actual PW119 refilling column, including both
        # periodic orientations, rather than choosing a flux label by hand.
        image = native_flux.pw.hopping_column(key, (terms[edge],))
        assert sum(s.conjugate(amplitude)*amplitude for amplitude in image.values()) == s.Rational(1, 16)
        for (_, links, _), amplitude in image.items():
            assert amplitude != 0 and links[edge] in supported[origin_endpoint]
            for mask in range(256):
                ec_flux = (mask & 1)-((mask >> arrival) & 1)
                assert native_flux.cut_flux_six(mask, links) == 6*reference_y*ec_flux
            assert native_flux.cut_flux_six(1 << arrival, links) == 3
            assert native_flux.cut_flux_six(1 << arrival, links) != -6


@pytest.mark.parametrize('seed', [11, 23, 37])
def test_complete_dyson_bound_on_independent_retained_carrier(seed):
    rng = np.random.default_rng(seed)
    h0 = np.diag([0., 1., 2., 4., 7.])
    raw = rng.normal(size=(5, 5))+1j*rng.normal(size=(5, 5))
    hopping = (raw+raw.conj().T)/2
    hopping /= np.linalg.norm(hopping, 2)
    coupling = s.Rational(1, 10000)
    phi = np.array([.6, 0, 0, 0, 0], dtype=complex)
    # Both retained initial components belong to S0; the arrival is coordinate1.
    epsilon = s.Rational(1, 100000)
    psi = phi+np.array([0, 0, float(epsilon), 0, 0], dtype=complex)
    moment = abs(hopping[1, 0]*phi[0])**2
    enclosed = (s.Rational(str(moment))-s.Rational(1, 10**15),
                s.Rational(str(moment))+s.Rational(1, 10**15))
    for time in (0, 1, 3, 6):
        result = detuned_arrival_bounds(enclosed, coupling, time, epsilon, 1, 1)
        amplitude = (expm(-1j*time*(h0+float(coupling)*hopping))@psi)[1]
        lo, hi = map(float, result['probability_interval'])
        assert lo-1e-20 <= abs(amplitude)**2 <= hi+1e-20
        if time:
            # Dropping the actual flux detuning gives the wrong first-order law.
            resonant = float(coupling)**2*time*time*moment
            if time == 6:
                assert abs(amplitude)**2 < resonant/20


def test_position_projection_retains_probability_and_exposes_energy_disturbance():
    theta = .02
    chi = np.array([np.cos(theta), np.sin(theta)])
    spectral = np.outer(chi, chi)
    s0 = np.diag([1., 0.])
    phi = np.array([1., 0.])
    epsilon = np.linalg.norm(chi-phi)
    localized = s0@chi
    success = np.linalg.norm(localized)**2
    lost = np.linalg.norm((np.eye(2)-spectral)@localized)**2/success
    assert 0 < lost <= epsilon**2/success
    assert np.linalg.norm((np.eye(2)-s0)@chi)**2 <= epsilon**2
    for time in (1, 3, 6):
        evolved = expm(-1j*time*(spectral+3*(np.eye(2)-spectral)))@localized
        assert abs(np.linalg.norm((np.eye(2)-spectral)@evolved)**2/success-lost) < 1e-14


@pytest.mark.parametrize('family', ['ec', 'L'])
def test_actual_prepared_records_certify_excursion_and_return(family):
    value = family_certificate(family)
    by_time = {int(row['source_time']): row for row in value['clock_records']}
    interval = lambda time: tuple(map(s.Rational, by_time[time]['six_route_probability_from_original_input']))
    assert interval(0) == (0, 0)
    assert interval(1)[0] > 0
    assert interval(3)[0] > interval(6)[1] > 0
    assert s.Rational(value['energy_bin_loss_after_position_conditioned_on_position_success_upper']) < s.Rational(1, 10**6)
    assert not value['exact_energy_bin_support_after_position_asserted']
    resources = value['independent_complete_preparation_sampling']
    assert int(resources['time_three_detection_trials']) > 10**21
    assert int(resources['comparison_total_original_input_trials']) == 2*int(resources['time_three_versus_six_trials_per_setting'])
    assert resources['independent_source_apparatus_repetitions_explicitly_required']
    assert not resources['physical_preparation_rate_or_feasible_runtime_inferred']
    records = value['complete_records']['records']
    assert len(records) == value['complete_records']['coarse_record_count'] == 11
    assert len([row for row in records if row['record'] == 'line_prepared_charged_arrival']) == 6
    q, f, outside, line, localized = s.symbols('q f o p a')
    arrivals = s.symbols('r:6')
    complete = (1-q)+q*f+q*outside+q*(line-localized)+q*sum(arrivals)+q*(localized-sum(arrivals))
    assert s.expand(complete.subs(outside, 1-f-line)) == 1
    assert s.expand((complete-q*(line-localized)).subs(outside, 1-f-line)) != 1


def test_detuned_bounds_reject_invalid_source_and_large_error_loses_positive_claim():
    moment = (s.Rational(1, 100),)*2
    assert detuned_arrival_bounds(moment, s.Rational(1, 1000), 3, 1)['probability_interval'][0] == 0
    for args in ((moment, 0, 1, 0), (moment, .01, 1, 0),
                 (moment, 1, -1, 0), (moment, 1, 7, 0),
                 (moment, 1, 1, -1), ((0, 1), 1, 1, 0)):
        with pytest.raises(ValueError):
            detuned_arrival_bounds(*args)


def test_exact_sampling_budget_includes_original_input_failures():
    early = (s.Rational(1, 10**6), s.Rational(2, 10**6))
    late = (s.Integer(0), s.Rational(1, 10**7))
    result = sampling_budget(early, late)
    count = result['detection_trials']
    assert count*early[0] >= 3 > (count-1)*early[0]
    assert result['exp_three_partial_sum_lower'] > 20
    comparison = result['comparison_trials_per_setting']
    contrast = early[0]-late[1]
    assert (early[1]+late[1])/(comparison*contrast**2) <= s.Rational(1, 20)
    assert (early[1]+late[1])/((comparison-1)*contrast**2) > s.Rational(1, 20)
    # Renormalizing by successful preparation would change the required count.
    conditional = sampling_budget(tuple(2*p for p in early), tuple(2*p for p in late))
    assert conditional['detection_trials'] < count
    for left, right in (((0, 1), (0, 0)), ((s.Rational(1, 4), s.Rational(1, 2)), (0, s.Rational(1, 3))),
                        ((.1, .2), (0, .01))):
        with pytest.raises(ValueError):
            sampling_budget(left, right)
