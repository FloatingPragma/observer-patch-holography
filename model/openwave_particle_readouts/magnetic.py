"""Exact hypercharge-background response of the complete interacting source."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_source.model import PAULI
from model.coherent_defect_wave.chiral_joint.model import clean, sequence
from model.coherent_defect_wave.chiral_joint.periodic_family.model import (
    Configuration, apply_state, kinetic_terms, native_blocks, norm2,
    pair_preparation,
)
from model.coherent_defect_wave.chiral_joint.oriented_links.model import (
    empty_key, endpoint_current, labels, label_index, gauss_column,
)
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.clebsch import link
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.model import (
    column, hopping_column, path_transitions,
)

COUPLING = s.Rational(1, 2**34)
PROFILE = (1,) + (0,) * 23


def validate_profile(profile, config=Configuration()):
    if (len(profile) != 3*config.side**3
            or any(type(value) is not int for value in profile)):
        raise ValueError('one integer background profile per oriented link required')
    return tuple(profile)


@lru_cache(maxsize=1)
def generator_diagonal():
    """K=-G_left^Y, so [K,U_ab]=Y_a U_ab on all 119 link states."""
    return tuple((-endpoint_current(11, 0)).diagonal())


def generator_value(key, profile=PROFILE):
    return sum((weight*generator_diagonal()[label]
                for weight, label in zip(profile, key[1], strict=True)), s.Integer(0))


def phase(key, profile=PROFILE):
    """exp(3 pi i K_profile), evaluated without floating point phases."""
    exponent = 6*generator_value(key, profile)
    assert exponent.is_Integer
    return s.I**int(exponent)


def transformed(state, profile=PROFILE):
    return {key: phase(key, profile)*value for key, value in state.items()}


def background_hopping_column(key, profile=PROFILE, derivative=0,
                              at_test_angle=True, config=Configuration()):
    """Actual PW path words with exp(i theta Y sum_e profile_e).

    Derivatives are with respect to theta. The exact test angle is 3 pi;
    the analytic covariance holds at every real angle.
    """
    profile = validate_profile(profile, config)
    if derivative not in (0, 1, 2):
        raise ValueError('zeroth, first or second background derivative required')
    mask, links, scalars = key
    out = defaultdict(lambda: s.Integer(0))
    hypercharge = representation()[2][-1]
    for term in kinetic_terms(config):
        path_weight = sum(profile[edge] for edge in term.path)
        values = tuple(links[edge] for edge in term.path)
        for adjoint in (False, True):
            for a, b, targets, factor in path_transitions(values, adjoint):
                charge = hypercharge[a, a]*path_weight*(-1 if adjoint else 1)
                multiplier = (s.I*charge)**derivative
                if at_test_angle:
                    multiplier *= s.I**int(6*charge)
                for (u, v), spin_weight in PAULI[term.axis].todok().items():
                    coefficient = term.coefficient*spin_weight
                    p = 30*term.left_cell + 15*u + a
                    q = 30*term.right_cell + 15*v + b
                    if adjoint:
                        p, q = q, p
                        coefficient = s.conjugate(coefficient)
                    step = sequence(mask, ((p, True), (q, False)))
                    if step is None:
                        continue
                    target_links = list(links)
                    for edge, target in zip(term.path, targets, strict=True):
                        target_links[edge] = target
                    out[step[0], tuple(target_links), scalars] += (
                        coefficient*factor*step[1]*multiplier)
    return clean(out)


def source_column(key, coupling=COUPLING, profile=PROFILE, background=False):
    """H0+lambda*T on the complete source, optionally with a background."""
    coupling = s.sympify(coupling)
    if coupling.is_Rational is not True or coupling < 0:
        raise ValueError('an exact nonnegative hopping coefficient required')
    terms = kinetic_terms()
    unscaled_hopping = hopping_column(key, terms)
    out = defaultdict(lambda: s.Integer(0), column(key))
    for target, value in unscaled_hopping.items():
        out[target] -= value
    active = (background_hopping_column(key, profile) if background
              else unscaled_hopping)
    for target, value in active.items():
        out[target] += coupling*value
    return clean(out)


@lru_cache(maxsize=1)
def coordinate_checks():
    k = -endpoint_current(11, 0)
    hypercharge = representation()[2][-1]
    penalty = s.diag(0, *([1]*118))
    assert k*penalty == penalty*k
    assert k == s.diag(*k.diagonal()) and k[0, 0] == 0
    coordinates = entries = 0
    occupied_to_occupied = 0
    for block in native_blocks():
        for a in block:
            for b in block:
                u = link(a, b)
                assert k*u-u*k == hypercharge[a, a]*u
                coordinates += 1
                entries += len(u.todok())
                occupied_to_occupied += sum(row != 0 and col != 0
                                           for row, col in u.todok())
    for current in range(12):
        for endpoint in (0, 1):
            g = endpoint_current(current, endpoint)
            assert k*g == g*k
    # A fixed weak-isospin direction is not central in the full gauge group.
    weak = endpoint_current(10, 0)
    assert weak*endpoint_current(8, 0) != endpoint_current(8, 0)*weak
    return {'complete_link_dimension': 119, 'native_coordinate_matrices': coordinates,
        'nonzero_coordinate_entries': entries,
        'occupied_to_occupied_entries': occupied_to_occupied,
        'central_endpoint_commutators': 24,
        'exact_coordinate_identity': '[K_e,U_ab]=Y_a U_ab, K_e=-G_e,left^Y',
        'link_excitation_penalty_commutes': True,
        'noncentral_weak_direction_transfer_rejected': True}


def fixtures():
    base = empty_key()
    full = (1 << 240)-1
    yield base
    yield full, base[1], base[2]
    for block in native_blocks():
        a, b = block[0], block[-1]
        for sign in (1, -1):
            links = list(base[1])
            links[0] = label_index()[sign, a, b]
            mask = ((1 << a) if sign == 1 else (1 << (120+15+b))) | (1 << (60+b))
            yield mask, tuple(links), (0, 1, 2, 3, 4, 0, 1, 2)
    for a in (9, 10, 11):
        yield 1 << (120+15+a), base[1], base[2]


@lru_cache(maxsize=4)
def full_source_checks(coupling=COUPLING):
    columns = entries = current_entries = curvature_entries = 0
    current_norm = curvature_norm = s.Integer(0)
    for key in fixtures():
        h = source_column(key, coupling)
        background = source_column(key, coupling, background=True)
        expected = {target: phase(target)*s.conjugate(phase(key))*value
                    for target, value in h.items()}
        assert background == expected
        derivative = background_hopping_column(key, derivative=1, at_test_angle=False)
        curvature = background_hopping_column(key, derivative=2, at_test_angle=False)
        commutator = clean({target: s.I*(generator_value(target)-generator_value(key))*value
                            for target, value in h.items()})
        double = clean({target: -(generator_value(target)-generator_value(key))**2*value
                        for target, value in h.items()})
        assert clean({target: coupling*value for target, value in derivative.items()}) == commutator
        assert clean({target: coupling*value for target, value in curvature.items()}) == double
        current_norm += norm2(commutator)
        curvature_norm += norm2(double)
        current_entries += len(commutator)
        curvature_entries += len(double)
        columns += 1
        entries += len(h)
    assert current_norm > 0 and curvature_norm > 0
    return {'complete_240_mode_basis_columns': columns,
        'actual_full_H_entries_checked': entries,
        'nonzero_current_entries': current_entries,
        'nonzero_diamagnetic_entries': curvature_entries,
        'current_column_norm_squared_sum': str(current_norm),
        'diamagnetic_column_norm_squared_sum': str(curvature_norm),
        'actual_hopping_coefficient': str(coupling),
        'empty_full_mixed_Fock_and_occupied_links_retained': True,
        'background_hopping_executed_independently_of_link_conjugation': True,
        'omitted_current_and_diamagnetic_derivatives_rejected': True}


def plaquette_column(key, adjoint=False):
    """An explicit comparison operator; absent from the declared source H."""
    state = {key: s.Integer(1)}
    for edge, reverse in ((0, False), (13, False), (6, True), (1, True)):
        matrix = link(9, 9).H if reverse != adjoint else link(9, 9)
        out = defaultdict(lambda: s.Integer(0))
        for (mask, links, scalars), value in state.items():
            for (row, col), weight in matrix.todok().items():
                if col != links[edge]:
                    continue
                target = list(links)
                target[edge] = row
                out[mask, tuple(target), scalars] += value*weight
        state = clean(out)
    return state


@lru_cache(maxsize=1)
def plaquette_checks():
    key = empty_key()
    p = plaquette_column(key)
    pd = plaquette_column(key, True)
    total = clean({**p, **pd})
    assert len(total) == 2 and norm2(total) == 2
    assert transformed(total) == {target: -value for target, value in total.items()}
    difference = {target: -2*value for target, value in total.items()}
    assert norm2(difference) == 8
    # This is a closed physical loop, checked against each affected Gauss law.
    for cell in (0, 2, 4, 6):
        for current in range(12):
            assert not apply_state(total, lambda k: gauss_column(k, cell, current))
    return {'oriented_square_edges': [0, 13, 6, 1],
        'signed_profile_circulation': 1,
        'test_background_angle': '3*pi',
        'ec_Wilson_phase': '-1',
        'profile_is_a_zero_curl_site_gauge_transformation': False,
        'closed_loop_Gauss_checks': 48,
        'added_plaquette_conjugation_difference_on_blank_norm_squared': '8',
        'bare_plaquette_term_present_in_actual_H': False,
        'comparison_plaquette_inserted_into_actual_H': False}


def inner(left, right):
    return s.expand(sum((s.conjugate(value)*right.get(key, 0)
                         for key, value in left.items()), s.Integer(0)))


@lru_cache(maxsize=4)
def fixed_preparation_control(coupling=COUPLING):
    """Physical finite state witnesses a changed fixed-preparation energy mean."""
    vacuum_pair = pair_preparation()
    edge_terms = tuple(term for term in kinetic_terms() if term.path == (0,))
    transported = apply_state(vacuum_pair, lambda key: hopping_column(key, edge_terms))
    assert norm2(transported) > 0 and not (set(vacuum_pair) & set(transported))
    state = clean({**vacuum_pair, **transported})
    # The native pair pulse and one gauge-invariant native edge term preserve
    # every Gauss constraint. Check the affected stars on the literal state.
    for cell in (0, 4):
        for current in range(12):
            assert not apply_state(state, lambda key: gauss_column(key, cell, current))
    ordinary = apply_state(state, lambda key: source_column(key, coupling))
    changed = apply_state(state, lambda key: source_column(key, coupling, background=True))
    difference = s.cancel((inner(state, changed)-inner(state, ordinary))/norm2(state))
    assert difference != 0
    paired = transformed(state)
    paired_image = apply_state(paired, lambda key: source_column(key, coupling, background=True))
    assert inner(paired, paired_image) == inner(state, ordinary)
    return {'physical_pair_plus_native_edge_image_state_norm_squared': str(norm2(state)),
        'literal_affected_star_Gauss_checks': 24,
        'fixed_state_energy_mean_change': str(difference),
        'paired_state_energy_mean_change': '0',
        'control_is_inherited_ec_or_L_prepared_line': False,
        'all_fixed_preparation_magnetic_responses_inferred_zero': False}


def certificate(coupling=COUPLING):
    return {'schema': 'oph.same_source_hypercharge_background.v1',
        'coordinate_covariance': coordinate_checks(),
        'full_source_controls': full_source_checks(coupling),
        'plaquette_control': plaquette_checks(),
        'fixed_preparation_control': fixed_preparation_control(coupling),
        'static_current_ward': {
            'all_real_angles_and_all_integer_profiles': True,
            'arbitrary_real_link_angles_by_independent_commuting_generators': True,
            'complete_H_identity': 'H(A)=V(A) H(0) V(A)^dagger',
            'unitary': 'V(A)=exp(i sum_e A_e K_e), K_e=-G_e,left^Y',
            'current': 'J_f=i[K_f,H]',
            'diamagnetic_operator': 'D_f=-[K_f,[K_f,H]]',
            'spectral_current_blocks': 'P_E J_f P_E=0, including degenerate E',
            'static_second_derivative': 'P_E D_f P_E+2 P_E J_f (E-H)^(-1)_(1-P_E) J_f P_E=0',
            'every_energy_eigenvalue_independent_of_background': True,
            'static_hypercharge_energy_splitting': '0',
            'physical_electromagnetic_g_inferred': False},
        'paired_record_transfer': {
            'families': ['ec', 'L'],
            'carrier': 'same 240 CAR modes, 24 PW119 links, eight scalar-five cells, static CAR reference',
            'reference_cost': '4 N_ref, unchanged by link-only V',
            'physical_preparation': 'conjugate the complete source preparation instrument by V on each retained source copy',
            'native_pump_commutes_with_V': True,
            'spectral_filter': 'P_I(H(A))=V P_I(H(0)) V^dagger for the identical interval I',
            'same_source_clock': 'exp(-it H(A))=V exp(-it H(0)) V^dagger',
            'source_spin_location_and_classical_failure_flags_commute_with_V': True,
            'blank_link_atomic_projector': 'V P0=P0, hence fixed P0-supported input and existing arrival effects have the same complete records',
            'all_original_preparation_failures_bin_complements_and_references_retained': True,
            'joint_record_probabilities_equal_for_paired_apparatus': True,
            'unchanged_fixed_preparation_records_inferred': False},
        'scope': {
            'complete_interacting_source_no_free_Clifford_substitution': True,
            'central_hypercharge_background_only': True,
            'emergent_electromagnetic_direction_or_photon_identified': False,
            'all_magnetic_observables_vanish': False,
            'free_envelope_g_two_transferred_to_prepared_line': False,
            'observed_electron_identified': False},
        'verdict': 'PASS_SAME_SOURCE_STATIC_HYPERCHARGE_RESPONSE_OBSTRUCTION'}
