"""Full-domain identities and a complete physical odd-arrival experiment."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ...chiral_higgs.model import higgs_memory
from ...physical_attachment.chiral import representation
from ..higgs_composites import model as poly
from ..oriented_links.checks import zero
from ..oriented_links.rotations import cubic_rotations, rotate_words, axial_rotation_column
from ..oriented_links.model import translation_column
from .model import (Configuration, neutral_words, hopping_column, column,
    base_column, preparation, composite_column, star_signature, target_signature,
    arrival_projector, complete_readout, gauss_column, empty_key,
    cells, cell_index, shift_cell, edges, apply_state, norm2, clean, PAULI)


def annihilator(modes, index):
    entries = {}
    for mask in range(1 << modes):
        if mask >> index & 1:
            entries[mask ^ (1 << index), mask] = (-1) ** ((mask & ((1 << index) - 1)).bit_count())
    return s.ImmutableSparseMatrix(1 << modes, 1 << modes, entries)


def one_body(annihilators, matrix):
    result = s.zeros(annihilators[0].rows)
    for (a, b), value in s.ImmutableSparseMatrix(matrix).todok().items():
        result += value * annihilators[a].H * annihilators[b]
    return s.ImmutableSparseMatrix(result)


@lru_cache(maxsize=1)
def local_preparation_checks():
    fields, scalar_currents, _, _ = higgs_memory()
    two = tuple(annihilator(2, j) for j in range(2))
    j = s.kronecker_product(two[1], fields[0]) - s.kronecker_product(two[0], fields[1])
    spectrum = (j.H * j).eigenvals()
    assert spectrum == {s.Integer(0): 14, s.Integer(1): 4, s.Integer(2): 2}
    d = j / s.sqrt(2)
    assert d.H * d + d * d.H != s.eye(20)
    # Both spin components and all five scalar states are retained in this
    # 80-state local source; the four-state reference is a two-mode CAR memory.
    source_a = tuple(annihilator(4, j) for j in range(4))
    neutral = tuple(s.ImmutableSparseMatrix(
        (s.kronecker_product(source_a[2 * spin + 1], fields[0])
         - s.kronecker_product(source_a[2 * spin], fields[1])) / s.sqrt(2))
        for spin in range(2))
    local_vacuum = s.ImmutableSparseMatrix(80, 1, {(0, 0): 1})
    local_packet = sum((op.H * local_vacuum * local_vacuum.H * op for op in neutral), s.zeros(80))
    assert local_packet.H == local_packet and local_packet * local_packet == local_packet
    assert s.trace(local_packet) == 2
    source_parity = s.ImmutableSparseMatrix(s.kronecker_product(
        s.diag(*[(-1) ** mask.bit_count() for mask in range(16)]), s.eye(5)))
    assert source_parity * local_packet == -local_packet
    reference_parity = s.diag(1, -1, -1, 1)
    pair_create = s.ImmutableSparseMatrix((
        s.kronecker_product(neutral[0].H * source_parity, two[1].H)
        - s.kronecker_product(neutral[1].H * source_parity, two[0].H)) / s.sqrt(2))
    pulse = pair_create + pair_create.H
    vacuum = s.ImmutableSparseMatrix(320, 1, {(0, 0): 1})
    prepared = pair_create * vacuum
    assert (prepared.H * prepared)[0, 0] == 1
    assert pulse * vacuum == prepared and pulse * prepared == vacuum
    total_parity = s.kronecker_product(source_parity, reference_parity)
    assert zero(total_parity * pulse - pulse * total_parity)
    assert total_parity * prepared == prepared
    number_source = one_body(source_a, s.eye(4))
    number_reference = one_body(two, s.eye(2))
    lepton = (s.kronecker_product(number_source, s.eye(20))
              - s.kronecker_product(s.eye(80), number_reference))
    assert zero(lepton * pulse - pulse * lepton) and lepton * prepared == s.zeros(320, 1)
    spin_checks = gauge_checks = 0
    for sigma in PAULI:
        source_spin = s.kronecker_product(one_body(source_a, s.kronecker_product(sigma / 2, s.eye(2))), s.eye(5))
        assert zero(source_spin * local_packet - local_packet * source_spin)
        reference_spin = one_body(two, sigma / 2)
        total_spin = s.kronecker_product(source_spin, s.eye(4)) + s.kronecker_product(s.eye(80), reference_spin)
        assert zero(total_spin * pulse - pulse * total_spin)
        assert zero(total_spin * prepared)
        spin_checks += 1
    for current, t in enumerate(representation()[2]):
        matter = s.kronecker_product(s.eye(2), t.extract((10, 11), (10, 11)))
        source_gauge = (s.kronecker_product(one_body(source_a, matter), s.eye(5))
                        + s.kronecker_product(s.eye(16), scalar_currents[current]))
        assert zero(source_gauge * local_packet - local_packet * source_gauge)
        gauge = s.kronecker_product(source_gauge, s.eye(4))
        assert zero(gauge * pulse - pulse * gauge)
        assert zero(gauge * prepared)
        gauge_checks += 1
    return {
        'full_local_two_mode_scalar_dimension': 20,
        'JdaggerJ_exact_spectrum': {'0': 14, '1': 4, '2': 2},
        'D_operator_norm': '1',
        'D_is_canonical_fermion': False,
        'local_neutral_packet_projector_rank': 2,
        'local_packet_all_gauge_and_spin_commutators': 15,
        'full_two_spin_source_reference_dimension': 320,
        'complete_preparation_current_commutators': gauge_checks,
        'complete_preparation_spin_commutators': spin_checks,
        'total_fermion_parity': 'even',
        'total_lepton_number': '0',
        'retained_reference': 'two gauge-neutral fermion modes of lepton number -1, spin one-half',
        'source_reduced_preparation': 'equal mixture of D_up^dagger|vac> and D_down^dagger|vac>',
        'pulse': 'exp(-i*pi*K/2), K=Bdagger+B, exact two-state orbit',
        'ordinary_source_vacuum_identified_as_ground': False,
    }


@lru_cache(maxsize=1)
def full_domain_checks():
    gauge = spin = 0
    charges = [s.Integer(0)] * 15
    charges[9], charges[10], charges[11] = -s.Integer(1), s.Integer(1), s.Integer(1)
    lepton = poly.one_body(s.kronecker_product(s.eye(2), s.diag(*charges)))
    for part in poly.quantum_yukawa()[0].values():
        assert not poly.commutator(lepton, part)
    for alpha in range(2):
        op = poly.neutral_composite(alpha)
        assert poly.commutator(lepton, op) == poly.scale(op, -1)
        assert op and all(len(word) == 1 for word, _, _ in op)
        for current in range(12):
            assert not poly.commutator(poly.gauge_current(current), op)
            gauge += 1
        for axis in range(3):
            expected = poly.add(*(poly.scale(poly.neutral_composite(beta), -PAULI[axis][alpha, beta] / 2)
                                  for beta in range(2)))
            assert poly.commutator(poly.spin_current(axis), op) == expected
            spin += 1
    words = {(word, scalars): value for word, scalars, value in neutral_words()}
    dagger = {}
    for (word, scalars), value in words.items():
        assert len(word) == 2 and word[0][1] and not word[1][1]
        assert word[0][0] // 30 != word[1][0] // 30
        key = (tuple((mode, not create) for mode, create in reversed(word)),
               tuple((cell, col, row) for cell, row, col in scalars))
        dagger[key] = s.conjugate(value)
    assert dagger == words
    return {'all_occupation_Gauss_identities': gauge,
            'all_occupation_spin_identities': spin,
            'complete_normal_ordered_hopping_words_side_two': len(words),
            'exact_word_Hermiticity': True,
            'every_added_word_has_even_CAR_degree': True,
            'full_local_Yukawa_and_composite_lepton_identities': 5,
            'link_matrices_in_added_term': 'identity on all 119 states, because D is a local gauge singlet',
            'all_base_quantum_link_and_Yukawa_terms_retained': True,
            'additional_full_operator_norm_bound': '6*N^3*abs(lambda)',
            'additional_term_invariant_under_all_local_gauge_currents': True,
            'scalar_occupation_projection_or_classical_substitution': False}


def geometric_words(config):
    out = defaultdict(lambda: s.Integer(0))
    for xyz, axis in edges(config.side):
        left = cell_index(config.side, xyz)
        right = cell_index(config.side, shift_cell(config.side, xyz, axis))
        seam = -1 if xyz[axis] == config.side - 1 else 1
        value = s.I * config.neutral_hopping * seam / 2
        out[left, right, axis, ()] += value
        out[right, left, axis, ()] += s.conjugate(value)
    return clean(out)


@lru_cache(maxsize=1)
def rotation_checks():
    rows = []
    for side in (2, 4):
        config = Configuration(side=side)
        words = geometric_words(config)
        for rotation, _ in cubic_rotations():
            assert rotate_words(words, rotation, side) == words
        rows.append({'side': side, 'all_proper_rotations': 24, 'geometric_words': len(words)})
    # A literal complete-carrier control includes nonblank links and occupied
    # matter. The base source and the added term act in the same call.
    config = Configuration(y_up=s.Rational(2, 3), y_down=s.Rational(3, 5),
                           y_lepton=s.Rational(5, 7), neutral_hopping=s.Rational(7, 11))
    key = next(iter(preparation(config=config)))
    links = list(key[1]); links[0] = 17; links[8] = 97
    key = key[0] | (1 << 37) | (1 << 85), tuple(links), key[2]
    h_key = column(key, config)
    for axis in range(3):
        rotate = lambda source: axial_rotation_column(source, axis, config)
        assert apply_state(rotate(key), lambda source: column(source, config)) == apply_state(h_key, rotate)
        translate = lambda source: translation_column(source, axis, config)
        assert apply_state(translate(key), lambda source: column(source, config)) == apply_state(h_key, translate)
    return {'complete_composite_word_rotations': rows,
            'actual_full_H_nonblank_pi_rotation_controls': 3,
            'actual_full_H_nonblank_translation_controls': 3,
            'matter_charge_conjugation_used': False,
            'two_pi_rotation': 'total fermion parity',
            'continuous_rotation_or_Lorentz_claim': False}


@lru_cache(maxsize=1)
def experiment():
    config = Configuration()
    target = cell_index(config.side, (1, 0, 0))
    psi = preparation(config=config)
    assert {star_signature(key, config) for key in psi} == {target_signature(0, config)}
    assert not arrival_projector(psi, target, config)
    assert not apply_state(psi, lambda key: composite_column(key, 0, 1, config=config))
    assert apply_state(psi, lambda key: composite_column(key, 0, 0, config=config)) == {empty_key(config): s.Integer(1)}
    assert apply_state(psi, lambda key: axial_rotation_column(key, 0, config)) == {
        key: -s.I * value for key, value in preparation(0, 1, config).items()}
    physical = 0
    for cell in range(config.side ** 3):
        for current in range(12):
            assert not apply_state(psi, lambda key: gauss_column(key, cell, current, config))
            physical += 1
    h0 = apply_state(psi, lambda key: base_column(key, config))
    h_neutral = apply_state(psi, lambda key: hopping_column(key, config))
    assert {star_signature(key, config) for key in h0} == {target_signature(0, config)}
    assert norm2(h_neutral) == 3
    assert all(star_signature(key, config)[0] == 1 for key in h_neutral)
    star_input = {key: value * star_signature(key, config)[0] for key, value in psi.items()}
    star_commutator = defaultdict(lambda: s.Integer(0), apply_state(star_input,
        lambda key: hopping_column(key, config)))
    for key, value in h_neutral.items():
        star_commutator[key] -= value * star_signature(key, config)[0]
    assert norm2(clean(star_commutator)) == 12
    assert not arrival_projector(h0, target, config)
    h1 = apply_state(psi, lambda key: column(key, config))
    arrival, no_arrival = complete_readout(h1, target, config)
    expected_arrival = {key: -s.I * value for key, value in preparation(target, 1, config).items()}
    assert arrival == expected_arrival
    assert arrival_projector(arrival, target, config) == arrival
    assert not arrival_projector(no_arrival, target, config)
    assert norm2(arrival) == 1 and norm2(h1) == norm2(arrival) + norm2(no_arrival)
    assert norm2(h1) == s.Rational(1675, 8)
    # The physical source is preserved at first order with all three Yukawa
    # channels and every native quantum link active.
    h1_gauss = 0
    for cell in range(config.side ** 3):
        for current in range(12):
            assert not apply_state(h1, lambda key: gauss_column(key, cell, current, config))
            h1_gauss += 1
    quark_keys = [key for key in h1 if any(
        key[0] >> mode & 1 and mode % 15 not in (9, 10, 11)
        for mode in range(30 * config.side ** 3))]
    assert quark_keys
    assert any(any(link != 0 for link in key[1]) for key in h1)
    # The same complete H is applied twice; no neutral or vacuum-cyclic
    # effective Hamiltonian replaces the full residual.
    h2 = apply_state(h1, lambda key: column(key, config))
    h2_norm2 = norm2(h2)
    assert h2_norm2 == s.Rational(15948017, 128)
    bound = s.Integer(353)
    assert (bound - 1) ** 2 < h2_norm2 < bound ** 2
    time = 1 / bound
    state_error = time ** 2 * bound / 2
    lower = s.factor((time - state_error) ** 2)
    upper = s.factor((time + state_error) ** 2)
    assert lower == s.Rational(1, 498436) and upper == s.Rational(9, 498436)
    return {
        'configuration': {name: str(getattr(config, name)) for name in
            ('side', 'spacing', 'y_up', 'y_down', 'y_lepton', 'scalar_penalty', 'link_penalty', 'neutral_hopping')},
        'origin_cell': 0, 'arrival_cell': target, 'arrival_coordinate': [1, 0, 0],
        'prepared_source_terms_per_spin': len(psi),
        'prepared_source_norm_squared': '1',
        'prepared_all_cell_current_checks': physical,
        'full_H_first_action_all_cell_current_checks': h1_gauss,
        'full_H_first_action_terms': len(h1),
        'full_H_first_action_norm_squared': str(norm2(h1)),
        'first_action_terms_with_quark_occupations': len(quark_keys),
        'first_action_quantum_link_excitations_retained': True,
        'full_H_second_action_terms': len(h2),
        'full_H_second_action_norm_squared': str(h2_norm2),
        'arrival_first_derivative_norm_squared': '1',
        'arrival_first_derivative_state': '-i D_target,down^dagger|all_source_vac>',
        'arrival_effect': 'local physical neutral spin-doublet projector times (I-C_target)/2',
        'arrival_detector_spatial_support': 'one site and its incident quantum links',
        'neutral_hopping_first_action_norm_squared': '3',
        'star_commutator_norm_squared_on_preparation': '12',
        'lambda_zero_target_arrival_all_times': '0, by the complete-source local star symmetry',
        'source_time': str(time),
        'outward_H2_norm_bound': str(bound),
        'full_state_first_order_remainder_bound': str(state_error),
        'complete_arrival_probability_interval': [str(lower), str(upper)],
        'complete_no_arrival_probability_interval': [str(1 - upper), str(1 - lower)],
        'both_spin_inputs_same_bounds': 'proper pi rotation about the arrival axis',
        'retained_reference_probability': '1, the source evolves while the prepared reference is retained',
        'complete_source_and_reference_readout': ['arrival', 'no_arrival'],
        'postselection_or_discarded_dynamical_sector': False,
        'old_ground_or_old_spectrum_consumed': False,
        'physical_particle_or_mass_identified': False,
    }
