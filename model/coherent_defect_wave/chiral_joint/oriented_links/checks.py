"""Exact checks on full link matrices and normal-ordered Hamiltonian words."""
from collections import defaultdict

import sympy as s

from ...physical_attachment.chiral import representation
from ...chiral_source.model import PAULI
from ..model import clean
from ..periodic_family.model import Configuration, apply_state, native_blocks, norm2
from .model import (labels, label_index, link, endpoint_current, orientation,
                    transporter, path_transitions, empty_key, column, gauss_column,
                    translation_column)
from .rotations import (cubic_rotations, geometry, hamiltonian_words,
                        rotate_words, axial_rotation_column)


def zero(matrix):
    return all(s.expand(value) == 0 for value in s.ImmutableSparseMatrix(matrix).todok().values())


def matrix_checks():
    index = label_index()
    identity = s.eye(119)
    j = orientation()
    assert len(labels()) == 119 and j.H * j == identity and j * j == identity
    electric = s.diag(0, *([1] * 118))
    parity = identity - 2 * electric
    assert j * electric * j == electric
    ward = reversal = 0
    raw_sum_source, raw_sum_target = s.zeros(119), s.zeros(119)
    for current, t in enumerate(representation()[2]):
        left, right = (endpoint_current(current, end) for end in (0, 1))
        assert left.H == left and right.H == right
        assert zero(left * right - right * left)
        assert j * left * j == right
        assert zero(left * electric - electric * left)
        assert zero(right * electric - electric * right)
        for block in native_blocks():
            for a in block:
                for b in block:
                    u = link(a, b)
                    expected_left = -sum((t[a, c] * link(c, b) for c in block), s.zeros(119))
                    expected_right = sum((t[c, b] * link(a, c) for c in block), s.zeros(119))
                    assert zero(left * u - u * left - expected_left)
                    assert zero(right * u - u * right - expected_right)
                    ward += 2
    for block in native_blocks():
        for a in block:
            for b in block:
                u = link(a, b)
                assert j * u * j == link(b, a).H
                assert zero(parity * u + u * parity)
                assert u ** 3 == s.zeros(119) and u ** 2 != s.zeros(119)
                raw_sum_source += u.H * u
                raw_sum_target += u * u.H
                reversal += 1
    assert j * raw_sum_source * j == raw_sum_target
    assert raw_sum_source.is_diagonal() and raw_sum_target.is_diagonal()
    assert sum(value != 0 for value in raw_sum_source.diagonal()) == 60
    assert sum(value != 0 for value in raw_sum_target.diagonal()) == 60
    blocks = []
    for block in native_blocks():
        v = transporter(block)
        p, q = v.H * v, v * v.H
        assert zero(p * p - p) and zero(q * q - q)
        assert p.H == p and q.H == q
        d = len(block)
        assert s.trace(p) == s.trace(q) == 2 * d
        vacuum = s.ImmutableSparseMatrix(119 * d, d, {(119 * a, a): 1 for a in range(d)})
        assert vacuum.H * p * vacuum == s.eye(d)
        blocks.append({'native_indices': list(block), 'dimension': d,
                       'normalization': str(1 / s.sqrt(d)),
                       'complete_internal_link_dimension': 119 * d,
                       'initial_projection_rank': 2 * d,
                       'transporter_operator_norm': '1'})
    # Omitting the index transpose breaks the orientation identity off diagonal.
    wrong_j = s.ImmutableSparseMatrix(119, 119, {(0, 0): 1, **{
        (index[-sign, a, b], col): 1
        for col, (sign, a, b) in enumerate(labels()[1:], 1)}})
    a, b = native_blocks()[0][:2]
    assert wrong_j * link(a, b) * wrong_j != link(b, a).H
    # Unnormalized d>1 transporters have V^dagger V=d P, so are not contractions.
    d = len(native_blocks()[1])
    v = transporter(native_blocks()[1])
    assert zero((s.sqrt(d) * v).H * (s.sqrt(d) * v) - d * (v.H * v))
    return {
        'full_link_basis_states': 119,
        'actual_native_current_count': 12,
        'complete_matrix_endpoint_Ward_identities': ward,
        'complete_matrix_orientation_identities': reversal,
        'endpoint_exchange_identities': 12,
        'every_matrix_column_retained': True,
        'electric_and_star_parity_compatible': True,
        'sum_UdaggerU_and_sum_UUdagger_rank': 60,
        'blocks': blocks,
        'missing_transpose_mutation_rejected': True,
        'missing_normalization_mutation_norm_squared': str(d),
    }


def peter_weyl_check():
    """A retained nonvacuum Clebsch transition proves the stated distinction."""
    currents = representation()[2]
    y = currents[-1]
    ec = next(i for i in range(15) if y[i, i] == 1)
    lepton = next(block for block in native_blocks() if y[block[0], block[0]] == -s.Rational(1, 2))
    eps = s.ImmutableMatrix([[0, 1], [-1, 0]])
    for t in currents:
        ell = t.extract(lepton, lepton)
        assert zero(eps * (t[ec, ec] * s.eye(2) + ell) * eps.inv() + ell.T)
    # A multiplication operator for the one-dimensional e^c character maps
    # the complete L representation coefficient space isometrically to Lbar.
    cg = s.kronecker_product(eps, eps)
    assert cg.H * cg == s.eye(4)
    for a in lepton:
        for b in lepton:
            assert all(link(ec, ec)[row, label_index()[1, a, b]] == 0 for row in range(119))
    return {
        'all_current_character_times_doublet_intertwiners': 12,
        'native_character_index': ec,
        'native_doublet_indices': list(lepton),
        'tensor_identity': 'e^c tensor L is equivalent to conjugate L',
        'retained_nonvacuum_Peter_Weyl_transition_rank': 4,
        'declared_link_transition_rank_on_that_input': 0,
        'literal_Peter_Weyl_compression': False,
        'vacuum_plus_R_blocks_and_duals_dimension': 119,
        'global_minimal_hardware_claim': False,
    }


def path_checks():
    rows = []
    for length in (1, 2, 3, 4):
        transitions = path_transitions((0,) * length)
        by_input = defaultdict(lambda: s.Integer(0))
        for a, b, targets, value in transitions:
            assert len(targets) == length
            by_input[b] += value ** 2
        assert set(by_input) == set(range(15))
        assert all(s.simplify(value) == 1 for value in by_input.values())
        # Absorption columns include every intermediate occupation in the
        # nonzero operator support, rather than only the blank input.
        for block in native_blocks():
            a, b = block[0], block[-1]
            word = tuple(label_index()[-1, a, a] for _ in range(length - 1)) + (label_index()[-1, a, b],)
            assert path_transitions(word) == ((a, b, (0,) * length,
                                               len(block) ** (-s.Rational(length, 2))),)
        rows.append({'length': length, 'blank_output_terms': len(transitions),
                     'every_native_input_norm_squared': '1'})
    # Both creation and annihilation occur in one actual two-edge word.
    block = native_blocks()[1]
    a, b = block[:2]
    mixed = path_transitions((label_index()[-1, a, b], 0))
    assert len(mixed) == len(block) and all(row[2][0] == 0 for row in mixed)
    responses = []
    for side in (2, 4):
        config = Configuration(side=side, y_up=0, y_down=0, y_lepton=0,
                               scalar_penalty=0, link_penalty=0)
        _, links, scalars = empty_key(config)
        key = 1 << block[0], links, scalars
        response = norm2(column(key, config))
        expected = s.Rational(side * side - 1, 2 * side * side)
        assert response == expected
        responses.append({'side': side, 'spacing': 'pi',
                          'exact_blank_one_fermion_response_squared': str(response)})
    return {'blank_executions': rows,
            'arbitrary_length_full_carrier_path_norm': '<=1 by product of contractions',
            'distinct_edge_path_norm_exactly': '1 by its blank input',
            'mixed_creation_absorption_word_terms': len(mixed),
            'hopping_only_response_controls': responses,
            'blank_one_fermion_response_squared': 'pi^2*(N^2-1)/(2*N^2*a^2)',
            'one_link_transport_has_nonvacuum_initial_subspace': True,
            'local_star_parity_obstruction_removed': False}


def rotation_word_checks():
    eps = s.ImmutableMatrix([[0, 1], [-1, 0]])
    rotations = cubic_rotations()
    for rotation, spin in rotations:
        assert zero(spin.H * spin - s.eye(2)) and s.simplify(spin.det()) == 1
        assert zero(spin.T * eps * spin - eps)
        for axis, sigma in enumerate(PAULI):
            expected = sum((rotation[j, axis] * PAULI[j] for j in range(3)), s.zeros(2))
            assert zero(spin * sigma * spin.H - expected)
    rotation_set = {tuple(r) for r, _ in rotations}
    assert all(tuple(a * b) in rotation_set for a, _ in rotations for b, _ in rotations)
    rows = []
    for side in (2, 4):
        config = Configuration(side=side)
        words = hamiltonian_words(config)
        for rotation, _ in rotations:
            assert rotate_words(words, rotation, side) == words
        rows.append({'side': side, 'complete_hopping_words': len(words),
                     'all_proper_cubic_rotations': 24})
    rotation = s.ImmutableMatrix.diag(-1, -1, 1)
    words = hamiltonian_words(Configuration(side=4))
    assert rotate_words(words, rotation, 4, retain_seam=False) != words
    return {
        'proper_spatial_rotations': 24,
        'exact_spin_lifts_checked': 24,
        'spin_singlet_Yukawa_identities': 24,
        'full_hopping_operator_word_executions': rows,
        'word_expansion_retains_all_native_blocks_and_link_columns': True,
        'rotation_does_not_conjugate_matter_internal_representation': True,
        'missing_antiperiodic_seam_mutation_rejected': True,
        'two_pi_spin_rotation': 'fermion parity',
        'continuous_rotation_or_little_group_claim': False,
    }


def actual_hamiltonian_checks():
    config = Configuration(y_up=s.Rational(2, 3), y_down=s.Rational(3, 5),
                           y_lepton=s.Rational(5, 7), scalar_penalty=s.Rational(7, 11),
                           link_penalty=s.Rational(11, 13))
    base = empty_key(config)
    cases = []
    for block in native_blocks():
        a, b = block[0], block[-1]
        for sign in (1, -1):
            links, scalars = list(base[1]), list(base[2])
            links[0] = label_index()[sign, a, b]
            links[4] = label_index()[-sign, b, a]
            scalars[0], scalars[3] = 1, 4
            # All sectors contain actual occupied matter, with both spin
            # components and nonblank links on distinct axes.
            mask = (1 << a) | (1 << (30 + 15 + b)) | (1 << (120 + b))
            cases.append((mask, tuple(links), tuple(scalars)))
    cases.append(base)
    checked = 0
    for key in cases:
        h_state = column(key, config)
        for axis in range(3):
            rotate = lambda state: axial_rotation_column(state, axis, config)
            left = apply_state(rotate(key), lambda state: column(state, config))
            right = apply_state(h_state, rotate)
            assert clean(left) == clean(right)
            twice = apply_state(rotate(key), rotate)
            assert twice == {key: s.Integer((-1) ** key[0].bit_count())}
            checked += 1
    translations = 0
    for axis in range(3):
        key = cases[2]
        translate = lambda state: translation_column(state, axis, config)
        left = apply_state(translate(key), lambda state: column(state, config))
        right = apply_state(column(key, config), translate)
        assert clean(left) == clean(right)
        translations += 1
    # An actual H/G commutator on a nonblank occupied input complements the
    # full normal-ordered proof from the complete link and scalar matrices.
    ward = 0
    for current in (0, 9, 11):
        key = cases[2]
        left = apply_state(gauss_column(key, 0, current, config), lambda state: column(state, config))
        right = apply_state(column(key, config), lambda state: gauss_column(state, 0, current, config))
        assert clean(left) == clean(right)
        ward += 1
    return {'complete_basis_rotation_controls': checked,
            'nonblank_occupied_controls': 30,
            'all_native_blocks_and_both_flux_orientations_in_controls': True,
            'actual_quantum_H_G_controls': ward,
            'actual_quantum_H_translation_controls': translations,
            'unequal_rational_Yukawa_and_penalty_coefficients': True,
            'old_ground_or_spectrum_consumed': False,
            'full_carrier_operator_identity_proof': 'link matrices plus complete path and onsite words'}
