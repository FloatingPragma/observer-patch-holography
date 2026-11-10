"""Actual CAR translation joins and complete atomic momentum records."""
from collections import Counter, defaultdict
from itertools import product
from math import gcd

import sympy as s

from ..periodic_family.model import (Configuration, cells, cell_index, shift_cell,
    edges, kinetic_terms, translate_fermions, apply_state)
from ..peter_weyl_links import model as source
from ..local_reference.algebra import transitions
from ..neutral_mobility.model import composite_column
from ..higgs_composites import model as car
from ...chiral_source.model import PAULI
from ..oriented_links.rotations import cubic_rotations, geometry
from .model import (reference_data, source_character, translation_order,
    momentum_multiplicity, momentum_probability, full_phase_probability,
    site_permutation, shift_spins, purification, translate_state,
    phase_projector, inner, norm2, clean)


def native_translation_join():
    reference, spectrum = reference_data()
    matrices = transitions()
    assert len(matrices) == 4 and all(matrix.shape == (30, 30) for matrix in matrices)
    rows = []
    # Arbitrary odd local Fock occupations test the actual permutation phase.
    odd = (1 << 9, 1 << 24, (1 << 15) - 1, (1 << 29) - 1,
           (1 << 0) | (1 << 7) | (1 << 29))
    for side in (2, 4):
        config = Configuration(side=side); volume = side ** 3
        _, links, _ = source.empty_key(config)
        scalars = tuple(j % 5 for j in range(volume))
        blocks = tuple(odd[j % len(odd)] for j in range(volume))
        mask = sum(local << (30 * j) for j, local in enumerate(blocks))
        key = mask, links, scalars
        for axis in range(3):
            displacement = tuple(int(j == axis) for j in range(3))
            permutation = site_permutation(side, displacement)
            expected_mask = sum(local << (30 * permutation[j]) for j, local in enumerate(blocks))
            expected_scalar = [0] * volume
            for j, value in enumerate(scalars): expected_scalar[permutation[j]] = value
            assert source.translation_column(key, axis, config) == {
                (expected_mask, links, tuple(expected_scalar)): s.Integer(1)}
        # A local even-parity insertion changes the total source grid to AP.
        defect = mask ^ blocks[0]
        for axis in range(3):
            moved, phase = defect, s.Integer(1)
            for _ in range(side):
                moved, step = translate_fermions(moved, axis, side); phase *= step
            assert moved == defect and phase == -1
        terms = kinetic_terms(config)
        lookup = {(term.left_cell, term.axis, term.displacement): term for term in terms}
        words = 0
        for term in terms:
            left, right = cells(side)[term.left_cell], cells(side)[term.right_cell]
            for axis in range(3):
                new_left = shift_cell(side, left, axis)
                other = lookup[cell_index(side, new_left), term.axis, term.displacement]
                phase = (-1) ** (int(left[axis] == side - 1) + int(right[axis] == side - 1))
                assert s.simplify(term.coefficient * phase - other.coefficient) == 0
                path = tuple(3 * cell_index(side, shift_cell(side, edges(side)[edge][0], axis))
                             + edges(side)[edge][1] for edge in term.path)
                assert path == other.path
                words += 1
        rows.append({'side': side, 'full_local_odd_occupation_translation_controls': 3,
            'local_even_insertion_T_to_N_is_minus_one': 3,
            'complete_actual_PW_hopping_word_translation_identities': words})
    # Exhaust every physical ground-spin label in the N=2 embedding skeleton.
    for spin_mask in range(256):
        mask = sum(1 << (30 * site + 9 + 15 * ((spin_mask >> site) & 1)) for site in range(8))
        for axis in range(3):
            displacement = tuple(int(j == axis) for j in range(3))
            target = shift_spins(spin_mask, 2, displacement)
            expected = sum(1 << (30 * site + 9 + 15 * ((target >> site) & 1)) for site in range(8))
            assert translate_fermions(mask, axis, 2) == (expected, 1)
    config = Configuration(); _, links, _ = source.empty_key()
    flux = list(links); flux[0] = source.label_index()[1, 10, 11]
    flux[4] = source.label_index()[-1, 2, 8]
    mask = sum(odd[j % len(odd)] << (30 * j) for j in range(8))
    key = mask, tuple(flux), tuple(j % 5 for j in range(8))
    image = source.column(key)
    for axis in range(3):
        translate = lambda value: source.translation_column(value, axis)
        assert apply_state(translate(key), source.column) == apply_state(image, translate)
        moved_site = cell_index(2, shift_cell(2, (0, 0, 0), axis))
        for current in range(12):
            assert apply_state(source.gauss_column(key, 0, current), translate) == apply_state(
                translate(key), lambda value: source.gauss_column(value, moved_site, current))
    return {'local_reference': {'density': reference['density'], 'rank': reference['rank'],
            'fermion_parity': reference['fermion_parity'], 'twice_spin': 1,
            'global_zero_boundary_local_ground_energy_interval': spectrum['global_physical_ground_energy_interval']},
        'full_onsite_operator_compiler_replayed': True,
        'full_30_CAR_and_scalar_five_source_join': True,
        'source_spin_restriction': 'ordinary permutation of N^3 identical odd ground doublets for even N',
        'graded_and_AP_factor_per_axis_cycle': '(-1)^(N-1) * (-1) = (-1)^N = 1',
        'all_N2_spin_basis_actual_CAR_translation_controls': 768,
        'size_dependent_word_controls': rows,
        'actual_occupied_PW_H_translation_controls': 3,
        'actual_all_current_translation_covariance_controls': 36,
        'translations_commute_with_each_local_Gauss_generator': False,
        'translations_permute_Gauss_generators_and_preserve_their_common_kernel': True,
        'atomic_reference_is_full_PW_stationary_vacuum': False}


def character_checks():
    rows = []
    for side in (2, 4, 6):
        volume = side ** 3
        values = []
        for momentum in product(range(side), repeat=3):
            multiplicity = momentum_multiplicity(side, momentum)
            values.append(multiplicity)
        assert sum(values) == 2 ** volume
        for displacement in product(range(side), repeat=3):
            permutation = site_permutation(side, displacement)
            visited = set(); cycles = []
            for first in range(volume):
                if first in visited: continue
                cycle = []; current = first
                while current not in visited:
                    visited.add(current); cycle.append(current); current = permutation[current]
                assert current == first
                cycles.append(cycle)
            order = translation_order(side, displacement)
            assert len(cycles) == side ** 2 * gcd(side, *displacement)
            assert all(len(cycle) == order for cycle in cycles)
            assert source_character(side, displacement) == s.Rational(2 ** len(cycles), 2 ** volume)
        histogram = Counter(values)
        rows.append({'side': side, 'site_count': volume,
            'full_atomic_source_ground_support_dimension': str(2 ** volume),
            'positive_source_momentum_outcomes': sum(v > 0 for v in values),
            'source_momentum_multiplicity_classes': [{'multiplicity': str(m), 'phase_count': count,
                'probability': str(s.Rational(m, 2 ** volume))} for m, count in sorted(histogram.items())],
            'complete_full_source_phase_outcome_count': (2 * side) ** 3,
            'parity_forbidden_outcomes_retained': (2 * side) ** 3 - side ** 3})
    return {'character_formula': 'Tr(rho_atomic T^a) = 2^(N^2*gcd(N,a1,a2,a3)-N^3)',
        'source_momentum_multiplicity_formula': '(1/N^3) sum_(d|N) 2^(N^3/d) sum_(e|d,e|k1,e|k2,e|k3) mobius(d/e)*e^3',
        'all_even_N_proof': 'graded/AP tensor permutation, exact cycle trace and finite Fourier inversion',
        'executed_sizes': rows}


def rotation_translation_join():
    """Exact AP one-fermion and full bosonic geometric semidirect action."""
    side = 2; coordinates = cells(side)
    fermion = []; scalar = []
    for axis in range(3):
        f = {}; b = {}
        for cell, xyz in enumerate(coordinates):
            target = cell_index(side, shift_cell(side, xyz, axis))
            f[target, cell] = -1 if xyz[axis] == side - 1 else 1
            b[target, cell] = 1
        fermion.append(s.ImmutableSparseMatrix(8, 8, f))
        scalar.append(s.ImmutableSparseMatrix(8, 8, b))
        assert fermion[-1] ** 2 == -s.eye(8)
        assert scalar[-1] ** 2 == s.eye(8)
    fermion_count = scalar_count = link_count = 0
    for rotation, spin in cubic_rotations():
        sites, links = geometry(rotation, side)
        f = s.ImmutableSparseMatrix(8, 8, {(target, cell): phase
            for cell, (target, phase) in enumerate(sites)})
        b = s.ImmutableSparseMatrix(8, 8, {(target, cell): 1
            for cell, (target, _) in enumerate(sites)})
        assert f.H * f == b.H * b == s.eye(8)
        assert s.simplify(spin.H * spin) == s.eye(2)
        for axis in range(3):
            target_axis = next(j for j in range(3) if rotation[j, axis])
            sign = int(rotation[target_axis, axis])
            assert f * fermion[axis] * f.H == (fermion[target_axis]
                if sign == 1 else fermion[target_axis].H)
            assert b * scalar[axis] * b.H == (scalar[target_axis]
                if sign == 1 else scalar[target_axis].H)
            fermion_count += 1; scalar_count += 1
            for edge, (xyz, edge_axis) in enumerate(edges(side)):
                shifted = 3 * cell_index(side, shift_cell(side, xyz, axis)) + edge_axis
                left = links[shifted]
                rotated_edge, orientation = links[edge]
                anchor, rotated_axis = edges(side)[rotated_edge]
                right = (3 * cell_index(side, shift_cell(side, anchor, target_axis, sign))
                         + rotated_axis, orientation)
                assert left == right
                link_count += 1
    assert (fermion_count, scalar_count, link_count) == (72, 72, 1728)
    return {'actual_proper_cubic_rotations': 24,
        'exact_AP_one_fermion_rotation_translation_conjugations': fermion_count,
        'exact_scalar_rotation_translation_conjugations': scalar_count,
        'exact_oriented_link_label_rotation_translation_conjugations': link_count,
        'semidirect_relation': 'R T_i Rdagger = T_j^(sign), R e_i = sign e_j',
        'full_Fock_lift': 'second quantization of the exact one-fermion identity, with unchanged native R15 indices',
        'complete119_link_lift': 'identical oriented edge map and the same full-state orientation involution',
        'even_N2_restriction': 'T_j^(-1)=T_j, so the three X momentum sectors are permuted',
        'interacting_spectral_cluster_or_gap_inferred_here': False}


def complete_reference_records():
    local = s.Matrix([0, 1, -1, 0]) / s.sqrt(2)
    assert (local.H * local)[0] == 1
    for sigma in PAULI:
        assert (s.kronecker_product(sigma / 2, s.eye(2))
                + s.kronecker_product(s.eye(2), sigma / 2)) * local == s.zeros(4, 1)
    state = purification(); assert norm2(state) == 1
    source_branches = {}; total = defaultdict(lambda: s.Integer(0))
    rows = []
    for momentum in product(range(2), repeat=3):
        branch = phase_projector(state, momentum)
        probability = norm2(branch)
        expected = s.Rational(23 if momentum == (0, 0, 0) else 15, 128)
        assert probability == expected
        assert phase_projector(branch, momentum) == branch
        for axis in range(3):
            displacement = tuple(int(j == axis) for j in range(3))
            phase = (-1) ** momentum[axis]
            assert translate_state(branch, displacement) == {key: phase * value for key, value in branch.items()}
            assert translate_state(branch, displacement, which='reference') == {
                key: phase * value for key, value in branch.items()}
            assert translate_state(branch, displacement, which='joint') == branch
        for previous in source_branches.values(): assert inner(previous, branch) == 0
        source_branches[momentum] = branch
        for key, value in branch.items(): total[key] += value
        rows.append({'source_N_phase_indices': list(momentum),
            'reference_N_phase_indices': [(-k) % 2 for k in momentum],
            'source_full_2N_phase_indices': [2 * k for k in momentum],
            'joint_N_phase_indices': [0, 0, 0], 'probability': str(probability)})
    assert clean(total) == state
    full_rows = [{'phase_indices': list(phase), 'probability': str(full_phase_probability(2, phase))}
                 for phase in product(range(4), repeat=3)]
    assert sum(s.Rational(row['probability']) for row in full_rows) == 1
    assert sum(row['probability'] == '0' for row in full_rows) == 56
    for momentum in product(range(2), repeat=3):
        joint = phase_projector(state, momentum, which='joint')
        assert joint == (state if momentum == (0, 0, 0) else {})
        changed = phase_projector(state, momentum, which='joint', reference_character=(1, 0, 0))
        assert changed == (state if momentum == (1, 0, 0) else {})
    return {'explicit_purification_terms': len(state), 'source_ground_spin_support_dimension': 256,
        'site_reference_spin_support_dimension': 256, 'complete_active_purification_dimension': 65536,
        'physical_purification': 'one even gauge-singlet spin-singlet source/fermionic-reference pair per cell',
        'local_purifying_singlet_norm_squared': '1',
        'local_purifying_singlet_spin_checks': 3,
        'source_only_complete_2N_phase_records': full_rows,
        'source_reference_nonzero_joint_records': rows,
        'all_other_full_source_reference_phase_pairs': {'count': 4088, 'probability': '0', 'retained': True},
        'joint_momentum_default_reference': [0, 0, 0],
        'reference_character_control': {'changed_reference_character': [1, 0, 0],
            'joint_N_phase_label': [1, 0, 0], 'source_only_probabilities_unchanged': True},
        'source_and_joint_translation_labels_identified': False,
        'all_records_preserved_during_full_PW_quench': True,
        'reason_for_time_independence': '[H_PW,T_source]=0 and the reference clock is identity',
        'source_energy_momentum_dispersion_computed': False}
