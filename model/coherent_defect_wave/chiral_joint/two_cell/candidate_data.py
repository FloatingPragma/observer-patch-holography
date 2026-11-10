"""Auditable candidate labels for the finite source, without species guesses."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ..model import fermion
from ...chiral_source.model import native_basis
from ...chiral_higgs.model import annihilator, higgs_memory, quantum_lepton_block
from .basis import current_data
from .hamiltonian import PARAMETERS, conserved_labels
from .instrument import SOURCE_MODES, branch_metadata
from .symmetries import (conjugating_mask, spin_reversal_mask, PERIODIC_DIMENSION)


@lru_cache(maxsize=None)
def _field(mask, mode, create):
    result = fermion(mask, mode, create)
    return None if result is None else (result[0], int(result[1]))


def _pair(mask, left, right):
    first = _field(mask, *right)
    if first is None:
        return None
    second = _field(first[0], *left)
    return None if second is None else (second[0], first[1] * second[1])


def _sum_terms(*terms):
    total = defaultdict(int)
    for term in terms:
        if term is not None:
            total[term[0]] += term[1]
    return {key: value for key, value in total.items() if value}


@lru_cache(maxsize=1)
def charged_operator_checks():
    matter, _, _ = current_data()
    ec = s.Matrix([1, 0, 0])
    weights = []
    for t in matter:
        weight = t[0, 0]
        assert t * ec == weight * ec
        weights.append(str(weight))
    assert weights == ['0'] * 11 + ['1']
    assert sum((matter[j] ** 2 for j in range(8)), s.zeros(3)) * ec == s.zeros(3, 1)
    assert sum((matter[j] ** 2 for j in (8, 9, 10)), s.zeros(3)) * ec == s.zeros(3, 1)
    car = conjugate = reversal = 0
    for mask in range(4096):
        for a in SOURCE_MODES:
            for b in SOURCE_MODES:
                result = _sum_terms(_pair(mask, (a, False), (b, True)),
                                    _pair(mask, (b, True), (a, False)))
                assert result == ({mask: 1} if a == b else {})
                assert not _sum_terms(_pair(mask, (a, False), (b, False)),
                                      _pair(mask, (b, False), (a, False)))
                car += 2
            # A c^dagger_(site,spin) A^-1=-i c_(other_site,spin).
            image, phase = conjugating_mask(mask)
            first = _field(mask, a, True)
            after = _field(image, (a + 6) % 12, False)
            assert (first is None) == (after is None)
            if first is not None:
                target, last = conjugating_mask(first[0])
                assert target == after[0]
                assert s.expand(first[1] * last + s.I * phase * after[1]) == 0
            conjugate += 1
        target, phase = spin_reversal_mask(mask)
        back, second = spin_reversal_mask(target)
        ell, spin = conserved_labels((mask, 0, 0, 0))
        assert back == mask and phase * second == (-1) ** mask.bit_count()
        assert (-1) ** mask.bit_count() == (-1) ** (ell % 2)
        assert conserved_labels((target, 0, 0, 0)) == (ell, -spin)
        reversal += 1
    return {'native_R15_index': 9, 'native_subset': list(native_basis()[9]),
            'native_cloud_degree': len(native_basis()[9]),
            'exact_gauge_generator_weights': weights,
            'SU3_representation': 'singlet', 'SU2_representation': 'singlet', 'hypercharge': '1',
            'algebraic_T3_plus_Y': '1', 'physical_EM_current_identified': False,
            'CAR_anticommutator_column_checks': car,
            'source_occupation_domain': 4096,
            'reflected_conjugate_operator_columns': conjugate,
            'spin_reversal_and_parity_columns': reversal,
            'CAR_grading_is_declared_source_architecture': True,
            'one_internal_ec_channel_repeated_over_two_sites_and_two_spin_components': True,
            'generation_count_or_four_distinct_species_inferred': False}


@lru_cache(maxsize=1)
def neutral_operator_checks():
    h, gauss, _, _, _ = quantum_lepton_block()
    _, _, scalar_number, _ = higgs_memory()
    number = s.ImmutableSparseMatrix(s.kronecker_product(s.eye(64), scalar_number))
    yukawa = s.ImmutableSparseMatrix(h - number)
    cs = [annihilator(j) for j in range(6)]
    ell = s.diag(*[sum((-1, 1, 1, -1, 1, 1)[j] for j in range(6) if mask >> j & 1)
                   for mask in range(64)])
    sz = s.diag(*[sum((1, 1, 1, -1, -1, -1)[j] for j in range(6) if mask >> j & 1)
                  for mask in range(64)])
    raising = sum((cs[j].H * cs[j + 3] for j in range(3)), s.zeros(64))
    parity = s.diag(*[(-1) ** mask.bit_count() for mask in range(64)])
    extra = [s.ImmutableSparseMatrix(s.kronecker_product(matrix, s.eye(5)))
             for matrix in (ell, sz, raising, parity)]
    checks = 0
    for operator in (number, yukawa):
        assert operator == operator.H
        for generator in [*gauss, *extra]:
            assert generator * operator == operator * generator
            checks += 1
    assert number * number == number
    assert number * h != h * number
    assert yukawa * h != h * yukawa
    assert yukawa[:, 0] != s.zeros(320, 1)
    return {'complete_local_domain': 320, 'actual_operator_current_commutators': checks,
            'scalar_occupation': {'operator': 'N_H at the specified cell',
                'complete_effects': ['N_H', 'I-N_H'], 'bounded_projector': True,
                'delta_L': 0, 'delta_twoSz': 0, 'gauge_charge': 0, 'fermion_parity': 'even'},
            'Yukawa_composite': {'operator': 'local H_Y on the same source',
                'complete_effects': ['(I+H_Y/8)/2', '(I-H_Y/8)/2'],
                'safe_operator_norm_upper_bound': '8',
                'bound_reason': 'four unit-norm pair/scalar terms and their four adjoints',
                'delta_L': 0, 'delta_twoSz': 0, 'gauge_charge': 0, 'fermion_parity': 'even'},
            'both_operators_have_nonzero_same_H_commutators': True,
            'local_source_spin_scalar_commutators_checked': True,
            'source_spin_scalar_implies_relativistic_scalar_particle': False,
            'neutral_spectral_weights_computed_by_this_module': False}


def native_family_matrix():
    matter, _, _ = current_data()
    weak_casimir = sum((matter[j] ** 2 for j in (8, 9, 10)), s.zeros(3))
    color_casimir = sum((matter[j] ** 2 for j in range(8)), s.zeros(3))
    assert weak_casimir == s.diag(0, s.Rational(3, 4), s.Rational(3, 4))
    assert color_casimir == s.zeros(3)
    rows = []
    for i, native in enumerate((9, 10, 11)):
        hypercharge = matter[11][i, i]
        weak_weight = matter[10][i, i]
        rows.append({'native_R15_index': native, 'native_subset': list(native_basis()[native]),
                     'channel': 'e^c' if i == 0 else 'L weak-doublet component',
                     'hypercharge': str(hypercharge), 'weak_T3': str(weak_weight),
                     'algebraic_T3_plus_Y': str(hypercharge + weak_weight),
                     'weak_Casimir': str(weak_casimir[i, i]), 'color_Casimir': '0',
                     'source_L': -1 if i == 0 else 1})
    return {'rows': rows,
            'actual_weak_doublet_generators': [[[str(v) for v in row] for row in t[1:3, 1:3].tolist()]
                                              for t in matter[8:11]],
            'quark_sector': 'retained empty reducing spectators with both quark Yukawa couplings zero',
            'quark_particle_spectra_calculated_here': False,
            'physical_EM_field_or_lepton_identities_inferred': False}


def candidate_data():
    charged = charged_operator_checks()
    neutral = neutral_operator_checks()
    probes = []
    for mode in SOURCE_MODES:
        probes.append({'id': f'ec-cell{mode // 6}-spin{mode % 6 // 3}',
                       'native_R15_index': charged['native_R15_index'],
                       'source_mode': mode,
                       'charged_add_branch': branch_metadata(mode, 1, 0),
                       'charged_remove_branch': branch_metadata(mode, 0, 1),
                       'status': 'exact source operator channel and complete compensated energy instrument',
                       'physical_species_assignment': None})
    return {'source': {'carrier': 'uniform open two-cell lepton/quantum-link/quantum-Higgs control',
                       'dimension': 614400, 'parameters': PARAMETERS,
                       'reference': 'unique physical neutral ground of this same H, L=2 and 2Sz=0',
                       'original_ground_preparation_probability': '1/1053',
                       'clock': 'exp(-it H) for the same declared dimensionless Hamiltonian'},
            'charged_channel_checks': charged, 'charged_probes': probes,
            'native_active_family': native_family_matrix(),
            'neutral_channels': neutral,
            'exact_symmetry_scope': {
                'odd_L_Kramers_pairs': 'same-boundary Theta symmetry; Theta^2=(-1)^L',
                'reflected_conjugate': 'right addition and left removal with reflected opposite boundary; not same-side equality',
                'physical_Lorentz_spin_or_CPT_proved': False},
            'separate_periodic_control': {'dimension': PERIODIC_DIMENSION,
                'two_independent_quantum_links': True,
                'finite_translation_phase_PVM': True,
                'open_ground_or_energy_spectrum_transferred': False,
                'larger_even_ring_menu_requirements': {'two_distinct_momentum_magnitudes': 'at least 4 cells',
                    'three_distinct_momentum_magnitudes': 'at least 6 cells',
                    'menu_counts_suffice_for_a_particle_dispersion': False},
                'periodic_translation_phase_is_3D_continuum_momentum': False},
            'unresolved_classification': ['elementary versus composite or boundary-bound excitation',
                'physical electron, muon, tau or their antiparticle identification',
                'three-dimensional momentum and relativistic little-group spin',
                'volume/cutoff-stable isolated particle pole and residue',
                'relativistic rest mass and physical energy calibration'],
            'empirical_fit_or_confidence_assigned': False}


def certificate_data():
    """Stable JSON-ready classification packet for the owning receipt."""
    return candidate_data()
