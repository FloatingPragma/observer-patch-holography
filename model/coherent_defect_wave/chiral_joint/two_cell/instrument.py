"""Complete compensated fermionic reservoir instrument on the same source.

All 16384 fermion columns are retained. Link/scalar factors are spectators
during the driven swap, and the retained boundary charges select the correct
source Gauss sector for each outcome. No apparatus work value is inferred.
"""
from functools import lru_cache

import sympy as s

from .basis import current_data, complete_metadata
from .spectrum import sector_certificate


SOURCE_MODES = (0, 3, 6, 9)
FERMION_DIMENSION = 1 << 14


def _probe(mode):
    if type(mode) is not int or mode not in SOURCE_MODES:
        raise ValueError('source e^c mode must be 0, 3, 6 or 9')
    return mode // 6, 1 if mode % 6 == 0 else -1


def _annihilate(mask, mode):
    if not mask >> mode & 1:
        return None
    return mask ^ (1 << mode), (-1) ** ((mask & ((1 << mode) - 1)).bit_count())


def _create(mask, mode):
    if mask >> mode & 1:
        return None
    return mask | (1 << mode), (-1) ** ((mask & ((1 << mode) - 1)).bit_count())


@lru_cache(maxsize=4)
def swap_table(source_mode):
    """Actual exterior lift of the source/reservoir transposition."""
    _probe(source_mode)
    permutation = list(range(14))
    permutation[source_mode], permutation[12] = 12, source_mode
    columns = []
    for mask in range(FERMION_DIMENSION):
        occupied = [permutation[i] for i in range(14) if mask >> i & 1]
        sign = (-1) ** sum(a > b for k, a in enumerate(occupied) for b in occupied[k + 1:])
        columns.append((sum(1 << i for i in occupied), sign))
    return tuple(columns)


def reduced_kraus_column(source_mask, source_mode, setting, reservoir_out):
    """Unnormalized branch column; setting probability 1/2 is separate."""
    if type(source_mask) is not int or not 0 <= source_mask < 4096:
        raise ValueError('source mask is an exact twelve-mode occupation')
    if setting not in (0, 1) or reservoir_out not in (0, 1):
        raise ValueError('setting and reservoir_out are binary')
    full = source_mask | (setting << 12) | (setting << 13)
    target, phase = swap_table(source_mode)[full]
    assert target >> 13 & 1 == setting
    if target >> 12 & 1 != reservoir_out:
        return None
    return target & 4095, phase


def branch_metadata(source_mode, setting, reservoir_out):
    site, spin = _probe(source_mode)
    charge = reservoir_out - setting
    boundary = [0, 0]
    boundary[site] = charge
    return {'source_mode': source_mode, 'source_site': site,
            'setting': 'add' if setting else 'remove', 'initial_r_b': [setting, setting],
            'final_r_b': [reservoir_out, setting], 'initial_setting_probability': '1/2',
            'boundary_hypercharges': boundary, 'source_delta_L': charge,
            'source_delta_twoSz': -spin * charge,
            'charged_transition': charge != 0,
            'ground_source_L_after': 2 + charge, 'ground_source_twoSz_after': -spin * charge}


@lru_cache(maxsize=4)
def swap_checks(source_mode):
    site, spin = _probe(source_mode)
    table = swap_table(source_mode)
    permutation = list(range(14))
    permutation[source_mode], permutation[12] = 12, source_mode
    intertwiners = 0
    for mask, (target, phase) in enumerate(table):
        reverse, reverse_phase = table[target]
        assert reverse == mask and phase * reverse_phase == 1
        assert target.bit_count() == mask.bit_count()
        assert target >> 13 & 1 == mask >> 13 & 1
        for mode in range(14):
            step = _annihilate(target, mode)
            if step is None:
                left = None
            else:
                output, last_phase = table[step[0]]
                left = output, phase * step[1] * last_phase
            assert left == _annihilate(mask, permutation[mode])
            intertwiners += 1
    matter, _, _ = current_data()
    gauss_checks = 0
    # Equality of actual one-body generators plus all annihilator
    # intertwiners proves the full exterior-current commutator exactly.
    for endpoint in (0, 1):
        for t in matter:
            generator = s.zeros(14)
            for source_spin in range(2):
                offset = 6 * endpoint + 3 * source_spin
                generator[offset:offset + 3, offset:offset + 3] = t
            if endpoint == site:
                generator[12, 12] = t[0, 0]
                generator[13, 13] = -t[0, 0]
            assert generator.extract(permutation, permutation) == generator
            gauss_checks += 1
    lepton = (-1, 1, 1) * 4 + (-1, 1)
    twice_spin = (1, 1, 1, -1, -1, -1) * 2 + (spin, -spin)
    for weights in (lepton, twice_spin):
        assert tuple(weights[j] for j in permutation) == weights
    return {'source_mode': source_mode, 'full_fermion_columns': FERMION_DIMENSION,
            'actual_annihilator_column_intertwiners': intertwiners,
            'actual_local_Gauss_current_checks': gauss_checks,
            'exact_unitary_and_compensator_invariance': True,
            'total_L_twoSz_and_parity_conserved': True,
            'source_plus_apparatus_dimension': FERMION_DIMENSION * 6 * 5 ** 2}


@lru_cache(maxsize=4)
def reduced_checks(source_mode):
    _probe(source_mode)
    even_columns = 0
    completeness = [0] * 4096
    for setting in (0, 1):
        for outcome in (0, 1):
            targets = set()
            for mask in range(4096):
                column = reduced_kraus_column(mask, source_mode, setting, outcome)
                if column is not None:
                    assert column[0] not in targets
                    targets.add(column[0])
                    assert column[1] in (-1, 1)
                    completeness[mask] += 1
                if mask.bit_count() % 2:
                    continue
                occupied = bool(mask >> source_mode & 1)
                if setting == 1 and outcome == 0:
                    expected = _create(mask, source_mode)
                elif setting == 1 and outcome == 1:
                    expected = (mask, -1) if occupied else None
                elif setting == 0 and outcome == 1:
                    removed = _annihilate(mask, source_mode)
                    expected = (removed[0], -removed[1]) if removed else None
                else:
                    expected = (mask, 1) if not occupied else None
                assert column == expected
                even_columns += 1
    assert all(value == 2 for value in completeness)
    return {'source_mode': source_mode, 'all_source_columns_in_completeness': 4096,
            'even_reference_branch_columns_checked': even_columns,
            'original_setting_probability': '1/2',
            'even_source_Kraus_in_add_success_add_failure_remove_success_remove_failure_order':
                ['a^dagger/sqrt(2)', '-n/sqrt(2)', '-a/sqrt(2)', '(I-n)/sqrt(2)'],
            'all_source_parity_instrument_complete': True,
            'unqualified_cross_parity_Kraus_formula_asserted': False,
            'branches': [branch_metadata(source_mode, setting, outcome)
                         for setting, outcome in ((1, 0), (1, 1), (0, 1), (0, 0))]}


@lru_cache(maxsize=1)
def ground_preparation():
    dimension = complete_metadata(0)['physical_dimension']
    ground = sector_certificate(0)
    assert dimension == 1053
    assert ground['ground_is_certifiably_simple']
    assert ground['ground_window_candidate_labels'] == [{'L': 2, 'twoSz': 0, 'eigen_index': 0}]
    return {'physical_dimension': dimension, 'source_input': 'P/1053',
            'preparation_ancilla': 'retained neutral graded register, L_anc=-L_j, twoSz_anc=-twoSz_j, parity_anc=parity_j',
            'total_preparation_L_twoSz_parity': [0, 0, 'even'],
            'ground_projector_rank': 1, 'ground_L_twoSz_parity': [2, 0, 'even'],
            'complete_preparation_effects': ['P_ground', 'I-P_ground'],
            'original_ground_success_probability': '1/1053',
            'original_preparation_failure_probability': '1052/1053',
            'conditional_success_state': 'the unique physical ground of the same two-cell Hamiltonian',
            'ground_interval': ground['global_ground_interval'],
            'post_ground_ancilla_L_twoSz_parity': [-2, 0, 'even'],
            'efficiency_or_native_idle_preparation_claim': False,
            'complete_energy_readout_after_swap': True,
            'apparatus_work_calculated': False}


def certificate_data():
    return {'physical_ground_preparation': ground_preparation(),
            'source_modes': list(SOURCE_MODES),
            'actual_swaps': [swap_checks(mode) for mode in SOURCE_MODES],
            'complete_reduced_instruments': [reduced_checks(mode) for mode in SOURCE_MODES],
            'all_failure_outcomes_retained': True,
            'source_energy_measurements_use_the_same_H': True,
            'phased_physical_embedding_required': 'spectrum.physical_block returns D B',
            'finite_charged_energy_is_relativistic_rest_mass': False}
