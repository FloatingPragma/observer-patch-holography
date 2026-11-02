"""Source checks and exhaustive cut certificates for twelve native ports."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb

import sympy as sp

from model.exchange_sector import many_body as native
from . import bounds


@lru_cache(maxsize=1)
def source_columns():
    """Check every native address, raw record and output, for both signs/laws."""
    checked = 0
    for law in ('C', 'E'):
        for a, b in native.SEAMS:
            for coin in (0, 1):
                for mask in range(native.DIMENSION):
                    left, right = (mask >> a) & 1, (mask >> b) & 1
                    moving = left != right and ((left == 1 and coin == 0) or (right == 1 and coin == 1))
                    expected = mask ^ ((1 << a) | (1 << b)) if moving else mask
                    raw = (1 if left and moving else 2 if moving else 0)
                    # Native raw1 is the a->b branch and raw2 is b->a.
                    if law == 'E':
                        raw += 4*(left != right)
                    for character in (-1, 1):
                        actual_raw, output, phase, _ = native.column(mask, law, (a, b), coin, character)
                        if actual_raw != raw or output != expected or phase not in (-1, 1):
                            raise ArithmeticError('native source differs from the complete endpoint table')
                        if character == 1 and phase != 1:
                            raise ArithmeticError('positive lift acquired an unaccounted phase')
                        checked += 1
    return {'complete_native_columns_checked': checked, 'all_addresses': native.DIMENSION,
            'laws': ['C', 'E'], 'characters': [-1, 1], 'reserved_abort_slots_have_no_support': True}


def local_overlap_table():
    a, b = native.SEAMS[0]
    masks = (0, 1 << a, 1 << b, (1 << a) | (1 << b))
    result = {}
    for law in ('C', 'E'):
        table = []
        for first in masks:
            row = []
            for second in masks:
                actual = F(sum(native.column(first, law, (a, b), c, 1)[0]
                               == native.column(second, law, (a, b), c, 1)[0] for c in (0, 1)), 2)
                da = ((first >> a) & 1)-((second >> a) & 1)
                db = ((first >> b) & 1)-((second >> b) & 1)
                expected = 1-(F(abs(da-db), 2) if law == 'C' else int(da != db))
                if actual != expected:
                    raise ArithmeticError('complete-record coefficient formula failed')
                row.append(str(actual))
            table.append(row)
        result[law] = table
    return result


def native_cuts():
    cuts = bounds.cut_table(native.MODES, native.SEAMS)
    floor = min(cuts[1:-1])
    if floor != F(1, 6):
        raise ArithmeticError('native weighted min-cut changed')
    # Every disjoint +/- pair represents the complete ternary difference field.
    minimum, fixed_minimum, count = F(2), F(2), 0
    full = native.FULL
    for positive in range(1 << native.MODES):
        negative = full ^ positive
        while True:
            count += 1
            if (positive, negative) not in ((0, 0), (full, 0), (0, full)):
                loss = (cuts[positive]+cuts[negative])/2
                minimum = min(minimum, loss)
                if positive.bit_count() == negative.bit_count():
                    fixed_minimum = min(fixed_minimum, loss)
            if negative == 0:
                break
            negative = (negative-1) & (full ^ positive)
    if count != 3**native.MODES or minimum != floor/2 or fixed_minimum != floor:
        raise ArithmeticError('exhaustive coarea certificate failed')
    return {'vertices': native.MODES, 'edges': len(native.SEAMS),
            'nontrivial_cuts_checked': len(cuts)-2, 'ternary_differences_checked': count,
            'minimum_cut': str(floor), 'fixed_number_loss_floor': str(fixed_minimum),
            'nonextremal_full_carrier_loss_floor': str(minimum),
            'fixed_number_per_repair_factor': str(1-floor),
            'full_carrier_nonextremal_per_repair_factor': str(1-floor/2)}


def diagonal_gap():
    A = sp.zeros(native.MODES)
    for a, b in native.SEAMS:
        A[a, b] = A[b, a] = 1
    x = sp.Symbol('x')
    expected = (x-5)*(x+1)**5*(x*x-5)**3
    actual = A.charpoly(x).as_expr()
    if sp.expand(actual-expected) != 0 or any(sum(A[i, j] for j in range(native.MODES)) != 5 for i in range(native.MODES)):
        raise ArithmeticError('exact native one-particle spectrum changed')
    # sqrt(5)<9/4, hence (5-sqrt(5))/60 > 11/240.
    if F(9, 4)**2 <= 5:
        raise ArithmeticError('rational spectral lower bound failed')
    return {'adjacency_characteristic_polynomial': str(sp.factor(actual)),
            'all_nontrivial_number_sectors': list(range(1, native.MODES)),
            'fixed_selector_diagonal_gap': '(5-sqrt(5))/60',
            'strict_rational_gap_floor': '11/240',
            'imported_theorem': 'Caputo-Liggett-Richthammer (2010), Theorem1.1 and Eq4.1',
            'primary_source': 'https://arxiv.org/abs/0906.1238',
            'imported_theorem_reproved_here': False,
            'generator_edge_rates': '1/60', 'clock': 'one uniformly selected edge per repair tick',
            'faithful_stationary_sector_state': 'I_N/binomial(12,N)',
            'ground_state_or_physical_temperature_selected': False}


def report():
    cuts = bounds.cut_table(native.MODES, native.SEAMS)
    rows = []
    for number in (1, 6, 11):
        ticks = 0
        while bounds.mixing_upper(12, number, ticks, F(11, 240), F(1, 6)) > F(1, 100):
            ticks += 1
        rows.append({'number': number, 'dimension': comb(12, number), 'ticks_for_trace_norm_below_0_01': ticks,
                     'exact_trace_norm_upper': str(bounds.mixing_upper(12, number, ticks, F(11, 240), F(1, 6)))})
    return {'source': source_columns(), 'record_overlap_tables': local_overlap_table(),
            'cuts': native_cuts(), 'diagonal': diagonal_gap(),
            'two_configuration_superposition_26_tick_offdiagonal_norm_upper': str(bounds.coherence_upper(1, cuts, 26)),
            'fixed_selector_full_sector_mixing': rows,
            'stationary_mixing_discards_classical_tape': True,
            'adaptive_classical_history_coherence_bound': True,
            'adaptive_controller_stationary_state_inferred': False,
            'coherent_active_recovery_excluded_by_theorem': False,
            'mass_or_Hamiltonian_gap_identified': False,
            'physical_particle_species_identified': False,
            'unrestricted_OPH_refinement_excluded': False}
