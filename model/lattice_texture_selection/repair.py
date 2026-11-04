"""Complete native repair records and the separate hopping countercontrol."""
from model.defect_field_bridge import full_domain
from model.native_repair_spin_limit.algebra import hopping_twice, protected

import numpy as np


def audit(sites, law):
    if type(sites) is not int or not 2 <= sites <= 3 or law not in ('C', 'E'):
        raise ValueError('the exact oracle covers two or three sites and C/E repair')
    seam = (0, sites-1)
    dimension = 1 << (2*sites)
    rows = full_domain.finite_branches(sites=sites, flavors=2, law=law, seam=seam)
    input_counts = [0]*dimension
    output_counts = [0]*dimension
    protected_count = held_collision = 0
    for coin in (0, 1):
        for mask in range(dimension):
            raw, output, sign, current = full_domain.branch_action(
                mask, law, coin, sites=sites, flavors=2, seam=seam)
            input_counts[mask] += 1
            output_counts[output] += 1
            is_protected = all(((mask >> (2*x)) & 3).bit_count() == 1 for x in range(sites))
            if is_protected:
                protected_count += 1
                if (raw, output, sign, current) != (0, mask, 1, 0):
                    raise ArithmeticError('the native protected branch is not identity with raw zero')
            held_collision += raw == full_domain.COLLISION_RAW
    if set(input_counts) != {2} or set(output_counts) != {2}:
        raise ArithmeticError('complete native branch completeness changed')
    if protected_count != 2*(1 << sites) or held_collision == 0 or len(rows) != 32:
        raise ArithmeticError('protected or complete nonbinary input domain changed')
    T2 = hopping_twice(sites, tuple((x, x+1, 1) for x in range(sites-1)))
    P = protected(sites)
    leakage = T2[np.ix_(~P, P)]
    if not np.any(leakage) or np.any(T2[np.ix_(P, P)]):
        raise ArithmeticError('source hopping must leave some bare protected spin states')
    return {'sites': sites, 'law': law, 'seam': list(seam),
            'full_Fock_input_dimension': dimension, 'retained_coin_raw_slots': len(rows),
            'complete_input_and_output_multiplicity': 2,
            'protected_coin_columns': protected_count, 'held_collision_columns': held_collision,
            'protected_complete_isometry': 'identity_source tensor |+>_coin tensor |0>_raw',
            'arbitrary_untouched_reference_retained': True,
            'hopping_outside_protected_Frobenius_squared_for_T2': int(np.sum(leakage*leakage)),
            'repair_only_records_select_spin_direction': False,
            'interleaved_Hamiltonian_evolution_is_identity': False}


def certificate():
    return {'full_native_audits': [audit(n, law) for n in (2, 3) for law in ('C', 'E')],
            'arbitrary_finite_repair_only_history_identity': True,
            'fresh_coin_raw_memory_is_required': True,
            'controlled_rotations_preserve_bare_protected_carrier': True,
            'controlled_rotations_are_selected_by_repair': False,
            'dressed_finite_c_ground_or_full_process_stationarity_claimed': False}
