"""Conserved CAR charge, complete donor readout and a finite source witness."""
from collections import Counter, defaultdict

import sympy as s

from model.centered_source_ground import dimer
from model.defect_field_bridge import full_domain


def creator(modes, mode):
    if type(modes) is not int or not 1 <= modes <= 8 or type(mode) is not int or not 0 <= mode < modes:
        raise ValueError('one to eight modes and a valid mode index required')
    out = s.zeros(1 << modes)
    for mask in range(1 << modes):
        if not (mask >> mode & 1):
            out[mask | (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
    return out


def donor_unitary(modes, mode):
    """Full charge-preserving signed permutation with a filled donor input."""
    addition = creator(modes, mode)
    dim = 1 << modes
    unitary = s.zeros(2*dim)
    used = set()
    for mask in range(dim):
        if mask >> mode & 1:
            output, sign = dim+mask, 1
        else:
            output = mask | (1 << mode)
            sign = addition[output, mask]
        unitary[output, dim+mask] = sign
        used.add(output)
    available = defaultdict(list)
    for address in range(2*dim):
        if address not in used:
            available[address.bit_count()].append(address)
    for mask in range(dim):
        output = available[mask.bit_count()].pop(0)
        unitary[output, mask] = 1
    if any(available.values()):
        raise ArithmeticError('charge-sector completion is not exhaustive')
    return unitary


def certificate():
    data = dimer.data()
    identity = s.eye(16)
    charge = data['N']-2*identity
    if charge*data['K'] != data['K']*charge or charge*data['ground'] != s.zeros(16, 1):
        raise ArithmeticError('centered dimer source or sharp reference charge changed')
    channels = []
    for mode in range(4):
        addition = creator(4, mode)
        removal = addition.T
        occupation = addition*removal
        if addition != dimer.annihilator(mode).T:
            raise ArithmeticError('source CAR convention differs')
        if charge*addition-addition*charge != addition or charge*removal-removal*charge != -removal:
            raise ArithmeticError('exact charge-shifting intertwiner failed')
        if addition.T*addition+occupation.T*occupation != identity:
            raise ArithmeticError('failure branch was lost from the donor instrument')
        unitary = donor_unitary(4, mode)
        total = s.diag(*(address.bit_count()-2 for address in range(32)))
        if unitary.T*unitary != s.eye(32) or total*unitary != unitary*total:
            raise ArithmeticError('full donor completion loses unitarity or total charge')
        if unitary[:16, 16:] != addition or unitary[16:, 16:] != occupation:
            raise ArithmeticError('donor pointer branches do not equal the declared instrument')
        g, norm = data['ground'], data['ground_norm_squared']
        weights = [s.simplify(((op*g).T*(op*g))[0]/norm) for op in (addition, removal, occupation)]
        if weights != [s.Rational(1, 2)]*3:
            raise ArithmeticError('actual centered dimer response weight changed')
        channels.append({'mode': mode, 'addition_weight': '1/2', 'removal_weight': '1/2',
                         'creation_instrument_failure_weight': '1/2'})
    branches = []
    for sites in (2, 3):
        for law in ('C', 'E'):
            rows = full_domain.finite_branches(sites=sites, flavors=2, law=law, seam=(0, sites-1))
            if set(rows) != {(coin, raw) for coin in (0, 1) for raw in range(full_domain.RAW_SLOTS)}:
                raise ArithmeticError('complete raw record slots were removed')
            columns = Counter()
            for values in rows.values():
                if len({out for _, out, _ in values}) != len(values):
                    raise ArithmeticError('a raw branch is not a partial isometry')
                for initial, final, sign in values:
                    if initial.bit_count() != final.bit_count() or sign not in (-1, 1):
                        raise ArithmeticError('source repair changed charge in a raw branch')
                    columns[initial] += 1
            if dict(columns) != {mask: 2 for mask in range(1 << (2*sites))}:
                raise ArithmeticError('complete fair-coin effects do not sum to identity')
            branches.append({'sites': sites, 'law': law, 'raw_slots': len(rows),
                             'all_input_columns': 1 << (2*sites), 'charge_preserved': True,
                             'complete_identity_effect_verified': True})
    spin_raise = creator(4, 0)*creator(4, 1).T
    if charge*spin_raise != spin_raise*charge:
        raise ArithmeticError('spin excitation changed occupation charge')
    rho = s.simplify(data['ground']*data['ground'].T/data['ground_norm_squared'])
    repaired = s.zeros(16)
    for values in full_domain.finite_branches(sites=2, flavors=2, law='C', seam=(0, 1)).values():
        branch = s.zeros(16)
        for initial, final, sign in values:
            branch[final, initial] = sign
        repaired += branch*rho*branch.T/2
    ground_return = s.simplify(s.trace(rho*repaired))
    if s.simplify(1-ground_return).is_positive is not True:
        raise ArithmeticError('complete repair ground-disturbance countercontrol lost')
    return {
        'source': 'centered Weyl/contact dimer and specified two-flavour full-domain C/E extension',
        'invariant': 'Q=N-sites; conserved occupation excess, not a topological invariant',
        'source_clock_commutes_with_Q': True,
        'actual_finite_reference_has_sharp_Q0': True,
        'sharp_charge_inferred_from_one_point_mean': False,
        'creator_and_remover_charge_shifts': [1, -1],
        'full_donor_dimension': 32,
        'full_donor_unitary_preserves_total_charge': True,
        'all_creation_success_and_failure_outcomes_retained': True,
        'actual_centered_reference_weights': channels,
        'complete_raw_branch_checks': branches,
        'protected_spin_excitations_have_occupation_charge_zero': True,
        'centered_ground_return_after_unread_complete_C_repair': str(ground_return),
        'charge_conservation_implies_reference_stationarity': False,
        'native_donor_pulse_word_executed': False,
        'real_space_topological_defect_identified': False,
        'momentum_Berry_charge_used_as_real_space_defect': False,
        'relativistic_particle_or_mass_identified': False,
    }
