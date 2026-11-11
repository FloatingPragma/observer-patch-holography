"""Complete finite readouts of the calibrated quadratic action."""
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement

import numpy as np
import sympy as s
from scipy import sparse

from .source import fermions, bosons, scalar, MOMENTA
from model.common_source import calibrated_matter as matter


def car_action(mask, mode, create, number_modes=180):
    """Literal CAR on the entire finite Fock basis, without a particle cap."""
    if not 0 <= mask < (1 << number_modes) or not 0 <= mode < number_modes:
        raise ValueError('occupation or mode outside full CAR source')
    occupied = (mask >> mode) & 1
    if occupied == create:
        return None
    sign = -1 if (mask & ((1 << mode)-1)).bit_count() % 2 else 1
    return mask ^ (1 << mode), sign


def fermion_energy(mask, mode_energies):
    if not 0 <= mask < 1 << len(mode_energies):
        raise ValueError('occupation outside source')
    return sum((mode_energies[j] for j in range(len(mode_energies)) if mask >> j & 1), s.S.Zero)


@lru_cache(maxsize=1)
def fermion_records():
    f = fermions()
    records = []
    for group in f['groups']:
        for create in (True, False):
            probability = s.Rational(len(group['modes']), 180) if create else s.S.Zero
            records.append({'species': group['species'], 'orientation': group['orientation'],
                'momentum_GeV': group['momentum_GeV'],
                'pointer': 'create' if create else 'remove',
                'source_energy_GeV': group['energy_GeV'] if create else scalar(0),
                'reference_energy_GeV': 1,
                'probability': scalar(probability),
                'coherent_reference_components': len(group['modes']),
                'source_electric_charge': group['electric_charge'] if create else '0',
                'reference_electric_charge_on_success': str(-s.sympify(group['electric_charge']) if create else s.sympify(group['electric_charge'])),
                'other_source_energy_projectors_probability': '0'})
    probabilities = [s.sympify(row['probability']['exact']) for row in records]
    assert sum(probabilities) == 1
    return {'source_CAR_modes': 180, 'source_Hilbert_dimension_exact': str(1 << 180),
        'Hamiltonian': 'H_F=sum_j omega_j a_j^dagger a_j, vacuum energy subtracted',
        'source_preparation': 'the Bogoliubov vacuum of the actual native momentum-pair quadratic source',
        'full_energy_projectors': 'P_E=sum_{occupation masks n: sum_j omega_j n_j=E}|n><n|; degeneracies are combined, no many-body state is discarded',
        'complete_instrument': 'K_{g,+}=sum_{j in g} a_j^dagger tensor |dual j>/sqrt(180), K_{g,-}=sum_{j in g} a_j tensor |j>/sqrt(180)',
        'all_input_identity': 'sum_g(K_g+^dagger K_g+ + K_g-^dagger K_g-)=sum_j{a_j,a_j^dagger}/180=I on all 2^180 source states',
        'reference_dimension': 361, 'reference_blank_energy_GeV': 0,
        'reference_nonblank_energy_GeV': 1,
        'reference_charge_rule': 'dual representation for creation and original representation for annihilation',
        'native_field_implementation': 'canonical transformation of the actual 360-dimensional Nambu symbol; the instrument is a declared Stinespring probe of its normalized quasiparticle fields',
        'pointer_record_count': len(records), 'positive_record_count': sum(bool(x > 0) for x in probabilities),
        'zero_record_count': sum(bool(x == 0) for x in probabilities),
        'records': records, 'probability_sum': '1',
        'zero_outside_projector_retained_per_pointer': True,
        'fermion_occupation_truncation': False}


@lru_cache(maxsize=1)
def boson_operators():
    """Project the oscillator creators to the total-number <=2 space."""
    n = len(bosons()['modes'])
    states = [()] + [(j,) for j in range(n)] + list(combinations_with_replacement(range(n), 2))
    index = {state: i for i, state in enumerate(states)}
    creators = []
    for mode in range(n):
        rows, columns, values = [], [], []
        for column, state in enumerate(states):
            if len(state) < 2:
                rows.append(index[tuple(sorted(state+(mode,)))])
                columns.append(column)
                values.append(np.sqrt(1+state.count(mode)))
        creators.append(sparse.csr_matrix((values, (rows, columns)), shape=(len(states), len(states))))
    summed = sum((a.T@a+a@a.T for a in creators), sparse.csr_matrix((len(states), len(states))))
    expected = np.array([56 if len(x) == 0 else 58 if len(x) == 1 else 2 for x in states])
    error = sparse.linalg.norm(summed-sparse.diags(expected))
    if error > 2e-12:
        raise ArithmeticError('global-cap oscillator identity differs')
    no_operation_squared = 1-expected/58
    completion = summed/58+sparse.diags(no_operation_squared)
    completion_error = sparse.linalg.norm(completion-sparse.eye(len(states)))
    if completion_error > 2e-12:
        raise ArithmeticError('complete boson instrument failed on an input column')
    return states, creators, no_operation_squared, float(error), float(completion_error)


@lru_cache(maxsize=1)
def boson_records():
    b = bosons()
    states, creators, no_operation_squared, cap_error, complete_error = boson_operators()
    labels = ('W', 'Z', 'Higgs', 'massless')
    category = {mode['index']: (mode['species'] if mode['species'] in labels[:3] else 'massless') for mode in b['modes']}
    energies = {name: s.sqrt(1+b['masses'][name]**2) for name in labels[:3]}
    energies['massless'] = s.S.One
    energy_counts = Counter(tuple(sum(category[j] == label for j in state) for label in labels) for state in states)
    sectors = []
    for index, (counts, rank) in enumerate(sorted(energy_counts.items(), key=lambda item: (sum(item[0]), item[0]))):
        energy = sum((n*energies[name] for n, name in zip(counts, labels)), s.S.Zero)
        sectors.append({'index': index, 'occupation_by_energy_class': dict(zip(labels, counts)),
                        'energy_GeV': scalar(energy), 'projector_rank': rank})
    if len(sectors) != 15 or sum(row['projector_rank'] for row in sectors) != 1653:
        raise ArithmeticError('full capped bosonic spectrum is incomplete')
    # Distinct sector energies for this finite empirical packet; retain exact labels.
    if min(abs(x['energy_GeV']['numeric']-y['energy_GeV']['numeric']) for i,x in enumerate(sectors) for y in sectors[i+1:]) < 1e-7:
        raise ArithmeticError('coincident boson sectors require combined projectors')
    pointers = []
    for group in b['groups']:
        category_name = category[group['modes'][0]]
        occupied = {label: int(label == category_name) for label in labels}
        for create in (True, False):
            pointers.append({'species': group['species'], 'orientation': group['orientation'],
                'momentum_GeV': group['momentum_GeV'], 'pointer': 'create' if create else 'remove',
                'modes': group['modes'], 'reference_electric_charge': str(-s.sympify(group['electric_charge']) if create else s.sympify(group['electric_charge'])),
                'target_occupation': occupied if create else None,
                'success_probability': s.Rational(len(group['modes']), 58) if create else s.S.Zero})
    pointers.append({'species': 'all', 'orientation': 'blank', 'momentum_GeV': None,
        'pointer': 'no-operation', 'modes': [], 'reference_electric_charge': '0',
        'target_occupation': dict.fromkeys(labels, 0), 'success_probability': s.Rational(1, 29)})
    records = []
    for pointer_index, pointer in enumerate(pointers):
        for sector in sectors:
            probability = pointer['success_probability'] if sector['occupation_by_energy_class'] == pointer['target_occupation'] else s.S.Zero
            records.append({'pointer_index': pointer_index, 'species': pointer['species'],
                'orientation': pointer['orientation'], 'momentum_GeV': pointer['momentum_GeV'],
                'pointer': pointer['pointer'], 'energy_projector_index': sector['index'],
                'source_energy_GeV': sector['energy_GeV'],
                'reference_energy_GeV': int(pointer['pointer'] != 'no-operation'),
                'probability': scalar(probability)})
    if sum(s.sympify(row['probability']['exact']) for row in records) != 1:
        raise ArithmeticError('boson pointer probabilities do not close')
    return {'source_Hilbert_dimension': 1653, 'source_oscillator_modes': 56,
        'occupation_cutoff': {'kind': 'total particle number', 'maximum': 2,
            'reason': 'finite closed source domain invariant under mode rotations within equal-energy and charge sectors'},
        'source_preparation': 'capped oscillator vacuum',
        'Hamiltonian': 'H_B=sum_j omega_j n_j on total occupation <=2, vacuum energy subtracted',
        'complete_energy_projectors': sectors,
        'field_truncation': 'a_j=P_{N<=2} b_j P_{N<=2}; vacuum and one-particle propagation are exact, top-occupation CCR differs',
        'vacuum_field_response': 'q_j=(a_j+a_j^dagger)/sqrt(2 omega_j); weight 1/(2 omega_j); 2 omega_j times weight=1',
        'complete_instrument': 'K_{g,+}=sum_{j in g} a_j^dagger tensor |dual j>/sqrt(58), K_{g,-}=sum_{j in g} a_j tensor |j>/sqrt(58), K_empty=sqrt(I-sum_j{a_j,a_j^dagger}/58)',
        'sum_anticommutators_by_total_occupation': {'0': 56, '1': 58, '2': 2},
        'no_operation_squared_by_total_occupation': {'0': '1/29', '1': '0', '2': '28/29'},
        'all_column_cap_identity_error': cap_error, 'all_column_completeness_error': complete_error,
        'all_column_count': len(states), 'reference_dimension': 113,
        'reference_blank_energy_GeV': 0, 'reference_nonblank_energy_GeV': 1,
        'pointer_record_count': len(pointers), 'joint_record_count': len(records),
        'positive_record_count': sum(row['probability']['numeric'] > 0 for row in records),
        'zero_record_count': sum(row['probability']['numeric'] == 0 for row in records),
        'records': records, 'probability_sum': '1'}


@lru_cache(maxsize=1)
def neutrino_flavour_records():
    """Mass-mode flavour analysers, without a weak production-rate claim."""
    u = matter.pmns_matrix()
    if not matter._zero(u.H*u-s.eye(3)):
        raise ArithmeticError('the three-flavour analyzers are not canonical')
    masses = fermions()['masses']['neutrino']
    weights = [[s.simplify(s.expand_complex(u[a,i]*s.conjugate(u[a,i]))) for i in range(3)] for a in range(3)]
    rows = []
    for flavour, name in enumerate(('electron', 'muon', 'tau')):
        for orientation in ('neutrino', 'antineutrino'):
            for momentum in MOMENTA:
                for create in (True, False):
                    for eigenstate, mass in enumerate(masses):
                        weight = weights[flavour][eigenstate]
                        probability = weight/12 if create else s.S.Zero
                        rows.append({'flavour': name, 'orientation': orientation,
                            'momentum_GeV': list(momentum), 'pointer': 'create' if create else 'remove',
                            'mass_eigenstate': eigenstate+1, 'helicity_twice': -1 if orientation == 'neutrino' else 1,
                            'source_energy_GeV': scalar(s.sqrt(1+mass*mass)),
                            'probability': scalar(probability)})
    total = s.simplify(sum(s.sympify(row['probability']['exact']) for row in rows))
    if total != 1:
        raise ArithmeticError('neutrino flavour instrument probabilities do not close')
    return {'source_modes': 12, 'added_source_modes': 0, 'mass_eigenfields': 3,
        'canonical_flavour_map': 'f_alpha,nu=sum_i U_alpha,i a_i,h=-; f_alpha,antinu=sum_i conjugate(U_alpha,i) a_i,h=+',
        'complete_instrument': 'f_alpha,orientation,p^dagger tensor |dual label>/sqrt(12) and its annihilation partner; sum of all effects=I',
        'reference_dimension': 25, 'reference_nonblank_energy_GeV': 1,
        'pointer_record_count': 24, 'joint_nonoutside_record_count': len(rows),
        'records': rows, 'outside_energy_projector_probability_per_pointer': '0',
        'total_record_count_including_zero_outside': len(rows)+24,
        'probability_sum': '1', 'anti_tau_neutrino_analyzer_included': True,
        'Majorana_antineutrino_adds_sterile_species': False,
        'scope': 'normalized opposite-helicity flavour analysers of the three mass eigenfields; no weak production cross-section or helicity suppression rate is assigned'}
