"""Exact protected-source work identities and repair-record sufficiency."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from model.defect_field_bridge import full_domain
from model.certified_reference import source as ports
from model.native_ordered_ground import source
from model.native_spin_observer import algebra
from model.protected_defect_catalog import charge
from model.spin_modular_geometry import regional


def expectation(matrix, vector, denominator=1):
    vector = tuple(map(int, vector))
    total = sum(vector[i]*sum(int(matrix[i, j])*vector[j] for j in range(len(vector)))
                for i in range(len(vector)))
    return F(total, denominator*sum(x*x for x in vector))


def physical_terms():
    result = {}
    edges, odd = source.cube()
    for x, y in edges:
        for axis in 'XYZ':
            word = ['I']*8
            word[x] = word[y] = axis
            result[''.join(word)] = F(1, 4)
    for x in range(8):
        word = ['I']*8
        word[x] = 'X'
        result[''.join(word)] = F(1 if x in odd else -1, 8)
    constant = F(-len(edges), 4)
    matrix = float(constant)*np.eye(256, dtype=complex)
    for word, coefficient in result.items():
        matrix += float(coefficient)*algebra.tensor(word)
    if not np.array_equal(matrix, regional.state()['energy_integer']/8):
        raise ArithmeticError('energy-score dictionary differs from the actual source clock')
    return constant, dict(sorted(result.items()))


def repair_blindness():
    checked = 0
    edges, _ = source.cube()
    for word in range(256):
        mask = sum(1 << (2*x+(0 if word >> x & 1 else 1)) for x in range(8))
        for seam in edges:
            for law in ('C', 'E'):
                for coin in (0, 1):
                    record, out, sign, current = full_domain.branch_action(
                        mask, law, coin, sites=8, flavors=2, seam=seam)
                    if (record, out, sign, current) != (0, mask, 1, 0):
                        raise ArithmeticError('protected repair-record blindness failed')
                    checked += 1
    h, scale = source.matrix()
    neel = sum(1 << x for x in source.cube()[1])
    if (int(h[0, 0]), int(h[neel, neel]), scale) != (0, -48, 8):
        raise ArithmeticError('exact equal-record/different-energy control changed')
    return {'complete_protected_basis_branches_checked': checked,
            'each_fine_branch_on_protected_code': 'identity/sqrt(2) with raw record0',
            'arbitrary_record_adaptive_C_E_schedule_energy_blind': True,
            'bare_repair_system_energy_change_exact': '0',
            'two_protected_product_input_energies': ['0', '-6'],
            'both_inputs_have_identical_complete_repair_history_laws': True,
            'repair_record_cost_recovers_noncentral_H_B_energy': False}


def primitive_action():
    """Signed integrated action from all twelve exact six-level source ports."""
    matrices = ports.specification()
    rows = []
    for kind, angle in (('Z', 0.25), ('Z', -0.25), ('X', 0.25), ('I', 0.25)):
        controls = ports.literal_controls(kind, angle)
        value = ports.ZERO
        for coefficient, matrix in zip(controls, matrices, strict=True):
            # iG has real diagonal -Im(G); the input is the full six-level |0>.
            value = ports.plus(value, ports.scale(matrix[0][0][1], -coefficient))
        enclosure = ports.field_interval(value)
        rows.append({'kind': kind, 'angle': str(F(angle)),
                     'twelve_literal_controls': [str(x) for x in controls],
                     'signed_action_Q_sqrt5': [str(x) for x in value],
                     'signed_action_interval': enclosure.json(),
                     'full_generator_error': str(ports.generator_error(kind, F(angle), controls)),
                     'full_source_levels': 6, 'input': '|0><0|',
                     'pulse_count': 1})
    if ports.literal_controls('Z', -0.25) != tuple(-x for x in ports.literal_controls('Z', 0.25)):
        raise ArithmeticError('literal signed-source inverse law changed')
    return {'primitive_records': rows,
            'signed_action_definition': 'sum_surviving_prefixes Tr(rho_prefix*i*G_prefix)',
            'G_prefix_uses_all_twelve_source_ports': True,
            'signed_action_is_nonnegative_repair_count': False,
            'constant_phase_may_be_discarded_in_controlled_clock': False}


@lru_cache(maxsize=1)
def work():
    op = charge.operators()
    v, epsilon, certificate = charge.ground_input()
    h = op['H']
    norm = sum(int(x)**2 for x in v)
    rows = []
    for q, a in sorted(op['twice_creation'].items()):
        p = a.T@a
        b = 2*a.T@h@a-p@h-h@p
        if not np.array_equal(b.T, b):
            raise ArithmeticError('source work observable is not Hermitian')
        weight = expectation(p, v, 4)
        energy = expectation(b, v, 64)
        bnorm = F(max(sum(abs(int(x)) for x in row) for row in b), 64)
        wlo, whi = weight-2*epsilon, weight+2*epsilon
        elo, ehi = energy-2*bnorm*epsilon, energy+2*bnorm*epsilon
        if min(wlo, elo) <= 0:
            raise ArithmeticError('created charge channel lost positive weight or energy')
        rows.append({'charge': q, 'probability_exact_proposal': weight,
                     'probability_lower': wlo, 'probability_upper': whi,
                     'work_observable_norm_upper': bnorm,
                     'unnormalized_excitation_exact_proposal': energy,
                     'unnormalized_excitation_lower': elo,
                     'unnormalized_excitation_upper': ehi,
                     'mean_excitation_lower': elo/whi,
                     'mean_excitation_upper': ehi/wlo})
    z = op['local_Z']
    change = z@h@z-h
    unitary_work = expectation(change, v, 8)
    norm_change = F(max(sum(abs(int(x)) for x in row) for row in change), 8)
    lower = unitary_work-2*norm_change*epsilon
    if lower <= 1:
        raise ArithmeticError('energy-blind repair witness lost its energy separation')
    # Exact ground charge0 and [H,Q]=0 make the two creation channels energy-orthogonal.
    difference = abs(sum(r['unnormalized_excitation_exact_proposal'] for r in rows)-unitary_work)
    if difference > 8*14*epsilon:
        raise ArithmeticError('charge-resolved work fails the same-ground unitary sum rule')
    data = regional.state()
    if not np.array_equal(data['energy_integer'],
                          np.array([(-1)**sum((mask >> j) & 1 for j in source.cube()[1])
                                    for mask in range(256)])[:, None]*h*
                          np.array([(-1)**sum((mask >> j) & 1 for j in source.cube()[1])
                                    for mask in range(256)])[None, :]):
        raise ArithmeticError('modular and energy references are different source actions')
    k = regional.logarithm()
    rz = tuple(tuple(data['rho'][i][j]*((-1)**((i&1)+(j&1))) for j in range(4)) for i in range(4))
    modular = sum((rz[i][j]-data['rho'][i][j])*k['matrix'][j][i]
                  for i in range(4) for j in range(4))
    modular_error = 2*k['operator_error_upper']+16*epsilon
    if modular <= modular_error:
        raise ArithmeticError('actual regional modular work lost positivity')
    return {'channels': rows, 'complete_ground_certificate': certificate,
            'ground_vector_error': epsilon,
            'unitary_Z_excitation_lower': lower,
            'unitary_Z_excitation_upper': unitary_work+2*norm_change*epsilon,
            'unitary_Z_regional_entropy_change_exact': F(0),
            'unitary_Z_modular_relative_entropy_lower': modular-modular_error,
            'unitary_Z_modular_relative_entropy_upper': modular+modular_error,
            'charge_channel_energy_sum_equals_unitary_creation_energy_on_exact_ground': True,
            'zero_repair_record_effort_implies_zero_spectral_energy': False,
            'Gibbs_reference_or_temperature_assumed': False}
