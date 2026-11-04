"""Exact reversal of all native proposal and thermal collision outcomes.

The reversal pairs complete recorded branches of the same stochastic law.
It does not reverse the implementation pulses or erase a stored tape.
"""
import hashlib
import json

import sympy as s

from .bounds import source_data
from ..native_equilibration import reciprocal_vacuum as source


def complete_labels(law='C'):
    levels, _, _, _ = source_data()
    native, proposals = source.proposal_records(law)
    rows = []
    for i in range(16):
        for address, (geometry, field) in enumerate(proposals[i//8]):
            j = int(8*geometry+field)
            allowed = i//8 != j//8
            if allowed:
                if levels[i].hi < levels[j].lo:
                    low, high = i, j
                elif levels[j].hi < levels[i].lo:
                    low, high = j, i
                else:
                    raise ArithmeticError('a cross-geometry bath orientation is unresolved')
            else:
                low = high = i
            for before in (0, 1):
                accepted = allowed and before == int(i == low)
                final = j if accepted else i
                after = 1-before if accepted else before
                row = (i, address, int(geometry), int(field), before, after, final, int(accepted))
                # Exact coefficients of beta E in Delta E_system+Delta E_bath.
                energy = [0]*16
                energy[final] += 1
                energy[i] -= 1
                energy[high] += after-before
                energy[low] -= after-before
                if any(energy):
                    raise ArithmeticError('a retained source record violates exact energy balance')
                rows.append(row)
    lookup = {row: index for index, row in enumerate(rows)}
    if len(lookup) != 512:
        raise ArithmeticError('complete proposal/bath records were omitted')
    inverse_proposal = {
        (initial, int(g), int(m)): address
        for initial in (0, 1)
        for address, (g, m) in enumerate(proposals[initial])}
    reverse = []
    for row in rows:
        i, address, geometry, field, before, after, final, accepted = row
        if accepted:
            inverse = (final, inverse_proposal[final//8, i//8, i % 8],
                       i//8, i % 8, after, before, i, 1)
        else:
            inverse = row
        reverse.append(lookup[inverse])
    if any(reverse[reverse[k]] != k for k in range(512)):
        raise ArithmeticError('record reversal is not an involution')
    return tuple(rows), tuple(reverse), native


def certificate():
    low, high = s.symbols('w_low w_high', positive=True)
    upward, downward = high/(16*(low+high)), low/(16*(low+high))
    if s.cancel(upward/downward-high/low) != 0:
        raise ArithmeticError('the exact thermal record likelihood ratio failed')
    laws = {}
    for law in ('C', 'E'):
        rows, reverse, native = complete_labels(law)
        payload = {'rows': rows, 'reverse': reverse, 'native_proposal_tape': native.tolist()}
        laws[law] = {
            'complete_records': len(rows),
            'accepted_records': sum(row[-1] for row in rows),
            'rejected_and_noop_records': sum(not row[-1] for row in rows),
            'reversal_fixed_records': sum(k == reverse[k] for k in range(len(rows))),
            'full_label_involution_sha256': hashlib.sha256(
                json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
            'every_record_energy_balance_exact': True,
            'all_native_proposal_tapes_retained': True,
        }
    return {
        'laws': laws,
        'branch_log_likelihood_ratio': 'log(kappa_r/kappa_reverse(r))=-beta Q_r',
        'path_entropy_production': 'log(q0(i0)/qn(in))-beta sum Q_r',
        'mean_path_entropy_production': 'D(q0||pi)-D(T^n q0||pi)>=0',
        'faithful_initial_populations_integral_fluctuation_identity': 'E_forward[exp(-Sigma)]=1',
        'quantum_dephasing_identity': 'D(sigma||rho)=D(sigma||Delta_E(sigma))+D(q||pi)',
        'coherent_input_entropy_heat_balance':
            'Delta S-beta Q=D(sigma||Delta_E(sigma))+D(q||pi)-D(T^n q||pi)',
        'record_entropy_is_discarded_or_memory_erased': False,
        'bath_log_likelihood_flow_is_actual_bath_von_Neumann_entropy_change': False,
        'pulse_controller_work_and_fresh_bath_preparation_derived': False,
    }
