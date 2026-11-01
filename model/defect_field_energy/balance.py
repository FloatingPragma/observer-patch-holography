"""Actual ordered source controls and complete repair energy accounting."""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import expm

from ..defect_field_bridge import joint
from ..defect_field_clock import clock
from ..defect_field_interactions import interaction
from ..observer.source_scalar import compiler
from .energy import small_energy_oracle

SOURCE_GUARD = 2e-7


def expectation(operator, state):
    return float(np.trace(state.conj().T@operator@state).real)


def source_segments(current=0, phase=math.pi/2, duration=1.):
    effective = float(duration*clock.drive(current))
    chart = interaction.bare_component_chart((27, 78))
    rword = compiler.remap(chart['source_word'], [0, 1])+compiler.remap(chart['source_word'], [2, 3])
    hop_word = (rword+joint.routed_hopping_word(.5j, effective, 0, 2)
                +joint.routed_hopping_word(-.5j, effective, 1, 3)+compiler.inverse(rword))
    pair = interaction.pair_word(chart['energy'], phase, effective)
    pair_word = compiler.remap(pair, [0, 1])+compiler.remap(pair, [2, 3])
    hopping, pairing = compiler.execute(hop_word, 4), compiler.execute(pair_word, 4)
    target = small_energy_oracle(phase)
    for unitary, H in ((hopping, target['H_hopping']), (pairing, target['H_onsite'])):
        if np.linalg.norm(unitary-expm(-1j*effective*H), 2) > SOURCE_GUARD:
            raise ArithmeticError('executed energy-control segment disagrees with its generator')
    return {'unitaries': (hopping, pairing), 'generators': (2*target['H_hopping'], 2*target['H_onsite']),
            'segment_durations': (effective/2, effective/2), 'effective_duration': effective,
            'source_words': (hop_word, pair_word),
            'control_interpolation': 'two equal proper-time subintervals; twice each selected generator reproduces the declared ordered factors; clock advances once'}


def protocol_balance(policy='coherent', law='C', current=0):
    reference = small_energy_oracle(0.)
    H = reference['H']
    initial, prep = joint.preparation()
    controls = source_segments(current)
    u1, u2 = controls['unitaries']
    k1, k2 = controls['generators']
    first, second = u1@initial, u2@u1@initial
    work = (expectation(k1-H, initial), expectation(k2-k1, first), expectation(H-k2, second))
    drift = (expectation(k1, first)-expectation(k1, initial),
             expectation(k2, second)-expectation(k2, first))
    energy_change = expectation(H, second)-expectation(H, initial)
    if abs(energy_change-sum(work)-sum(drift)) > 2e-10:
        raise ArithmeticError('source control switching-work identity failed')
    K, _, repair_resources = joint.repair(law, policy)
    records = [k@second for k in K]
    repaired = np.concatenate(records, axis=1)
    local_records = []
    for index, state in enumerate(records):
        local_records.append({'coin': index//16, 'raw': index % 16,
                              'probability': float(np.linalg.norm(state)**2),
                              'energy_contribution': expectation(H, state),
                              'local_energy_contributions': [expectation(h, state) for h in reference['h_sites']]})
    transfer = expectation(H, repaired)-expectation(H, second)
    total_change = expectation(H, repaired)-expectation(H, initial)
    if abs(total_change-energy_change-transfer) > 2e-10:
        raise ArithmeticError('complete drive-plus-repair energy balance failed')
    return {'states': {'initial': initial, 'after_hopping': first, 'after_drive': second, 'after_repair': repaired},
            'report': {'policy': policy, 'law': law, 'current': current,
                       'effective_duration': controls['effective_duration'],
                       'control_interpolation': controls['control_interpolation'],
                       'switching_work': list(work), 'segment_energy_drift': list(drift),
                       'reference_energy_change_before_repair': energy_change,
                       'complete_repair_energy_transfer': transfer,
                       'total_reference_energy_change': total_change,
                       'complete_native_records': local_records,
                       'repair_transfer_assigned_to_heat_or_controller_work': False,
                       'preparation_resources': prep,
                       'segment_resources': [compiler.resources(word) for word in controls['source_words']],
                       'repair_resources': repair_resources}}
