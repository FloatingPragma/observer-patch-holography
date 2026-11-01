"""Same-carrier energy continuity, neutral readout and source work receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ..defect_field_benchmark.benchmark import compare_receipt
from ..defect_field_benchmark.quasiparticles import quasiparticles
from ..defect_field_bridge import joint
from ..defect_field_propagation import propagation
from ..observer.source_scalar import compiler
from . import balance, energy, readout

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def continuity_certificate():
    rows = []
    total_modes, total_scalar = 0, 0.
    for indices in propagation.components():
        for phase in (0., math.pi/2):
            data = energy.component_energy(indices, phase, reference='common_static')
            incoming = [np.zeros_like(data['B']) for _ in range(12)]
            supports = [set(np.flatnonzero(np.any(abs(B) > 0, axis=0))) for B in data['local_B']]
            current_count = 0
            for source in range(12):
                for destination in range(source+1, 12):
                    distance = min(destination-source, 12-(destination-source))
                    if distance > 2:
                        if supports[source] & supports[destination]:
                            raise ArithmeticError('nominally disjoint local energies share support')
                        continue
                    current = energy.energy_current(data, source, destination)
                    if np.linalg.norm(current-current.conj().T) > energy.OPERATOR_GUARD:
                        raise ArithmeticError('energy current is not Hermitian')
                    incoming[destination] += current
                    incoming[source] -= current
                    current_count += 2
            for site, local in enumerate(data['local_B']):
                direct = 1j*(data['B']@local-local@data['B'])
                if np.linalg.norm(direct-incoming[site]) > energy.OPERATOR_GUARD:
                    raise ArithmeticError('same-ring all-sector local continuity identity failed')
            if np.linalg.norm(sum(incoming)) > energy.OPERATOR_GUARD:
                raise ArithmeticError('internal energy currents did not cancel globally')
            rows.append({'internal_component': list(indices), 'phase': phase,
                         'spatial_CAR_modes': data['modes'], 'Nambu_directions': 2*data['modes'],
                         'local_energies': 12, 'oriented_range_at_most_two_current_slots': current_count,
                         'all_other_current_slots_zero_by_disjoint_support': True,
                         'identity_absolute_matrix_guard': energy.OPERATOR_GUARD})
            if phase == 0:
                total_modes += data['modes']
                total_scalar += data['reference_scalar']
    expected_scalar = -11*quasiparticles()['vacuum_energy']
    if total_modes != 1080 or abs(total_scalar-expected_scalar) > 2e-9:
        raise ArithmeticError('energy decomposition lost spatial modes or the common reference scalar')
    # A scalar cannot be inferred from a Majorana/Gaussian map. Execute it.
    word = [compiler.pulse('I', -total_scalar, [0])]
    actual = compiler.execute(word, 1)
    if np.linalg.norm(actual-np.exp(-1j*total_scalar)*np.eye(2), 2) > readout.SOURCE_GUARD:
        raise ArithmeticError('explicit common-reference source phase failed')
    return {'component_rows': rows, 'total_CAR_modes': total_modes,
            'complete_centered_quadratic_identity': '[F(A),F(B)]=F([A,B]), F(B)=Psi† B Psi/2',
            'normal_order_scalar': 'F(B)+c=c† h c+(c† Delta c†+h.c.)/2+c-Tr(h)/2',
            'local_energy': 'all onsite Higgs pairing plus one half of each incident kinetic bond and the stated scalar share',
            'current_convention': 'J[y,x]=F(i[B_y,B_x]), current from y to x',
            'work_power': 'dot(phi) F(partial_phi B_x); fixed reference has zero derivative',
            'global_balance': 'd<H>/dt_eff=dot(phi)<partial_phi H> under the simultaneous Hamiltonian',
            'actual_ordered_product_equals_simultaneous_exponential': False,
            'bare_propagation_reference_scalar': 0,
            'optional_common_static_reference_scalar': total_scalar,
            'common_reference_phase_source_resources': compiler.resources(word)}, actual


def prepare_probe(kind):
    if kind not in ('current', 'repair_sensitive'):
        raise ValueError('one specified neutral energy probe required')
    second = 4 if kind == 'current' else 7
    target = np.zeros(64)
    target[1+32] = target[second+32] = 1/math.sqrt(2)
    word = compiler.compile_real_state(target)
    if kind == 'current':
        word += [compiler.pulse('Z', -math.pi/4, [2]), compiler.pulse('I', math.pi/4, [2])]
        target = target.astype(complex); target[second+32] *= 1j
    state = np.zeros(64, complex); state[0] = 1
    actual = compiler.execute(word, 6, state)
    if np.linalg.norm(actual-target) > readout.SOURCE_GUARD:
        raise ArithmeticError('complete neutral source probe preparation failed')
    matrix = actual.reshape(4, 16).T
    spectator_Q = np.diag([0., 1., -1., 0.])
    if np.linalg.norm(joint.charge()@matrix+matrix@spectator_Q) > readout.SOURCE_GUARD:
        raise ArithmeticError('active probe and retained spectator are not neutral')
    return matrix, compiler.resources(word)


def build():
    continuity, reference_phase = continuity_certificate()
    oracle = energy.small_energy_oracle()
    current = energy.ring_current_oracle()
    observables = {'H': oracle['H'], 'h0': oracle['h_sites'][0], 'h1': oracle['h_sites'][1],
                   'J0to1': oracle['J'][0, 1], 'drive_power_per_phase_rate': oracle['dH_dphi'],
                   'ring_J0to2': current['operator']}
    arrays = {'common_reference_scalar_phase': reference_phase, **observables}
    definitions = readout.exact_quadratic_terms()
    pointers, pointer_rows = {}, []
    for definition in definitions.values():
        for term in definition['terms']:
            key = (term['kind'], tuple(term['modes']), tuple(term['coefficient']))
            if key in pointers:
                continue
            pointer = readout.grouped_pointer(*key)
            tag = f'pointer{len(pointers)}'
            pointers[key] = (tag, pointer)
            arrays[tag] = pointer['unitary']
            pointer_rows.append({'tag': tag, **term, 'eigenvalues': pointer['eigenvalues'],
                                 'full_source_input_columns': 64, 'outcomes': 4,
                                 'each_field_Kraus_commutes_with_Q': True,
                                 'source_resources': compiler.resources(pointer['word'])})
    states, preparations = {}, {}
    for name in ('current', 'repair_sensitive'):
        states[name], preparations[name] = prepare_probe(name)
    repair_energy = []
    for policy in ('coherent', 'flagged'):
        for law in ('C', 'E'):
            K, _, _ = joint.repair(law, policy)
            name = policy+'_'+law
            before = states['repair_sensitive']
            after = np.concatenate([k@before for k in K], axis=1)
            states[name] = after
            defect = sum(k.conj().T@oracle['H']@k for k in K)-oracle['H']
            repair_energy.append({'policy': policy, 'law': law,
                                  'energy_transfer_on_sensitive_source_state': balance.expectation(oracle['H'], after)-balance.expectation(oracle['H'], before),
                                  'adjoint_energy_defect_norm': float(np.linalg.norm(defect, 2)),
                                  'complete_native_records': 32,
                                  'native_records_retained_as_reference_columns': after.shape[1]})
    measurements = []
    for name, state in states.items():
        arrays['state_'+name] = state
        observed = {}
        for key, (tag, pointer) in pointers.items():
            outputs = np.stack([K@state for K in pointer['K']])
            probabilities = np.sum(abs(outputs)**2, axis=(1, 2))
            observed[key] = float(probabilities@pointer['eigenvalues'])
            arrays[name+'_'+tag+'_probabilities'] = probabilities
        values = {}
        for label, definition in definitions.items():
            estimate = definition['constant']*float(np.linalg.norm(state)**2)
            bound = 0.
            for term in definition['terms']:
                key = (term['kind'], tuple(term['modes']), tuple(term['coefficient']))
                estimate += observed[key]
                bound += 4*readout.SOURCE_GUARD*np.linalg.norm(pointers[key][1]['operator'], 2)
            direct = balance.expectation(observables[label], state)
            if abs(estimate-direct) > max(bound, 2e-10):
                raise ArithmeticError('charge-preserving source energy readout differs from the full observable')
            values[label] = {'source_readout': estimate, 'direct_full_operator': direct,
                             'numerical_readout_error_bound': float(max(bound, 2e-10))}
        measurements.append({'state': name, 'reference_columns': state.shape[1], 'observables': values})
    protocols = []
    for policy in ('coherent', 'flagged'):
        for law in ('C', 'E'):
            result = balance.protocol_balance(policy, law)
            protocols.append(result['report'])
    files = [HERE/'energy.py', HERE/'balance.py', HERE/'readout.py', HERE/'packet.py',
             HERE/'test_energy.py', HERE/'README.md', ROOT/'model/defect_field_bridge/joint.py',
             ROOT/'model/defect_field_bridge/full_domain.py', ROOT/'model/defect_field_bridge/coherent_extension.py',
             ROOT/'model/defect_field_interactions/interaction.py', ROOT/'model/defect_field_propagation/propagation.py',
             ROOT/'model/defect_field_benchmark/quasiparticles.py', ROOT/'model/common_source/calibrated_matter.py',
             ROOT/'model/common_source/calibrated_matter_parameters.json', ROOT/'model/defect_field_clock/clock.py',
             ROOT/'model/observer/source_scalar/compiler.py', ROOT/'model/common_source/chiral_candidate/gaussian_source.py']
    report = {'schema': 'oph.defect_field_energy.same_carrier_balance.v1',
              'source_sha256': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
              'continuity': continuity, 'grouped_neutral_pointer_definitions': pointer_rows,
              'exact_quadratic_observable_decomposition': definitions,
              'separate_preparations_required_for_noncommuting_readout_groups': True,
              'initial_neutral_probe_source_resources': preparations,
              'source_measurements': measurements, 'repair_energy_transfers': repair_energy,
              'actual_ordered_control_balances': protocols,
              'small_two_site_mutual_current': 'exactly zero for this allocation; no nonzero flow claimed',
              'nonzero_ring_current': 'J_0to2 on the actual12-ring, its complete four-endpoint-mode CAR operator, not a projected state',
              'kinetic_rest_consistency': 'Same supplied ring symbol E(k)^2=sin(k)^2/a^2+E0^2: E_second(k=0)=1/(a^2 E0); for p=k/a and c=1, E_second(p=0)=1/E0. This is the analytic symbol, not a continuum limit inferred from12 samples.',
              'physical_scope': ['all-sector finite energy balance for the specified1080-mode ring and Higgs action',
                                 'explicit switching work for a declared interpolation of the executed ordered factors',
                                 'complete repair energy transfer, without assuming a hidden thermodynamic controller',
                                 'charge-preserving grouped readout on the selected finite source factor'],
              'not_claimed': ['independent mass prediction', 'a3+1-dimensional conserved stress tensor',
                              'Einstein dynamics or a continuum limit', 'physical detector shots or hardware precision']}
    return report, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    report, arrays = build()
    output, archive = HERE/'results.json', HERE/'results.npz'
    if args.write:
        output.write_text(json.dumps(report, indent=2)+'\n')
        np.savez_compressed(archive, **arrays)
    elif args.check:
        compare_receipt(report, json.loads(output.read_text()))
        with np.load(archive, allow_pickle=False) as old:
            if set(old.files) != set(arrays):
                raise ValueError('complete energy-source array domain changed')
            for key, value in arrays.items():
                if value.shape != old[key].shape or np.max(abs(value-old[key])) > 2e-7:
                    raise ValueError('complete energy-source output changed: '+key)
    print('same-carrier energy balance and neutral source readouts verified')


if __name__ == '__main__':
    main()
