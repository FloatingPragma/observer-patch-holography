"""Full-occupation source certificate for a driven spatial CAR action."""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ..defect_field_benchmark.quasiparticles import quasiparticles
from . import interaction as model

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
ARRAYS = HERE/'results.npz'


def build():
    chart = model.paired_chart()
    source = quasiparticles()
    frozen = json.loads((ROOT/'model/defect_field_bridge/dynamics_controls.json').read_text())
    programmed = np.array(list(map(float, frozen['energies_GeV'])))
    if programmed.shape != (90,) or np.max(abs(programmed-source['energies'])) > 1e-8:
        raise ArithmeticError('interaction program lost the same frozen static spectrum')
    groups = {}
    for pair in chart['pairs']:
        energy = float(programmed[pair['left']])
        key = (pair['left_charge6'], energy)
        groups.setdefault(key, []).append(pair)
    rows, arrays = [], {}
    bare_join = model.bare_component_chart((27, 78))
    arrays['complete_bare_QP_chart_lepton27_78'] = bare_join['source_unitary']
    swap = model.compiler.execute(model.fermionic_swap_word(), 2)
    target_swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, -1]], complex)
    if np.linalg.norm(swap-target_swap, 2) > model.SOURCE_GUARD:
        raise ArithmeticError('complete fermionic routing swap failed')
    arrays['complete_fermionic_swap'] = swap
    routing_swaps = 24*sum(pair['right']-pair['left']-1 for pair in chart['pairs'])
    for number, ((charge6, energy), pairs) in enumerate(sorted(groups.items())):
        for current in (-1, 0, 1):
            duration = 1+current/4
            gate, word = model.pair_gate(energy, math.pi/2, duration)
            counts = model.compiler.resources(word)
            probabilities = []
            for initial in range(4):
                observation = model.observed_pair(energy, math.pi/2, duration, initial)
                arrays[f'group{number}_current{current}_initial{initial}'] = observation['amplitudes']
                probabilities.append(observation['probabilities'].tolist())
            arrays[f'group{number}_current{current}_gate'] = gate
            probability = math.sin(energy*duration)**2
            rows.append({'energy_GeV': energy, 'negative_charge_sixths': charge6,
                         'pairs_per_site': len(pairs), 'sites': 12, 'current': current,
                         'effective_duration': duration, 'all_four_input_columns_executed': True,
                         'occupation_probabilities_by_input': probabilities,
                         'ideal_vacuum_to_pair_probability': probability,
                         'ideal_work_after_return_to_static_action_GeV': 2*energy*probability,
                         'source_operator_error_guard': model.SOURCE_GUARD,
                         'source_resources': counts})
    # The common reference subtraction is implemented as an actual scalar
    # source phase, once for the whole spatial state, not once per block.
    vacuum = float(frozen['vacuum_energy_GeV'])
    scalar_rows = []
    for current in (-1, 0, 1):
        duration = 1+current/4
        word = [model.compiler.pulse('I', 11*vacuum*duration, [0])]
        gate = model.compiler.execute(word, 1)
        expected = np.exp(1j*11*vacuum*duration)*np.eye(2)
        if np.linalg.norm(gate-expected, 2) > model.SOURCE_GUARD:
            raise ArithmeticError('spatial vacuum reference source phase changed')
        arrays[f'global_reference_current{current}'] = gate
        scalar_rows.append({'current': current, 'effective_duration': duration,
                            'source_resources': model.compiler.resources(word)})
    energy = float(programmed[source['ramsey_indices'][1]])
    echoes = model.repair_echo(energy)
    echo_rows = {}
    for law, echo in echoes.items():
        arrays[f'complete_{law}_repair_echo'] = echo.pop('all_record_field_amplitudes')
        echo_rows[law] = echo
    number_laws = {str(j): model.extended_number_law(1., math.pi/2, j) for j in (-1, 0, 1)}
    files = [HERE/'interaction.py', HERE/'packet.py', HERE/'test_interaction.py', HERE/'README.md',
             ROOT/'model/defect_field_benchmark/quasiparticles.py',
             ROOT/'model/defect_field_bridge/dynamics_controls.json',
             ROOT/'model/defect_field_bridge/full_domain.py',
             ROOT/'model/common_source/calibrated_matter.py',
             ROOT/'model/common_source/calibrated_matter_parameters.json',
             ROOT/'model/common_source/chiral_candidate/action.py',
             ROOT/'model/observer/source_scalar/compiler.py']
    result = {
        'schema': 'oph.defect_field_interactions.full_occupation_higgs_drive.v1',
        'source_sha256': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        'design': 'external neutral Higgs phase drive H=(0,v exp(i phi)), fixed gauge background, phi=pi/2',
        'creation_phase_signs_from_actual_invariant_flags': {n: row['creation_phase_sign'] for n, row in model.channel_pairings().items()},
        'domain': {'ports': 12, 'CAR_modes': 1080, 'all_occupation_sectors': 1081,
                   'neutral_two_mode_blocks': 504, 'massless_spectator_modes': 72,
                   'local_complete_Fock_states_per_block': 4,
                   'all_global_Fock_states_retained': '2^1080 by the proved complete tensor product; never materialized',
                   'same_port_collisions_retained': True},
        'action': {'local_pair': 'E cos(phi)(n_a+n_b-1)+E sin(phi)(a^dagger b^dagger+b a)',
                   'same_static_restriction': 'phi=0 gives the same diagonal onsite action; restrict to hard-core N=2 to recover the prior U tensor U phase',
                   'global_reference': 'subtract 11 E_vac I from twelve onsite actions, with the same fixed reference at every phi',
                   'clock': 'the entire Hamiltonian uses tau_eff=tau*(1+j/4)',
                   'charge_and_parity': 'every pair has opposite electric charges; global charge and fermion parity are conserved',
                   'number': 'quasiparticle number changes by two',
                   'drive_switching_work': 'sudden external switch from phi0 to phi and back supplies work; H_phi is conserved during its constant segment',
                   'full_quantum_gauge_interaction': False},
        'explicit_bare_QP_source_join': {
            'original_component': [27, 78], 'creation_phase_sign': bare_join['creation_phase_sign'],
            'energy_GeV': bare_join['energy'], 'all_four_source_columns_retained': True,
            'static_and_driven_action_intertwined': True,
            'commuting_QP_phase_included_to_fix_driven_term': True,
            'source_resources': model.compiler.resources(bare_join['source_word'])},
        'source_program': {'frozen_spectrum_admission_GeV': 1e-8,
                           'distinct_energy_charge_blocks': len(groups),
                           'source_executed_block_and_current_cases': len(rows),
                           'rows': rows, 'global_scalar_phase_rows': scalar_rows,
                           'all_local_source_array_outputs_retained': True,
                           'adjacent_fermionic_swaps_in_full_routed_layer': routing_swaps,
                           'fermionic_swap_source_resources': model.compiler.resources(model.fermionic_swap_word()),
                           'all_1080_qubits_materialized': False,
                           'global_tensor_operator_error_bound_from_local_guards': (1+model.SOURCE_GUARD)**(505+routing_swaps)-1,
                           'physical_hardware_or_interval_certificate': False,
                           'arbitrary_reference_extension': 'complete local input-column bounds extend under tensor identity; all block factors are even CAR operations'},
        'complete_number_laws_for_unrounded_action': number_laws,
        'number_law_precision': '70 decimal digits on factors evaluated in binary floating point; decimal tails retained below double underflow',
        'drive_then_full_domain_repair_echo': echo_rows,
        'full_spatial_projector_echo': {str(j): model.full_repair_echo(current=j) for j in (-1, 0, 1)},
        'echo_initial_state': 'driven seam empty, two distinct occupied spectator ports; local |00> versus |11> means global N2 versus N4',
        'echo_without_repair_return_probability': 1,
        'echo_with_repair_law': '1-2 p(1-p), p=sin(E tau)^2; all32 coin/raw records and16 seam occupations retained',
        'remaining_joins': ['This is a specified externally driven quadratic background, not a quantized Higgs or gauge scattering calculation.',
                            'A nonstatic drive can change the energy basis; static quasiparticle particle labels are not asserted at all times.',
                            'The explicit split-step composition with propagation and repair is a model law, not an unbounded-continuum or observed-particle proof.',
                            'Calibrated Yukawas and the common clock normalization remain stated model inputs.'],
    }
    return result, arrays


def compare(actual, expected, path='receipt'):
    if isinstance(actual, (tuple, list)) and isinstance(expected, (tuple, list)):
        if len(actual) != len(expected):
            raise ValueError(path+': lengths differ')
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, path+f'[{i}]')
    elif isinstance(actual, dict) and isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError(path+': fields differ')
        for key in actual:
            compare(actual[key], expected[key], path+'.'+key)
    elif type(actual) is float and type(expected) is float:
        if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=2e-9, abs_tol=2e-8):
            raise ValueError(path+': numeric output changed')
    elif type(actual) is str and type(expected) is str and ('probabilities_decimal[' in path or path.endswith('probability_staying_N2_decimal')):
        a, b = Decimal(actual), Decimal(expected)
        if (a == 0) != (b == 0) or (a and abs((a-b)/a) > Decimal('1e-8')):
            raise ValueError(path+': small positive number-sector probability changed')
    elif actual != expected or type(actual) is not type(expected):
        raise ValueError(path+': exact output changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result, arrays = build()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
        np.savez_compressed(ARRAYS, **arrays)
    elif args.check:
        compare(result, json.loads(OUTPUT.read_text()))
        with np.load(ARRAYS, allow_pickle=False) as saved:
            if set(saved.files) != set(arrays):
                raise ValueError('complete source-array domain changed')
            for key, value in arrays.items():
                if value.shape != saved[key].shape or np.max(abs(value-saved[key])) > 2e-7:
                    raise ValueError('complete source-array output changed: '+key)
    print('complete occupation-space Higgs drive and repair echo verified')


if __name__ == '__main__':
    main()
