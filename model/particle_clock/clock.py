"""The same native corner, interacting energy pair and complete clock readout.

The supplied field dictionary sends native |1,2> to the zero-contact singlet
r and |7,2> to the onsite eta pair d.  At total momentum (pi,pi,pi), the
periodic Weyl/contact action has H r=0 and H d=U d.  This is a faithful Gibbs
reference on their two-dimensional invariant span, not on the full field.
Every corner rejection, native raw outcome and source-address abort is kept.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..local_repair_family import acquisition, inference, learned_flow, record_clock
from ..local_repair_family import modular_clock as native_clock, local_bias as native
from ..observer.source_scalar import compiler
from ..observer.source_discrimination import circuit
from ..common_source import native_routing
from ..bound_pair.fiber import clock_subspace

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PORTS = (1, 7)
SPECTATOR = 2
POINTER_GATE = F(1, 10**8)
PHASE_CAP = F(3, 160)


def corner_indices():
    states = tuple(native.basis())
    return tuple(states.index(tuple(int(p in (port, SPECTATOR)) for p in range(12)))
                 for port in PORTS)


def _corner_codes():
    return tuple(circuit.encode(native.basis()[i]) for i in corner_indices())


def infer_clock(counts):
    """Record-only estimate; no true weights, reference or field fit is read."""
    result = inference.infer_counts(counts)
    if result['status'] != 'accept':
        return {'status': 'stop', 'reason': result['stop_reason']}
    _, tree, paths, _ = inference._design()
    coefficients = tuple(a-b for a, b in zip(paths[PORTS[0]], paths[PORTS[1]], strict=True))
    distance = sum(map(abs, coefficients))
    if distance != 2:
        raise ArithmeticError('the declared particle-clock tree path changed')
    intervals = []
    for a, b in tree:
        lo, hi = result['probability_intervals'][tuple(sorted((a, b)))]
        odds = (lo/(1-lo), hi/(1-hi))
        intervals.append(odds if a < b else (1/odds[1], 1/odds[0]))
    lo, hi = inference._path_interval(coefficients, tuple(intervals))
    weights = result['relative_weights']
    ratio = weights[PORTS[0]]/weights[PORTS[1]]
    # These are intervals for beta*U. U=1 is the declared action convention.
    # Half the two-eigenvalue oscillation is the centered operator bound.
    potential = F(25, 4)*distance*result['budget']['probability_radius']
    return {'status': 'accept', 'relative_weights': weights,
            'energy_ratio_estimate': str(ratio), 'energy_ratio_interval': [str(lo), str(hi)],
            'beta_times_U_estimate_diagnostic': math.log(ratio.numerator)-math.log(ratio.denominator),
            'beta_times_U_error_on_confidence_event': str(potential),
            'clock_halfdiamond_error_per_unit_parameter': str(potential/2),
            'statistical_failure_upper': str(result['budget']['statistical_failure_upper']),
            'tree_distance': distance, 'physical_temperature_or_seconds_inferred': False}


def exact_survival():
    """Exact reduced CP map P Phi(P rho P) P; no survivor normalization."""
    indices = corner_indices(); weights = native_clock.parameters().port_weights
    rows = {}
    for law in ('C', 'E'):
        coefficients = []
        for i in indices:
            for j in indices:
                out = [[F() for _ in indices] for _ in indices]
                for edge in native.SEAMS:
                    for coin, _, pairs in native.raw_maps(law, edge):
                        p = F(weights[edge[0]] if coin else weights[edge[1]],
                              30*(weights[edge[0]]+weights[edge[1]]))
                        mapping = dict(pairs)
                        if i in mapping and j in mapping and mapping[i] in indices and mapping[j] in indices:
                            out[indices.index(mapping[i])][indices.index(mapping[j])] += p
                coefficients.append(out)
        rows[law] = coefficients
    expected_diagonal = ((F(859, 990), F(4, 495)), (F(5, 198), F(404, 495)))
    for law, alpha in (('C', F(1519, 1980)), ('E', F(37, 60))):
        value = rows[law]
        if (tuple(value[0][k][k] for k in range(2)) != expected_diagonal[0]
                or tuple(value[3][k][k] for k in range(2)) != expected_diagonal[1]
                or value[1] != [[0, alpha], [0, 0]]
                or value[2] != [[0, 0], [alpha, 0]]):
            raise ArithmeticError('the complete native continuation no longer has the claimed corner map')
    return {'population_kernel_rows': [['859/990', '5/198'], ['4/495', '404/495']],
            'coherence_multiplier': {'C': '1519/1980', 'E': '37/60'},
            'initial_corner_mass': '264/7505', 'corner_reference': ['25/33', '8/33'],
            'post_probe_return_probability': '85/99', 'post_probe_outside_probability': '14/99',
            'energy_r_record_probability_given_initial_corner_formula':
                '85/198 + (17/66)*alpha_law*cos(s*log(25/8))',
            'complete_probability_multiplies_initial_corner_mass': True,
            'reported_channel_conditioned_on_return': False}


def _pointer_plan():
    """Map the two accepted addresses to one free bit and eight fixed bits.

    The final pointer is qubit9. The checked 9-qubit controlled-X macro acts
    on eight work controls and that pointer, leaving the free work bit alone.
    This measures a rank-two projector without resolving its energy basis.
    """
    a, b = _corner_codes(); difference = a ^ b
    pivot = (difference & -difference).bit_length()-1
    fan = [pulse for bit in range(9) if bit != pivot and difference >> bit & 1
           for pulse in compiler.cnot(pivot, bit)]
    transformed = a ^ ((difference ^ (1 << pivot)) if a >> pivot & 1 else 0)
    controls = tuple(bit for bit in range(9) if bit != pivot)
    before = fan+[compiler.pulse('X', math.pi/2, [bit]) for bit in controls if transformed >> bit & 1]
    wires = (*controls, 9)
    macro = native_routing._canonical_word()
    after = compiler.inverse(before)
    literal = before+compiler.remap(macro, wires)+after
    return before, macro, wires, pivot, after, literal


@lru_cache(maxsize=4)
def _pointer_payload(codes, word_digest):
    if codes != _corner_codes():
        raise ArithmeticError('source corner chart changed')
    before, macro, wires, free, after, literal = _pointer_plan()
    if word_digest != hashlib.sha256(json.dumps(literal, sort_keys=True).encode()).hexdigest():
        raise ArithmeticError('pointer literal word changed during snapshot')
    raw_macro = compiler.execute(macro, 9)
    target = np.eye(512, dtype=complex)
    target[[0, 256]] = target[[256, 0]]
    if np.linalg.norm(raw_macro-target, 'fro') > float(POINTER_GATE)/4:
        raise ArithmeticError('complete source pointer macro exceeds its gate')
    # All512 original inputs with a blank pointer; no padding column omitted.
    initial = np.zeros((1024, 512), complex); initial[:512] = np.eye(512)
    out = compiler.execute(before, 10, initial)
    # Tensor chart of the actual macro word; not a freely supplied operation.
    out = compiler._apply(out, raw_macro, list(reversed(wires)), 10)
    out = compiler.execute(after, 10, out)
    mask = np.zeros(512); mask[list(codes)] = 1
    ideal = np.vstack((np.diag(1-mask), np.diag(mask)))
    # Exact rational squared Frobenius difference from original float entries.
    squared = F()
    for value in (out-ideal).flat:
        squared += F(float(value.real))**2+F(float(value.imag))**2
    if squared > POINTER_GATE**2:
        raise ArithmeticError('complete source corner instrument exceeds its gate')
    resources = compiler.resources(literal)
    info = {'source_pulses_per_herald': len(literal), 'full_input_columns': 512,
            'pointer_output_dimension': 1024, 'operator_error_upper': str(POINTER_GATE),
            'full_isometry_squared_frobenius_upper': str(POINTER_GATE**2),
            'padding_inputs_checked': 446, 'outcome_probabilities_renormalized': False,
            'macro_word_sha256': hashlib.sha256(json.dumps(macro, sort_keys=True).encode()).hexdigest(),
            'literal_word_sha256': hashlib.sha256(json.dumps(literal, sort_keys=True).encode()).hexdigest(),
            'resources': {key: value for key, value in resources.items() if type(value) in (int, bool)}}
    return out.dtype.str, out.tobytes(), json.dumps(info, sort_keys=True)


def pointer():
    literal = _pointer_plan()[-1]
    digest = hashlib.sha256(json.dumps(literal, sort_keys=True).encode()).hexdigest()
    dtype, payload, info = _pointer_payload(_corner_codes(), digest)
    out = np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(1024, 512)
    return (out[:512].copy(), out[512:].copy()), json.loads(info)


def _ideal_pointer():
    mask = np.zeros(512); mask[list(_corner_codes())] = 1
    return np.diag(1-mask), np.diag(mask)


@lru_cache(maxsize=4)
def _zero(size):
    return np.frombuffer(np.zeros((size, size), complex).tobytes(), dtype=complex).reshape(size, size)


def _put(outputs, key, matrix):
    # Sparse representation of exact-zero outcome blocks. The declared raw
    # alphabet and accessor below distinguish these from missing positive maps.
    outputs[key] = matrix if np.any(matrix) else _zero(len(matrix))


def execute(density, relative_weights, law, *, raw=True, clock_on=True):
    """Complete two-herald history on an arbitrary source/reference input.

    Input is H66 tensor R, with R dimension1 or2 for numerical execution.
    Operator certificates cover arbitrary finite R. Rejections halt with the
    original quantum output. All source padding escapes retain the R marginal.
    """
    density = native_clock._validate_input(density)
    law = native.records.candidate(law)
    if type(raw) is not bool or type(clock_on) is not bool:
        raise ValueError('original Boolean execution choices required')
    weights = learned_flow._weights(relative_weights)
    arrays = learned_flow.source_arrays(weights) if raw else learned_flow._ideal_arrays(weights)
    herald = pointer()[0] if raw else _ideal_pointer()
    outputs = {}; r = len(density)//66
    outside, abort = native_clock._export(herald[0], density)
    _put(outputs, ('initial_outside',), outside); _put(outputs, ('abort', 'initial_outside'), abort)
    inside, abort = native_clock._export(herald[1], density)
    _put(outputs, ('abort', 'initial_inside'), abort)
    inside, abort = native_clock._export(arrays['probe'], inside)
    _put(outputs, ('abort', 'probe'), abort)
    for edge in native.SEAMS:
        for label, pairs, root, probability in native_clock.native_branches(law, edge):
            prefix = (edge, label[0], label[1])
            if not pairs:
                continue
            state = native_clock._native(inside, pairs, (root*root if raw else probability)/30)
            outside, abort = native_clock._export(herald[0], state)
            _put(outputs, ('outside', *prefix), outside)
            _put(outputs, ('abort', 'outside', *prefix), abort)
            state, abort = native_clock._export(herald[1], state)
            _put(outputs, ('abort', 'inside', *prefix), abort)
            if clock_on:
                state, abort = native_clock._export(arrays['clock'], state)
            else:
                abort = np.zeros((r, r), complex)
            _put(outputs, ('abort', 'clock', *prefix), abort)
            state, abort = native_clock._export(arrays['analyzer'], state)
            _put(outputs, ('abort', 'analyzer', *prefix), abort)
            for final, final_pairs, root2, probability2 in native_clock.native_branches(law, PORTS):
                if not final_pairs:
                    continue
                output = native_clock._native(state, final_pairs, root2*root2 if raw else probability2)
                _put(outputs, ('read', *prefix, final[0], final[1]), output)
    return outputs


def record_probability(outputs, law, energy=0):
    """Decode pre-detector energy from coin+current; no basis postselection."""
    if type(energy) is not int or energy not in (0, 1):
        raise ValueError('energy index0(r) or1(d) required')
    total = 0.
    for key, state in outputs.items():
        if key[0] != 'read':
            continue
        coin, raw = key[-2:]
        semantic = circuit.record_values(law, coin).get(raw)
        if semantic is None:
            continue
        current = native.records.record_current(law, (*PORTS, coin, semantic))
        # Occupation of port1 before repair equals coin+current on q=1.
        if coin+current == 1-energy:
            total += float(np.trace(state).real)
    return total


def raw_output(outputs, law, edge, coin, raw, final_coin, final_raw):
    """Every complete two-record raw word, including impossible/zero slots."""
    law = native.records.candidate(law); edge = native.records.seam(edge)
    if edge not in native.SEAMS:
        raise ValueError('native continuation seam required')
    width = 1 << circuit.RECORD_BITS[law]
    for c, r in ((coin, raw), (final_coin, final_raw)):
        if type(c) is not int or c not in (0,1) or type(r) is not int or not 0 <= r < width:
            raise ValueError('original complete native coin/raw slot required')
    first = next(row for row in native_clock.native_branches(law, edge) if row[0][:2] == (coin, raw))
    last = next(row for row in native_clock.native_branches(law, PORTS) if row[0][:2] == (final_coin, final_raw))
    key = ('read', edge, coin, raw, final_coin, final_raw)
    if not first[1] or not last[1]:
        return _zero(len(outputs['initial_outside',]))
    if key not in outputs:
        raise ValueError('missing live complete particle-clock raw block')
    return outputs[key]


def source_budget(weights):
    _, pointer_info = pointer()
    inherited = learned_flow.source_budget(weights)
    # Two heralds, two rotations, two native frames, and one phase. Source
    # failures have their own disjoint flags; physical hardware errors absent.
    error, gain = F(), F(1)
    epsilons = (POINTER_GATE, native_clock.ROTATION_GATE, native_clock.NATIVE_GATE,
                POINTER_GATE, native_clock.WORD_GATE, native_clock.ROTATION_GATE,
                native_clock.NATIVE_GATE)
    for eps in epsilons:
        error = (1+eps)**2*error+eps*(2+eps)/2
        gain *= (1+eps)**2
    program = learned_flow.source_program(weights)
    pulses = 2*pointer_info['source_pulses_per_herald']+sum(len(program[x]) for x in ('probe', 'clock', 'analyzer'))
    pulse_gate, failure_mass = native_clock._processor_failure_bound(pulses)
    all_words = (*_pointer_plan()[-1], *(p for name in ('probe', 'clock', 'analyzer') for p in program[name]))
    for kind, angle, bits in {(p['kind'], p['angle'], len(p['targets'])) for p in all_words}:
        leakage, excess = learned_flow._processor_bounds(compiler.native_pulse(kind, angle), 1 << bits)
        if max(leakage, excess) > pulse_gate:
            raise ArithmeticError('particle-clock processor failure/gain gate failed')
    native_costs = acquisition.native_word_costs()
    native_pulses = max(v[law] for v in native_costs.values() for law in ('C', 'E'))
    detector = acquisition.native_word_costs((PORTS,))
    native_pulses += max(detector[PORTS].values())
    return {'coherent_source_pulses_upper': pulses, 'coherent_load_export_calls_upper': 2*pulses,
            'native_continuation_attempts': 1, 'native_detector_attempts': 1,
            'native_routed_word_pulses_upper_count_only': native_pulses,
            'total_supplied_word_pulses_upper': pulses+native_pulses,
            'native_raw_record_bits_upper': 20, 'corner_bits': 2,
            'coherent_register_qubits_including_pointer': 10,
            'peak_native_source_and_raw_qubits': 12,
            'continuation_seam_selector_calls_upper': 1,
            'blank_pointer_initializations_and_reads_upper': 2,
            'clock_status_and_first_processor_failure_metadata': 'additional bounded counters and flags',
            'source_instrument_halfdiamond_upper': str(error+failure_mass/2),
            'source_channel_norm_upper': str(gain+failure_mass),
            'proper_code_first_failure_mass_upper': str(failure_mass),
            'proper_code_failure_mass_included_but_outputs_not_materialized': True,
            'source_address_export_aborts_materialized': True,
            'physical_source_error_included': False,
            'field_encoding_pulse_realization_claimed': False,
            'pointer': pointer_info}


def admission(counts):
    """Record-only, fixed-cap decision; all errors/stops are unconditional."""
    estimate = infer_clock(counts)
    if estimate['status'] != 'accept':
        return estimate
    preparation = record_clock.preparation_certificate(counts)
    if preparation['status'] != 'accept' or F(preparation['half_diamond_upper']) > record_clock.PREPARATION_CAP:
        return {'status': 'stop', 'reason': 'record_preparation_cap_failed'}
    try:
        budget = source_budget(estimate['relative_weights'])
    except ArithmeticError:
        return {'status': 'stop', 'reason': 'source_certification_failed'}
    if (F(budget['source_instrument_halfdiamond_upper']) > F(1, 1000000)
            or F(budget['source_channel_norm_upper']) > F(1000001, 1000000)):
        return {'status': 'stop', 'reason': 'fixed_source_cap_failed'}
    return {'status': 'accept', 'inference': estimate, 'preparation': preparation, 'source_budget': budget}


def continue_from_records(counts, terminal_index, law):
    """Simulator continuation; inference/compilation never reads terminal_index."""
    law = native.records.candidate(law)
    if type(terminal_index) is not int or not 0 <= terminal_index < 66:
        raise ValueError('original terminal binary-source index required')
    decision = admission(counts)
    if decision['status'] != 'accept':
        density = np.zeros((66,66), complex); density[terminal_index, terminal_index] = 1
        return decision, {('stop', decision['reason']): density}
    density = record_clock.continued_source_density(terminal_index)
    return decision, execute(density, decision['inference']['relative_weights'], law)


def bindings():
    paths = set(learned_flow.source_sha256()) | set(inference.source_bindings())
    paths.update(('model/particle_clock/clock.py', 'model/particle_clock/test_clock.py',
                  'model/particle_clock/README.md', 'model/particle_clock/__init__.py',
                  'model/particle_clock/__main__.py', 'model/local_repair_family/record_clock.py',
                  'model/local_repair_family/acquisition.py', 'model/common_source/native_routing.py',
                  'model/bound_pair/fiber.py', 'model/interacting_field/lattice.py'))
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def report():
    field = clock_subspace()
    if (np.linalg.norm(field['restriction']-np.diag([0.,1.])) > 1e-12
            or any(v > 1e-12 or not math.isfinite(v) for v in field['residuals'].values())):
        raise ArithmeticError('the same interacting field no longer supplies the declared clock subspace')
    source, counts = acquisition.acquire(native_clock.parameters())
    estimate = infer_clock(counts)
    if estimate['status'] != 'accept':
        raise ArithmeticError('retained native training did not admit the particle clock')
    weights = estimate.pop('relative_weights')
    estimate.pop('beta_times_U_estimate_diagnostic')
    estimate['beta_times_U_estimate'] = 'log('+estimate['energy_ratio_estimate']+')'
    preparation = record_clock.preparation_certificate(counts)
    if preparation['status'] != 'accept' or F(preparation['half_diamond_upper']) > record_clock.PREPARATION_CAP:
        raise ArithmeticError('record-only preparation gate failed')
    budget = source_budget(weights)
    actual_input = record_clock.continued_source_density(source['private_replay']['terminal_state_index'])
    target = np.diag([float(p) for p in native.reference(native_clock.parameters())]).astype(complex)
    diagnostics = {}
    q = F(264, 7505)
    for law in ('C', 'E'):
        actual = execute(actual_input, weights, law)
        ideal = execute(target, tuple(F(w) for w in native_clock.parameters().port_weights), law, raw=False)
        if actual.keys() != ideal.keys():
            raise ArithmeticError('complete particle-clock alphabets differ')
        distance = sum(np.linalg.norm(actual[k]-ideal[k], 'nuc') for k in actual)/2
        mass = sum(float(np.trace(x).real) for x in actual.values())
        alpha = F(exact_survival()['coherence_multiplier'][law])
        probability = float(q)*(85/198+float(F(17,66)*alpha)*math.cos(math.log(25/8)))
        if abs(record_probability(ideal, law)-probability) > 1e-12:
            raise ArithmeticError('complete quantum instrument disagrees with the exact particle energy response')
        on_event = (F(budget['source_channel_norm_upper'])*record_clock.PREPARATION_CAP
                    + F(budget['source_instrument_halfdiamond_upper'])+PHASE_CAP)
        if not math.isfinite(distance) or distance > float(on_event) or abs(mass-1) > 1e-6:
            raise ArithmeticError('complete particle-clock instrument exceeds its fixed error gate')
        energy_r = record_probability(actual, law)
        energy_d = record_probability(actual, law, 1)
        initial_outside = float(np.trace(actual['initial_outside',]).real)
        if not (.017 < energy_r < .019 and .011 < energy_d < .014 and .96 < initial_outside < .97):
            raise ArithmeticError('complete particle-clock probabilities left their fixed diagnostic intervals')
        diagnostics[law] = {'complete_stored_blocks': len(actual),
            'source_vs_true_reference_clock_half_trace_upper': str(on_event),
            'complete_output_mass_error_upper': str(2*F(budget['source_instrument_halfdiamond_upper'])),
            'energy_r_record_probability_interval': ['17/1000', '19/1000'],
            'energy_d_record_probability_interval': ['11/1000', '14/1000'],
            'initial_outside_probability_interval': ['24/25', '97/100'],
            'exact_ideal_response_and_complete_quantum_outputs_checked': True}
    alpha = inference.fixed_budget()['statistical_failure_upper']
    compiler_cap = F(1, 1000000); gain_cap = F(1000001, 1000000)
    if F(budget['source_instrument_halfdiamond_upper']) > compiler_cap or F(budget['source_channel_norm_upper']) > gain_cap:
        raise ArithmeticError('fixed all-transcript source caps exceeded')
    total = gain_cap*record_clock.PREPARATION_CAP+compiler_cap+PHASE_CAP+alpha*(gain_cap+1)/2
    return {'schema': 'oph.particle_clock.v1', 'source_sha256': bindings(),
            'field_dictionary': {'source_states': [[1,2], [7,2]], 'field_states': ['r', 'd'],
                'energy_matrix': [['0','0'],['0','1']], 'total_momentum': ['pi','pi','pi'],
                'domain': 'even periodic L>=4; two fermions; singlet relative amplitudes sqrt(2)sin(k_x) and1',
                'kinetic_identity': 'h(pi-k)=h(k); (sigma_a tensor I + I tensor sigma_a) singlet=0',
                'contact_identity': 'D d=d; D r=0 by zero mean',
                'same_field_source_module': 'model/bound_pair/fiber.py',
                'executed_fiber_side': field['side'],
                'executed_fiber_dimension': len(field['basis']),
                'field_operator_map_residual_gate': '1/1000000000000',
                'beta_times_U': 'log(25/8)', 'faithful_algebra': 'B(span{r,d}) = M2',
                'Heisenberg_modular_parameter_relation': 't=-beta*s',
                'supplied_source_to_field_isometry': True,
                'faithful_on_full_field_or_local_regional_algebras': False},
            'record_inference': estimate, 'preparation': preparation,
            'exact_continuation': exact_survival(), 'source_budget': budget,
            'complete_readouts': diagnostics,
            'complete_average_half_trace_upper': str(total),
            'error_bound_conditioned_on_herald_or_training_acceptance': False,
            'training_stream_sha256': source['stream_sha256'],
            'training_attempts': source['protocol']['total_events'],
            'scope': {'one_same_excitation_energy_and_clock_operator_map': True,
                'native_repair_equals_interacting_Hamiltonian_flow': False,
                'full_interacting_Gibbs_reference_prepared': False,
                'physical_temperature_time_energy_normalization': False,
                'mass_in_electronvolts_or_physical_particle_identification': False,
                'interacting_continuum_limit_or_integrated_gate_accepted': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args(); value = report()
    if args.check is not None:
        if strict_json.load(args.check) != value:
            raise ValueError('particle-clock receipt differs from fresh complete source execution')
    elif args.write:
        OUTPUT.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n', encoding='utf8')
    else:
        print(json.dumps(value, sort_keys=True, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
