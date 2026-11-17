"""Spatially resolved neutral coherence in native attempt and movement counts.

The moving coordinate is a common occupation, while the two signed-difference
endpoints stay pinned. Phase preparation, a complete joint phase/occupation
readout and all-seam current access are declared services. No physical particle
mass or Lorentzian continuation is inferred from these dissipative responses.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
from itertools import combinations
from pathlib import Path

import numpy as np
import sympy as sp

from ..federated_quantum_sector import construction as quantum
from ..native_relational_readout import construction as stopped
from ..native_quantum_history import construction as history
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..primitive_observer import runtime
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*stopped.SOURCE_PATHS,
    'model/federated_quantum_sector/construction.py',
    'model/native_action/construction.py',
    'model/native_effective_fields/neutral_transport.py'))))


def endpoints(pair):
    if (type(pair) is not tuple or len(pair) != 2
            or any(type(x) is not int for x in pair) or not 0 <= pair[0] < pair[1] < 12):
        raise ValueError('two distinct ordered source ports required')
    return pair


def common_ports(pair):
    pair = endpoints(pair)
    return tuple(j for j in range(12) if j not in pair)


def addresses(pair, common):
    pair = endpoints(pair)
    if type(common) is not int or common not in common_ports(pair):
        raise ValueError('a common occupied port outside the two endpoints required')
    return tuple(history.INDEX[tuple(int(i in (p, common)) for i in range(12))] for p in pair)


def adjacency(pair):
    ports, edges = common_ports(pair), set(history.mixing.seams())
    return sp.ImmutableMatrix([[int(tuple(sorted((i, j))) in edges) for j in ports] for i in ports])


def killing(pair, law):
    endpoints(pair); history.candidate(law)
    return 5 if law == 'C' else 10-int(pair in history.mixing.seams())


def response(pair, law, *, clock='attempt'):
    """Exact response of |a,j><b,j|; columns label its common initial port.

    Movement clock counts every carrier movement, not one port's incident
    events. Matrix entries are phase contrasts, not exclusive event masses.
    """
    a, v = adjacency(pair), killing(pair, law)
    degrees = [sum(a[:, j]) for j in range(10)]
    if clock == 'attempt':
        return sp.ImmutableMatrix(sp.eye(10)-sp.diag(*degrees)/60+a/60-sp.eye(10)*sp.Rational(v, 30))
    if clock == 'movement':
        return sp.ImmutableMatrix(a*sp.diag(*(sp.Rational(1, 2*v+d) for d in degrees)))
    raise ValueError('clock must be attempt or movement')


def native_response(pair, law, current=None):
    """Independently enumerate actual endpoint instruments, preserving labels."""
    ports = common_ports(pair)
    units = {((1 << pair[0]) | (1 << j), (1 << pair[1]) | (1 << j)): n
             for n, j in enumerate(ports)}
    out = sp.zeros(10)
    for (left, right), col in units.items():
        for seam in history.mixing.seams():
            for key, value in quantum.native_unit(law, left, right, seam, current).items():
                if key not in units:
                    raise ArithmeticError('native signed-difference endpoints moved')
                out[units[key], col] += sp.Rational(value.numerator, 30*value.denominator)
    return sp.ImmutableMatrix(out)


def preparation(pair, common, sign=1):
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('an exact plus or minus preparation sign required')
    i, j = addresses(pair, common)
    return {(i, i): F(1, 2), (j, j): F(1, 2),
            (i, j): F(sign, 2), (j, i): F(sign, 2)}


def readout(pair):
    """Twenty orthogonal +/- projectors and their entire 46-D complement."""
    projectors, used = [], np.zeros((66, 66), complex)
    for common in common_ports(pair):
        i, j = addresses(pair, common)
        for sign in (1, -1):
            p = np.zeros((66, 66), complex)
            p[i, i] = p[j, j] = .5
            p[i, j] = p[j, i] = sign*.5
            projectors.append(p); used += p
    projectors.append(np.eye(66)-used)
    return runtime.Instrument(kraus=tuple((p,) for p in projectors), visible=True,
        resources=(('spatial_phase_readouts', 1),),
        binding={'endpoints': list(pair), 'common_ports': list(common_ports(pair)),
                 'outcomes': 'common-port plus/minus pairs, then full complement',
                 'declared_coherent_analyzer': True, 'native_pulse_words_emitted': False})


def exact_probabilities(state, pair):
    values, used = [], set()
    for common in common_ports(pair):
        i, j = addresses(pair, common); used.update((i, j))
        diagonal = state.get((i, i), F(0))+state.get((j, j), F(0))
        phase = state.get((i, j), F(0))+state.get((j, i), F(0))
        values.extend(((diagonal+phase)/2, (diagonal-phase)/2))
    values.append(sum((v for (i, j), v in state.items() if i == j and i not in used), F(0)))
    if min(values) < 0:
        raise ArithmeticError('phase contrasts were substituted for probabilities')
    return tuple(values)


def exact_first_current(marked, initial, *, first=None, tail=None):
    """Actual CP Kraus series; no bounds means its complete infinite sum."""
    if first is not None and tail is not None:
        raise ValueError('choose a first-current count or an unresolved tail')
    if first is not None:
        stopped.horizon(first)
    if tail is not None:
        stopped.horizon(tail)
        return {(i, j): value*F(marked.no_current_counts[i][j], 60)**tail
                for (i, j), value in initial.items()}
    out = defaultdict(F)
    for row in marked.moving_maps:
        mapping = dict(row)
        for (i, j), value in initial.items():
            if i in mapping and j in mapping:
                w = F(marked.no_current_counts[i][j], 60)
                factor = w**(first-1)/60 if first is not None else F(1, 60)/(1-w)
                out[mapping[i], mapping[j]] += value*factor
    return dict(out)


def program(marked, pair, common, sign=1, cap=3):
    stopped.horizon(cap, maximum=3)
    latch = runtime.Instrument(kraus=((np.zeros((66, 66), complex),), (np.eye(66),)),
                              visible=True, resources=(('latched_flag_copies', 1),))
    events = [runtime.Event('native_current_1', instruments=(marked.instrument,))]
    events.extend(runtime.Event(f'native_current_or_hold_{n+1}',
                  instruments=(marked.instrument, latch), control=n-1) for n in range(1, cap))
    events.append(runtime.Event('complete_spatial_phase_readout', instruments=(readout(pair),)))
    return runtime.Program(f'neutral_transport_{marked.law}_{pair}_{common}_{sign}_{cap}',
        66, stopped.dense(preparation(pair, common, sign)), tuple(events),
        binding={'clock': 'first global moving current with all-seam access',
                 'initial_common_port': common, 'endpoints': list(pair), 'phase_sign': sign,
                 'attempt_cap': cap, 'coherent_preparation_and_analyzer_supplied': True})


def execution(marked, pair=(0, 1), common=2, sign=1, cap=3):
    p = program(marked, pair, common, sign, cap)
    result = runtime.execute(p)
    initial = preparation(pair, common, sign)
    exact = {n: exact_first_current(marked, initial, first=n) for n in range(1, cap+1)}
    exact[None] = exact_first_current(marked, initial, tail=cap)
    rows, error, mass = [], 0., F(0)
    complete = [F(0)]*21; unresolved = [F(0)]*21
    for branch in result.branches:
        flags, outcome = branch.record[:-1], branch.record[-1]
        possible = all(a <= b for a, b in zip(flags, flags[1:]))
        first = flags.index(1)+1 if 1 in flags else None
        state = exact[first] if possible else {}
        probability = exact_probabilities(state, pair)[outcome]
        k = p.events[-1].instruments[0].kraus[outcome][0]
        wanted = k@stopped.dense(state)@k.conj().T
        error = max(error, float(np.linalg.norm(branch.density-wanted, 'fro')),
                    abs(branch.probability-float(probability)))
        rows.append([list(branch.record), str(probability)])
        mass += probability
        if possible:
            ticks = cap if first is None else first
            if dict(branch.resources).get('native_attempts') != ticks:
                raise ArithmeticError('native cost includes a halted repair')
            (unresolved if first is None else complete)[outcome] += probability
    if mass != 1 or len(rows) != 21*2**cap:
        raise ArithmeticError('the full spatial readout lost outcomes or its timeout')
    prefix = runtime.execute(p, stop_after=1)
    replay = runtime.execute(p, checkpoint=runtime.checkpoint_from_json(p, runtime.checkpoint_to_json(prefix)))
    return {'candidate': marked.law, 'phase_sign': sign, 'initial_common_port': common,
        'endpoints': list(pair), 'attempt_cap': cap, 'complete_records': len(rows),
        'terminal_dimension_per_record': 66,
        'exact_complete_record_probability_sha256': hashlib.sha256(canonical(rows).encode()).hexdigest(),
        'first_movement_joint_outcomes': list(map(str, complete)),
        'timeout_joint_outcomes': list(map(str, unresolved)),
        'timeout_mass': str(sum(unresolved)), 'uniform_timeout_bound': str(F(13, 15)**cap),
        'terminal_state_and_probability_error_diagnostic': stopped.bounded(error, 'neutral transport'),
        'checkpoint_complete_distance_diagnostic': stopped.bounded(runtime.complete_distance(result, replay), 'neutral checkpoint')}


def exact_audit():
    checks, spectra = 0, []
    for pair in combinations(range(12), 2):
        a = adjacency(pair)
        for law in ('C', 'E'):
            actual = native_response(pair, law)
            zero = native_response(pair, law, 0)
            if actual != response(pair, law) or zero != sp.diag(*zero.diagonal()):
                raise ArithmeticError('source transfer disagrees with the punctured-graph action')
            moving = actual-zero
            if moving != a/60 or moving*(sp.eye(10)-zero).inv() != response(pair, law, clock='movement'):
                raise ArithmeticError('stopped spatial response lost its native CP derivation')
            checks += 10
        if pair in ((0, 1), (0, 2), (0, 3)):
            laplacian = sp.diag(*(sum(a[:, j]) for j in range(10)))-a
            spectra.append({'endpoints': list(pair),
                'punctured_laplacian_characteristic_polynomial': str(sp.factor(laplacian.charpoly().as_expr())),
                'attempt_eigenvalues': {law: f'1-{killing(pair, law)}/30-lambda/60' for law in ('C', 'E')}})
    return {'all_endpoint_pairs': 66, 'source_columns_checked_across_both_laws': checks,
        'all_marked_zero_current_columns_are_diagonal': True,
        'every_coherent_movement_uses_one_punctured_source_edge': True,
        'shell_spectra': spectra}


def first_tick_predictions(marked, pair=(0, 1)):
    matrix = response(pair, marked.law, clock='movement')
    rows = []
    for col, common in enumerate(common_ports(pair)):
        signs = []
        for sign in (1, -1):
            state = exact_first_current(marked, preparation(pair, common, sign))
            probs = exact_probabilities(state, pair)
            if sum(probs) != 1:
                raise ArithmeticError('the ideal first-tick spatial readout is incomplete')
            contrast = tuple(probs[2*k]-probs[2*k+1] for k in range(10))
            if contrast != tuple(sign*F(matrix[k, col]) for k in range(10)):
                raise ArithmeticError('the CP readout does not recover the spatial contrast')
            signs.append(list(map(str, probs)))
        rows.append({'initial_common_port': common, 'plus_preparation_complete_probabilities': signs[0],
                     'minus_preparation_complete_probabilities': signs[1]})
    return {'candidate': marked.law, 'rows': rows,
            'first_tick_contrast_matrix': [[str(v) for v in matrix.row(i)] for i in range(10)],
            'normalization_or_postselection_of_contrasts': False}


def discrimination(sources):
    """The complete joint analyzer attains the actual state discrimination."""
    pair, common = (0, 1), 2
    states = {law: exact_first_current(source, preparation(pair, common))
              for law, source in sources.items()}
    difference = {key: states['C'].get(key, F(0))-states['E'].get(key, F(0))
                  for key in states['C'].keys() | states['E'].keys()}
    pairs = [addresses(pair, k) for k in common_ports(pair)]
    allowed = set(pairs) | {(j, i) for i, j in pairs}
    if any(value and key not in allowed for key, value in difference.items()):
        raise ArithmeticError('C/E state difference leaves the orthogonal phase blocks')
    if any(difference.get((i, j), F(0)) != difference.get((j, i), F(0)) for i, j in pairs):
        raise ArithmeticError('the exact C/E phase contrast is not Hermitian')
    # Each disjoint 2x2 block is [[0,c],[c,0]], with half trace norm |c|.
    distance = sum((abs(difference.get(key, F(0))) for key in pairs), F(0))
    p, q = (exact_probabilities(states[law], pair) for law in ('C', 'E'))
    variation = sum((abs(x-y) for x, y in zip(p, q, strict=True)), F(0))/2
    if distance != variation or distance != F(4, 69):
        raise ArithmeticError('the complete analyzer fails to attain exact discrimination')
    return {'ideal_first_tick_output_trace_distance': str(distance),
            'complete_measurement_total_variation': str(variation),
            'initial_common_port': common, 'coherent_endpoints': list(pair),
            'same_diagonal_populations_and_current_waiting_law': True}


def run():
    sources = {law: stopped.source(law) for law in ('C', 'E')}
    return {'schema': 'oph.native_neutral_transport.v1', 'all_pass': True,
        'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'exact_source_attachment': exact_audit(),
        'first_tick_predictions': [first_tick_predictions(s) for s in sources.values()],
        'candidate_discrimination': discrimination(sources),
        'finite_executions': [execution(s, sign=sign) for s in sources.values() for sign in (1, -1)],
        'scope': {'source_selected_nonflat_neutral_coherence_transfer': True,
            'moving_common_occupation_not_motion_of_signed_difference_endpoints': True,
            'same_native_C_or_E_instrument_during_evolution': True,
            'uniform_scheduler_and_all_seam_current_access_supplied': True,
            'coherent_preparation_joint_analyzer_and_halt_controller_supplied': True,
            'time_parameter': 'integer native attempts or global recorded movements',
            'physical_spatial_coordinate_clock_or_particle_mass_identified': False,
            'charged_particle_propagation_or_Lorentzian_continuum_derived': False,
            'all_finite_theorem_premises_discharged': False,
            'native_pulse_words_emitted': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path); mode.add_argument('--check', type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native neutral transport receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
