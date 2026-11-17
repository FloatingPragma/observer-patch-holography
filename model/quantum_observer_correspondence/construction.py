"""Exact coherent stopped instruments and their transported matrix cones.

The same native C carrier supplies every scheduled and feedback repair.
The declared union apparatus is an intervention, not a replacement dynamics.
All coefficients in the retained experiment are rational.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
from pathlib import Path

import sympy as sp

from ..observer import native, quantum_mixing as mixing
from ..native_quantum_history import construction as quantum
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
STATES = quantum.STATES
INDEX = {state: i for i, state in enumerate(STATES)}
REGIONS = ((0, 1), (0, 5))
UNION = (0, 1, 5)
P = F(2, 5)
SOURCE_PATHS = tuple(sorted(set((*quantum.SOURCE_PATHS,
    'model/quantum_observer_correspondence/__init__.py',
    'model/quantum_observer_correspondence/__main__.py',
    'model/quantum_observer_correspondence/construction.py',
    'model/quantum_observer_correspondence/clock_response.py'))))


def clean(values):
    return {key: value for key, value in values.items() if value}


def add(left, right, scale=F(1)):
    result = defaultdict(F, left)
    for key, value in right.items():
        result[key] += scale*value
    return clean(result)


def trace(state):
    return sum((value for (i, j, r, s), value in state.items()
                if i == j and r == s), F(0))


def mark(seam, current):
    return (tuple(int(bool(set(seam).intersection(region))) for region in REGIONS),
            tuple(current*(int(v == seam[1])-int(v == seam[0])) for v in UNION))


def address(ports):
    return INDEX[tuple(int(v in ports) for v in range(12))]


def preparation(reference=False, dephased=False):
    if type(reference) is not bool or type(dephased) is not bool:
        raise ValueError('explicit Boolean preparation choices required')
    a, b = address((0, 1)), address((0, 5))
    ket = ((a, 0), (b, int(reference)))
    return {(i, j, r, s): F(1, 2) for (i, r), (j, s) in product(ket, repeat=2)
            if not dephased or i == j}


def entangled_preparation():
    """A source/reference Bell pair whose coherence survives the final read."""
    a, b = address((0, 2)), address((0, 3))
    return {(a, a, 0, 0): F(1, 2), (b, b, 1, 1): F(1, 2),
            (a, b, 0, 1): F(1, 2), (b, a, 1, 0): F(1, 2)}


@lru_cache(maxsize=1)
def source():
    """Actual C Kraus numerator columns, then all matrix-unit transitions.

    Coin and current index a Kraus operator on each seam. Summing operators
    with different records before forming their products would change C.
    """
    seams = mixing.seams()
    branches = {}
    for seam in seams:
        groups = {}
        for coin in (0, 1):
            for i, state in enumerate(STATES):
                output, record = native.native(state, seam, coin)
                if record != (coin, state[seam[0]]-output[seam[0]]):
                    raise ArithmeticError('native current semantics changed')
                groups.setdefault(record, {})[i] = INDEX[output]
        branches[seam] = tuple(sorted(groups.items()))
        if any(len(set(columns.values())) != len(columns) for columns in groups.values()):
            raise ArithmeticError('a current Kraus numerator is not a partial permutation')
    remote, incidents, feedback = defaultdict(Counter), {}, {}
    for seam, groups in branches.items():
        incident = bool(set(seam).intersection(UNION))
        for (_, current), columns in groups:
            label = mark(seam, current)
            for (i, a), (j, b) in product(columns.items(), repeat=2):
                before, after = (i, j), (a, b)
                if tuple(x-y for x, y in zip(STATES[i], STATES[j], strict=True)) != tuple(
                        x-y for x, y in zip(STATES[a], STATES[b], strict=True)):
                    raise ArithmeticError('native current branch changes its difference sector')
                if incident:
                    incidents.setdefault(label, defaultdict(Counter))[before][after] += 1
                else:
                    remote[before][after] += 1
                if seam in ((0, 1), (0, 5)):
                    feedback.setdefault(seam, {}).setdefault(label, defaultdict(Counter))[before][after] += 1
    if len(STATES) != 66 or len(seams) != 30:
        raise ArithmeticError('the source sector or schedule changed')
    for i, j in product(range(66), repeat=2):
        remote_trace = sum(count for (a, b), count in remote.get((i, j), {}).items() if a == b)
        incident_trace = sum(count for matrix in incidents.values()
            for (a, b), count in matrix.get((i, j), {}).items() if a == b)
        if remote_trace != 36*int(i == j):
            raise ArithmeticError('state-independent remote survival changed')
        if incident_trace != 24*int(i == j):
            raise ArithmeticError('incident effect differs from (2/5) identity')
    return {'branches': branches, 'remote': remote, 'incidents': incidents,
            'feedback': feedback, 'marks': tuple(sorted(incidents))}


def apply_count(matrix, state, denominator):
    result = defaultdict(F)
    for (i, j, r, s), value in state.items():
        for (a, b), count in matrix.get((i, j), {}).items():
            result[a, b, r, s] += value*F(count, denominator)
    return clean(result)


def remote(state):
    return apply_count(source()['remote'], state, 60)


@lru_cache(maxsize=1)
def resolvent():
    """Exact CP encoding R=p(I-A)^-1 in connected matrix-unit blocks.

    Only the inverse computation uses these blocks. Positive input/output
    claims follow from the geometric sum of actual CP remote maps.
    """
    transitions = source()['remote']
    unseen = {(i, j) for i in range(66) for j in range(66)}
    blocks, lookup, digest = [], {}, hashlib.sha256()
    while unseen:
        first = min(unseen)
        reached, pending = {first}, [first]
        while pending:
            pair = pending.pop()
            for other in transitions.get(pair, {}):
                if other not in reached:
                    reached.add(other)
                    pending.append(other)
        unseen.difference_update(reached)
        pairs = tuple(sorted(reached))
        position = {pair: index for index, pair in enumerate(pairs)}
        n = len(pairs)
        d = sp.eye(n)*60
        for column, pair in enumerate(pairs):
            for output, count in transitions.get(pair, {}).items():
                d[position[output], column] -= count
        if d != d.T:
            raise ArithmeticError('remote quantum transfer is not self-adjoint')
        inverse = d.inv(method='DM')*24
        if d*inverse != sp.eye(n)*24:
            raise ArithmeticError('quantum resolvent inversion failed')
        rows = tuple(tuple(F(value) for value in row) for row in inverse.tolist())
        block = len(blocks)
        lookup.update({pair: (block, index) for index, pair in enumerate(pairs)})
        blocks.append((pairs, rows))
        digest.update(canonical([pairs, [[str(x) for x in row] for row in rows]]).encode())
    return tuple(blocks), lookup, digest.hexdigest()


def lift(state):
    """The completely positive encoding R, extended by reference identity."""
    blocks, lookup, _ = resolvent()
    result = defaultdict(F)
    for (i, j, r, s), value in state.items():
        block, column = lookup[i, j]
        pairs, matrix = blocks[block]
        for row, (a, b) in enumerate(pairs):
            coefficient = matrix[row][column]
            if coefficient:
                result[a, b, r, s] += coefficient*value
    return clean(result)


def decode(state):
    """G=(I-A)/p: positive on the transported cone, not arbitrary matrices."""
    return {key: value/P for key, value in add(state, remote(state), -1).items()}


def interval(state, cap=2, *, effective=False):
    if type(cap) is not int or not 1 <= cap <= 3 or type(effective) is not bool:
        raise ValueError('cap one through three and Boolean representation required')
    waiting = state
    total = {}
    for _ in range(cap):
        total = add(total, waiting)
        waiting = remote(waiting)
    matrices = source()['incidents']
    if effective:
        # F_m,h = (I-A)^-1 B_m (I-A^h); no hidden reset or native-state oracle.
        difference = add(state, waiting, -1)
        marked = {label: {key: value/P for key, value in
                    lift(apply_count(matrix, difference, 60)).items()}
                  for label, matrix in matrices.items()}
    else:
        marked = {label: apply_count(matrix, total, 60)
                  for label, matrix in matrices.items()}
    return marked, waiting


@lru_cache(maxsize=4)
def apparatus(gamma=1, pair=(1, 5)):
    """Declared occupation interaction followed by a rational exchange rotation.

    Columns describe U=R_pair (-1)^(gamma*n0*n1). The interaction is supplied
    apparatus, not a new native repair law or a derived Hamiltonian interval.
    """
    if type(gamma) is not int or gamma not in (0, 1) or pair not in ((1, 5), (1, 2)):
        raise ValueError('binary interaction switch and declared pair required')
    a, b = pair
    result = {}
    for i, state in enumerate(STATES):
        phase = F((-1)**(gamma*state[0]*state[1]))
        if state[a] == state[b]:
            result[i] = ((i, phase),)
        else:
            swapped = list(state)
            swapped[a], swapped[b] = swapped[b], swapped[a]
            j = INDEX[tuple(swapped)]
            result[i] = ((i, phase*F(3, 5)),
                         (j, phase*F(4 if state[a] else -4, 5)))
    return result


def conjugate(columns, state):
    result = defaultdict(F)
    for (i, j, r, s), value in state.items():
        for (a, x), (b, y) in product(columns.get(i, ()), columns.get(j, ())):
            result[a, b, r, s] += value*x*y  # all retained coefficients are real
    return clean(result)


def compose_columns(left, right):
    """Sparse ket-operator composition, used to check actual service locality."""
    result = {}
    for column, values in right.items():
        outputs = defaultdict(F)
        for middle, x in values:
            for row, y in left.get(middle, ()):
                outputs[row] += x*y
        if outputs := clean(outputs):
            result[column] = tuple(sorted(outputs.items()))
    return result


@lru_cache(maxsize=1)
def structural_certificate():
    """All operator units and all local apparatus/remote Kraus commutators."""
    data = source()
    checked = 0
    for sector in mixing.sectors(binary=True):
        pairs = tuple((quantum.FULL_TO_BINARY[i], quantum.FULL_TO_BINARY[j])
                      for i, j in sector.pairs)
        positions = {pair: i for i, pair in enumerate(pairs)}
        expected = mixing.binary_reduction('C', sector.difference)
        for column, pair in enumerate(pairs):
            counts = Counter(data['remote'].get(pair, {}))
            for matrix in data['incidents'].values():
                counts.update(matrix.get(pair, {}))
            if any(F(counts.get(output, 0), 60) != expected[row][column]
                   for output, row in positions.items()):
                raise ArithmeticError('marked quantum source differs from the analytic native reduction')
            checked += 1
    identity = {i: ((i, F(1)),) for i in range(66)}
    services = [apparatus(gamma) for gamma in (0, 1)]
    for gate in services:
        transpose = defaultdict(list)
        for column, values in gate.items():
            for row, value in values:
                transpose[row].append((column, value))
        if compose_columns(dict(transpose), gate) != identity:
            raise ArithmeticError('the declared apparatus is not exactly unitary')
    services += [{i: ((i, F(1)),) for i, state in enumerate(STATES) if state[1] == bit}
                 for bit in (0, 1)]
    for seam in ((0, 1), (0, 5)):
        services += [{i: ((j, F(1)),) for i, j in columns.items()}
                     for _, columns in data['branches'][seam]]
    commutators = 0
    for seam, branches in data['branches'].items():
        if set(seam).intersection(UNION):
            continue
        for _, columns in branches:
            operator = {i: ((j, F(1)),) for i, j in columns.items()}
            for service in services:
                if compose_columns(operator, service) != compose_columns(service, operator):
                    raise ArithmeticError('an admitted local service fails remote commutation')
                commutators += 1
    # A one-body diagonal phase exp(i sum_v theta_v n_v) has rectangle ratio +1.
    # The supplied n0*n1 phase instead gives -1 on this full-sector rectangle.
    rectangle = ((0, 1), (2, 5), (0, 5), (1, 2))
    phases = [(-1)**(STATES[address(ports)][0]*STATES[address(ports)][1]) for ports in rectangle]
    ratio = phases[0]*phases[1]*phases[2]*phases[3]
    if ratio != -1:
        raise ArithmeticError('the declared two-occupation interaction became one-body')
    return {'matrix_units_checked_against_analytic_native_reduction': checked,
            'matrix_units_checked_for_remote_and_incident_effects': 4356,
            'remote_trace_on_matrix_unit': '(3/5) delta_ij',
            'incident_trace_on_matrix_unit': '(2/5) delta_ij',
            'exact_remote_Kraus_local_service_commutators': commutators,
            'apparatus_unitarity_checked_on_all_66_input_columns': True,
            'occupation_phase_rectangle': rectangle,
            'occupation_phase_rectangle_ratio': str(ratio),
            'one_body_phase_rectangle_ratio': '1',
            'interaction_is_supplied_apparatus': True}


def detector(state, bit):
    return {key: value for key, value in state.items()
            if STATES[key[0]][1] == bit == STATES[key[1]][1]}


def execute(gamma=1, cap=2, reference=False, effective=False, *, initial=None):
    """Complete source/reference matrices for two capped overlapping ticks.

    Output key: (word, status, detector_bit). A word retains scheduled and
    feedback marks; the coherent apparatus is charged but adds no native tick.
    The optional initial operator permits linear matrix-unit verification.
    """
    if type(reference) is not bool or type(effective) is not bool:
        raise ValueError('explicit Boolean reference and representation choices required')
    gate = apparatus(gamma)
    initial = preparation(reference) if initial is None else dict(initial)
    initial = lift(initial) if effective else initial
    result = {}

    def finish(word, status, state):
        for bit in (0, 1):
            output = detector(state, bit)
            result[word, status, bit] = decode(output) if effective else output

    first, timeout = interval(initial, cap, effective=effective)
    finish((), 'timeout_first', timeout)
    for label, state in sorted(first.items()):
        prefix = (('scheduled', label),)
        seam = (0, 1) if label[1][0] > 0 else (0, 5)
        for feedback, matrix in sorted(source()['feedback'][seam].items()):
            controlled = conjugate(gate, apply_count(matrix, state, 2))
            word = prefix+(('feedback', feedback),)
            second, timeout = interval(controlled, cap, effective=effective)
            finish(word, 'timeout_second', timeout)
            for last, terminal in sorted(second.items()):
                finish(word+(('scheduled', last),), 'complete', terminal)
    return result


def output_digest(history):
    return hashlib.sha256(canonical([[key, [[entry, str(value)] for entry, value in sorted(state.items())]]
        for key, state in sorted(history.items())]).encode()).hexdigest()


def compare(left, right):
    if set(left) != set(right):
        raise ArithmeticError('complete output menu changed')
    for key in left:
        if left[key] != right[key]:
            raise ArithmeticError(f'conditional source/reference output changed at {key}')
    mass = sum((trace(state) for state in left.values()), F(0))
    if mass != 1:
        raise ArithmeticError('complete output does not preserve probability')
    timeout = sum((trace(state) for (_, status, _), state in left.items()
                   if status.startswith('timeout')), F(0))
    return {'complete_output_words': len(left),
            'positive_probability_words': sum(trace(state) > 0 for state in left.values()),
            'total_probability': str(mass), 'timeout_probability': str(timeout),
            'exact_complete_output_sha256': output_digest(left)}


def cone_control():
    """An exact positive effective input outside the admitted encoded cone."""
    i = address((0, 2))
    state = {(i, i, 0, 0): F(1)}
    decoded = decode(state)
    negative = [value for (a, b, _, _), value in decoded.items() if a == b and value < 0]
    marked, _ = interval(state, 1, effective=True)
    label = ((0, 1), (0, 0, 1))
    j = address((0, 5))
    value = marked[label].get((j, j, 0, 0), F(0))
    if not negative or value != -F(1, 720):
        raise ArithmeticError('transported-cone negative control did not discriminate')
    return {'positive_input_occupied_ports': [0, 2], 'positive_input_basis_index': i,
            'decoder_negative_population': str(min(negative)),
            'cap_one_mark': label, 'cap_one_output_population_index': j,
            'cap_one_output_occupied_ports': [0, 5],
            'cap_one_negative_population': str(value),
            'ordinary_positive_map_on_all_effective_density_matrices': False}


def boundary_control():
    """A cross-boundary operation cannot keep the untransformed dictionary."""
    gate = apparatus(0, (1, 2))
    for i in range(66):
        state = {(i, i, 0, 0): F(1)}
        left, right = remote(conjugate(gate, state)), conjugate(gate, remote(state))
        difference = add(left, right, -1)
        if difference:
            key = min(difference)
            return {'input_basis_index': i, 'matrix_entry': key,
                    'exact_AU_minus_UA_entry': str(difference[key]),
                    'untransformed_boundary_crossing_apparatus_admitted': False}
    raise ArithmeticError('cross-boundary intervention unexpectedly commutes')


def entangled_reference_certificate():
    initial = entangled_preparation()
    direct = execute(initial=initial, reference=True)
    effective = execute(initial=initial, reference=True, effective=True)
    result = compare(direct, effective)
    reference = defaultdict(F)
    for state in direct.values():
        for (i, j, r, s), value in state.items():
            if i == j:
                reference[r, s] += value
    if clean(reference) != {(0, 0): F(1, 2), (1, 1): F(1, 2)}:
        raise ArithmeticError('the untouched reference marginal changed')
    state = direct[(), 'timeout_first', 0]
    a, b = address((0, 2)), address((0, 3))
    witness = (state.get((a, a, 1, 1), 0)+state.get((b, b, 0, 0), 0)
               -state.get((a, b, 0, 1), 0)-state.get((b, a, 1, 0), 0))/2
    if witness >= 0:
        raise ArithmeticError('the retained timeout entanglement witness failed')
    result.update({'input': '(|ports0,2>|0> + |ports0,3>|1>)/sqrt(2)',
        'untouched_reference_marginal': 'I2/2',
        'timeout_first_detector_zero_probability': str(trace(state)),
        'partial_transpose_vector': '(|ports0,2>|1> - |ports0,3>|0>)/sqrt(2)',
        'unnormalized_partial_transpose_expectation': str(witness),
        'conditional_partial_transpose_expectation': str(witness/trace(state)),
        'entanglement_survives_in_retained_timeout': True})
    return result


def run():
    from . import clock_response

    blocks, _, digest = resolvent()
    comparisons = {}
    for reference, gamma in product((False, True), (0, 1)):
        direct = execute(gamma, reference=reference)
        effective = execute(gamma, reference=reference, effective=True)
        result = compare(direct, effective)
        if result['timeout_probability'] != '369/625':
            raise ArithmeticError('native timeout mass changed')
        comparisons[f'reference_{int(reference)}_interaction_{gamma}'] = result
    receipt = {'schema': 'oph.quantum_observer_correspondence.v1',
        'source_sha256': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                          for path in SOURCE_PATHS}),
        'contract': {'native_law': 'C', 'binary_total_two_dimension': 66,
            'regions': REGIONS, 'union': UNION, 'scheduler': 'uniform thirty seams and fair native coin',
            'attempt_cap_per_interval': 2, 'scheduled_union_incidents_requested': 2,
            'coherent_apparatus': 'R_(1,5) (-1)^(gamma*n0*n1), supplied',
            'initial_field_disclosure': False, 'terminal_detector': 'complete n1 projection',
            'same_initial_disclosure_contract_as_diagonal_5768_word_receipt': False,
            'untouched_reference_dimension_executed': 2,
            'parameter': 'native attempts and the two incident clocks, not physical duration'},
        'declared_resources': {'maximum_scheduled_native_attempts': 4,
            'maximum_additional_native_feedback_repairs': 1,
            'ideal_phase_operations_after_first_completion_including_identity_control': 1,
            'ideal_exchange_rotations_after_first_completion': 1,
            'terminal_complete_measurements_on_every_branch': 1,
            'expected_scheduled_native_attempts': '328/125',
            'expected_native_feedback_repairs': '16/25',
            'apparatus_contracted_from_emitted_source_pulses': False,
            'apparatus_physical_service_error_certificate': False},
        'exact_quantum_resolvent': {'matrix_units': 66**2,
            'connected_block_dimensions': {str(k): v for k, v in sorted(Counter(len(p) for p, _ in blocks).items())},
            'resolvent_sha256': digest, 'remote_effect': '(3/5) I',
            'encoding': 'R=(2/5)(I-A)^-1, CPTP', 'decoder': 'G=(I-A)/(2/5)',
            'matrix_cone_with_reference': '(R tensor id)(positive source/reference operators)',
            'finite_cap_map': '(I-A)^-1 B_m (I-A^h)',
            'uniform_two_interval_timeout_bound': '1-[1-(3/5)^h]^2'},
        'all_input_structural_certificate': structural_certificate(),
        'complete_coherent_reference_comparisons': comparisons,
        'surviving_reference_entanglement': entangled_reference_certificate(),
        'independent_joint_clock_response': clock_response.audit(),
        'transported_cone_control': cone_control(),
        'boundary_intervention_control': boundary_control(),
        'acceptance': {'arbitrary_reference_extension_of_exact_map_identity': True,
            'number_preserving_union_interventions': True,
            'arbitrary_inaccessible_or_cross_boundary_interventions': False,
            'derived_interacting_Lorentzian_QFT_or_physical_clock': False,
            'physical_action_or_interaction_selected_by_native_repair': False,
            'observed_Universe_identification': False, 'federation_required': False}}
    return strict_json.loads(canonical(receipt))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    destination = parser.add_mutually_exclusive_group()
    destination.add_argument('--check', type=Path)
    destination.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('--output requires a new destination')
    receipt = run()
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text()), receipt)
        print('quantum observer correspondence receipt matches')
    elif args.output:
        with args.output.open('x') as stream:
            stream.write(canonical(receipt))
    else:
        print(canonical(receipt))
