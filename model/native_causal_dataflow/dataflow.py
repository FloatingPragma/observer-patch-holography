"""Classical-branch provenance for the actual bounded native interpreter.

All raw branches are retained. Authenticated value logs are constructed only
for branches with nonzero probability from the supplied basis preparation.
The current single FIFO makes full storage provenance sequential; this is
not a claim that the physical spacetime order is a chain.
"""
from copy import deepcopy
from dataclasses import asdict
from fractions import Fraction as F
import hashlib
import json

from ..native_program import interpreter as program
from ..observer import native


SOURCE_FILES = (
    'model/native_causal_dataflow/__init__.py',
    'model/native_causal_dataflow/dataflow.py',
    'model/native_causal_dataflow/test_dataflow.py',
    'model/native_causal_dataflow/README.md',
    'model/native_program/interpreter.py',
    'model/local_repair_family/local_bias.py',
    'model/observer/native.py',
    'model/observer/source_discrimination/circuit.py',
    'Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean',
    'Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean',
)


def _initial(value, loads):
    value = program.admit(value)
    if value['carrier']['kind'] != 'local':
        raise ValueError('this provenance adapter admits one local native carrier')
    if (type(loads) not in (tuple, list) or len(loads) != 12
            or any(type(x) is not int or x < 0 for x in loads) or sum(loads) != 2):
        raise ValueError('a classical total-two native basis preparation is required')
    basis = program.local.basis(value['carrier']['shell'] == 'H66')
    try:
        index = basis.index(tuple(loads))
    except ValueError as error:
        raise ValueError('the preparation is outside the admitted native shell') from error
    return value, tuple(loads), basis, index


def _commands(value):
    if sum(command['repeat'] for command in value['commands']) > program.MAX_EVALUATED_STEPS:
        raise ValueError('the provenance adapter requires a bounded fully executed program')
    return tuple((command, repetition) for command in value['commands']
                 for repetition in range(command['repeat']))


def _record_value(record):
    return (record.kind, record.selector, record.seam, record.law, record.coin,
            record.raw, record.label, record.current)


def independent(first, second):
    """The exact write/write and write/support predicate used by Lean."""
    a, b = set(first['write_set']), set(second['write_set'])
    return not (a & b or a & set(second['causal_support'])
                or b & set(first['causal_support']))


def replay_trace(value, loads, records):
    """Authenticate all values and writers against actual interpreter branches.

    Supplied records are not trusted as a proof of reachability. Every next
    record must have positive exact weight in `_outcomes` on this actual
    basis state and FIFO. This reuses the interpreter's selectors, mixtures,
    coin weights, padding and overflow behavior without rebuilding a law.
    """
    value, loads, basis, index = _initial(value, loads)
    commands = _commands(value)
    if type(records) not in (list, tuple) or len(records) != len(commands):
        raise ValueError('a complete nonterminal record per command repetition is required')
    values = {f'port:{i}': x for i, x in enumerate(loads)}
    values.update(fifo=(), archive=())
    writers = {name: None for name in values}
    density, memory, archive = {(index, 0, index, 0): F(1)}, (), ()
    rows = []
    for event_id, ((command, repetition), record) in enumerate(zip(commands, records, strict=True)):
        if not isinstance(record, program.Record) or record.kind not in ('hold', 'repair'):
            raise ValueError('only completed native repair/hold records are admitted')
        choices = program._outcomes(value, command, density, memory, archive, repetition=repetition)
        matches = [outcome for outcome in choices if outcome.record == record and program.trace(outcome.density) > 0]
        if len(matches) != 1:
            raise ValueError('record is not one unique positive-weight interpreter branch')
        outcome, = matches
        mass = program.trace(outcome.density)
        normalized = {key: item/mass for key, item in outcome.density.items()}
        before = dict(values)
        before_writers = dict(writers)
        ports = () if record.kind == 'hold' else tuple(f'port:{p}' for p in record.seam)
        # `_commit` constructs the archive tuple on every call, even while
        # it remains empty; only overflow changes/stamps that register.
        reads, writes = set(ports)|{'fifo', 'archive'}, set(ports)|{'fifo'}
        overflow = len(memory) == value['memory']['capacity']
        if overflow:
            reads.add('archive'); writes.add('archive')
        control_read = command['op'] == 'repair' and command['selector']['kind'] in ('last_current', 'observer')
        # `_choices` is evaluated for the repair menu even on its hold branch.
        # This is an executed-read support, not a claim of minimal influence.
        if record.kind == 'repair':
            after_loads, raw = native.native(loads, record.seam, record.coin)
            if raw[1] != record.current:
                raise ValueError('current record differs from the actual endpoint transition')
        else:
            after_loads = loads
        next_index = basis.index(after_loads)
        if normalized != {(next_index, 0, next_index, 0): F(1)}:
            raise ValueError('recorded branch is not the replayed native basis output')
        for register in ports:
            values[register] = after_loads[int(register.split(':')[1])]
        values['fifo'] = tuple(_record_value(r) for r in outcome.memory)
        values['archive'] = tuple(_record_value(r) for r in outcome.archive)
        endpoint_parents = sorted({writers[name] for name in ports if writers[name] is not None})
        controller_parents = ([] if not control_read or writers['fifo'] is None else [writers['fifo']])
        storage_parents = sorted({writers[name] for name in reads if name in ('fifo', 'archive')
                                  and writers[name] is not None})
        parents = sorted(set(endpoint_parents+controller_parents+storage_parents))
        for parent in parents:
            if not parent < event_id:
                raise ArithmeticError('authenticated parent is not earlier in the emitted run')
            if not any(name in rows[parent]['write_set'] and before_writers[name] == parent
                       and before[name] == rows[parent]['after_values'][name] for name in reads):
                raise ArithmeticError('parent has no authenticated read-after-write witness')
        for register in writes:
            writers[register] = event_id
        ancestors = sorted(set(parents).union(*(set(rows[parent]['ancestors']) for parent in parents)))
        height = 0 if not parents else 1+max(rows[parent]['height'] for parent in parents)
        rows.append({'event_id': event_id, 'record': asdict(record),
                     'before_values': before, 'after_values': dict(values),
                     'before_writers': before_writers, 'after_writers': dict(writers),
                     'read_set': sorted(reads), 'write_set': sorted(writes), 'causal_support': sorted(reads),
                     'endpoint_parents': endpoint_parents, 'controller_parents': controller_parents,
                     'storage_parents': storage_parents, 'parents': parents,
                     'controller_last_current': memory[-1].current if control_read and memory else None,
                     'controller_read_executed': control_read,
                     'ancestors': ancestors, 'height': height})
        density, memory, archive, loads = normalized, outcome.memory, outcome.archive, after_loads
    return rows


def verify_trace(value, loads, records, rows):
    if replay_trace(value, loads, records) != rows:
        raise ValueError('declared provenance differs from authenticated execution')
    return True


def execute(value, loads):
    """Return every actual raw interpreter branch, including all zero outcomes."""
    value, loads, basis, index = _initial(value, loads)
    history = program.execute(value, {(index, 0, index, 0): F(1)})
    if history.status != 'COMPLETE':
        raise ValueError('the provenance theorem requires a complete bounded native execution')
    rows = []
    for branch in history.branches:
        records = branch.archive+branch.memory
        probability = program.trace(branch.density)
        events = replay_trace(value, loads, records) if probability else None
        final = {i: coefficient for (i, a, j, b), coefficient in branch.density.items()
                 if i == j and a == b}
        rows.append({'records': records, 'probability': probability,
                     'events': events, 'final_populations': final})
    if sum(row['probability'] for row in rows) != 1:
        raise ArithmeticError('complete raw branches do not normalize')
    return {'program': value, 'loads': loads, 'basis': basis, 'branches': rows,
            'complete_raw_branch_count': len(rows),
            'positive_branch_count': sum(row['probability'] > 0 for row in rows)}


def receiver_probability(result, port):
    """Probability that the specified final native port is nonempty."""
    if type(port) is not int or not 0 <= port < 12:
        raise ValueError('one native port is required')
    return sum(coefficient*int(result['basis'][i][port] > 0) for row in result['branches']
               for i, coefficient in row['final_populations'].items())


def feedback_program(*, adaptive=True, capacity=1):
    if type(adaptive) is not bool:
        raise ValueError('adaptive must be an explicit Boolean')
    selector = ({'kind': 'last_current', 'positive': [2, 3], 'otherwise': [2, 10]}
                if adaptive else {'kind': 'fixed', 'seam': [2, 3]})
    return program.program(biased=False, capacity=capacity, commands=[
        {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [0, 1]}, 'law': 'C', 'repeat': 1},
        {'op': 'repair', 'selector': selector, 'law': 'C', 'repeat': 1}])


def certificate():
    value = feedback_program()
    inputs = [tuple(int(i in occupied) for i in range(12)) for occupied in ((0, 2), (1, 2))]
    results = [execute(value, loads) for loads in inputs]
    probabilities = tuple(receiver_probability(result, 3) for result in results)
    if probabilities != (F(1, 4), F()):
        raise ArithmeticError('native controller influence witness changed')
    for result in results:
        for row in result['branches']:
            if not row['probability']:
                continue
            first, second = row['events']
            if (second['endpoint_parents'] or second['controller_parents'] != [0]
                    or second['storage_parents'] != [0] or second['ancestors'] != [0]
                    or independent(first, second)):
                raise ArithmeticError('controller/storage causal edge was lost')
    raw = [[str(row['probability']) for row in result['branches']] for result in results]
    return {'schema': 'oph.native_causal_dataflow.v1', 'program': deepcopy(value),
            'preparations': [list(loads) for loads in inputs], 'receiver': 3,
            'receiver_probabilities': [str(x) for x in probabilities],
            'complete_raw_branch_counts': [r['complete_raw_branch_count'] for r in results],
            'positive_branch_counts': [r['positive_branch_count'] for r in results],
            'complete_probability_vectors_sha256': hashlib.sha256(json.dumps(raw).encode()).hexdigest(),
            'positive_branches_have_no_endpoint_parent_at_second_event': True,
            'positive_branches_have_authenticated_controller_and_storage_parent': True,
            'full_single_FIFO_dataflow_is_a_chain': True,
            'full_footprint_is_operational_dependency_overapproximation': True,
            'adaptive_receiver_probability_witnesses_influence': True,
            'storage_chain_alone_excludes_physical_Lorentzian_limit': False,
            'physical_spacetime_order_or_Lorentz_cone_derived': False,
            'quantum_coherent_programs_or_federations_covered': False}
