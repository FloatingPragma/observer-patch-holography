"""Admission and complete instruments for a bounded native program grammar.

Admission fixes the operations before a state is supplied. A numerical
evaluator resource stop is an identity channel with a retained stop record;
it is not exclusion of a member of the declared federation growth grammar.
The six-level coherent service has its own explicitly supplied register.
No interpreter operation accepts a field Hamiltonian, a target mass or a
callback that inspects the native state.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math

import numpy as np

from ..local_repair_family import local_bias as local
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..federated_observer import source as federation_source


SCHEMA = 'oph.native_program.v1'
WORD_SCHEMA = 'oph.native_program.source_word.v1'
MAX_NODES = 4096
MAX_REPEAT = 2**32
MAX_ATTEMPTS = 2**32
MAX_WORD = 65536
MAX_WORD_REPEATS = 2**32
MAX_DESCRIPTION_BYTES = 2**24
MAX_EVALUATED_LEVEL = 2
MAX_EVALUATED_BRANCHES = 65536
MAX_EVALUATED_STEPS = 64
MAX_EVALUATED_QUBITS = 6
MAX_EVALUATED_WORD_PULSES = 256


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def _keys(value, names):
    if type(value) is not dict or set(value) != set(names):
        raise ValueError('exact declared fields required: '+','.join(names))


def _integer(value, lower=0, upper=None):
    if type(value) is not int or value < lower or (upper is not None and value > upper):
        raise ValueError('integer outside the declared resource/domain boundary')
    return value


def _rational(value, *, probability=False):
    if type(value) is not str or len(value) > 700:
        raise ValueError('a bounded canonical rational string required')
    try:
        result = F(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('invalid rational') from error
    if str(result) != value:
        raise ValueError('canonical rational spelling required')
    if probability and (not 0 <= result <= 1 or result.denominator > 256):
        raise ValueError('probability outside the fixed denominator boundary')
    return result


def carrier(value):
    """Finite state/algebra declaration, independently of repair compatibility."""
    if type(value) is not dict:
        raise ValueError('typed carrier required')
    if value.get('kind') == 'local':
        _keys(value, ('kind', 'shell'))
        if value['shell'] not in ('H78', 'H66'):
            raise ValueError('the existing total-two local shell required')
    elif value.get('kind') == 'independent_native_product':
        _keys(value, ('kind', 'shell', 'cells'))
        if value['shell'] not in ('H78', 'H66'):
            raise ValueError('the existing total-two local shell required')
        _integer(value['cells'], 1)
    elif value.get('kind') == 'federation':
        _keys(value, ('kind', 'level', 'total'))
        level, total = _integer(value['level']), _integer(value['total'])
        # Compare with 240*4**level without allocating that enormous integer.
        quotient = total >> (2*level)
        if quotient > 240 or (quotient == 240 and total != (quotient << (2*level))):
            raise ValueError('global binary occupation exceeds the port count')
    else:
        raise ValueError('only native local or geodesic port_pair carriers admitted')
    return deepcopy(value)


def carrier_interface(value):
    c = carrier(value)
    if c['kind'] in ('local', 'independent_native_product'):
        return {'state': 'density on '+c['shell']+' tensor power '+str(c.get('cells', 1)),
                'algebra': 'M'+c['shell'][1:]+' tensor power '+str(c.get('cells', 1)),
                'compatibility': 'nonnegative integer loads sum to two; spectators fixed; repaired endpoint difference at most one',
                'ports': 12*c.get('cells', 1),
                'refinement': 'append a declared local preparation; partial trace forgets appended factors; H66 is an invariant shell of H78',
                'coherent_coupling': 'a separate counted service needs a bound proper-code encoding; no automatic native load is asserted'}
    level = c['level']
    return {'state': 'density on the full global binary fixed-number sector',
            'algebra': 'M_binomial(240*4^level,total)', 'total': c['total'],
            'ports_formula': '240*4^level', 'level': level,
            'compatibility': 'binary occupation, fixed global total, unchanged spectators, native endpoint conservation',
            'graph_grammar': 'oph_exact.federation.build_federation(level,port_pair)',
            'refinement': 'cell ancestry with unchanged local port address; coordinate R/J is separate from occupation partial trace',
            'dynamic_intertwining_claimed': False,
            'local_total_two_assumed': False}


def _parameters(value, c):
    _keys(value, ('port_weights', 'seam_weights', 'c_weight', 'repair_probability', 'exploration_probability'))
    if type(value['port_weights']) is not list or type(value['seam_weights']) is not list:
        raise ValueError('literal weight lists required')
    p = local.Parameters(tuple(value['port_weights']), tuple(value['seam_weights']),
                         _rational(value['c_weight'], probability=True),
                         _rational(value['repair_probability'], probability=True))
    epsilon = _rational(value['exploration_probability'], probability=True)
    if epsilon < F(1, 16):
        raise ValueError('positive exploration floor required')
    if c['kind'] == 'federation' and (p.port_weights != (1,)*12 or p.seam_weights != (1,)*30):
        raise ValueError('federation admits the existing fair native law, not an invented bias lift')
    return p, epsilon


def _seam(value, c):
    if type(value) is not list or len(value) != 2:
        raise ValueError('literal ordered endpoint pair required')
    a, b = (_integer(x) for x in value)
    if a >= b:
        raise ValueError('ordered distinct endpoints required')
    if c['kind'] in ('local', 'independent_native_product'):
        if (a, b) not in local.SEAMS:
            raise ValueError('seam is not in the native graph')
    elif c['level'] <= MAX_EVALUATED_LEVEL:
        if (a, b) not in federation_source.federation(c['level']).seams:
            raise ValueError('seam is not in the declared federation')
    else:
        # A finite native address may only use a known local cell seam at
        # unevaluated levels. Intercell selection uses the whole graph menu.
        if a//12 != b//12 or (a % 12, b % 12) not in local.SEAMS:
            raise ValueError('unevaluated intercell addresses require the uniform graph selector')
        if b >> (2*c['level']) >= 240:
            raise ValueError('endpoint outside the growth carrier')
    return a, b


def _selector(value, c):
    if type(value) is not dict:
        raise ValueError('literal record-only selector required')
    kind = value.get('kind')
    if kind == 'uniform':
        _keys(value, ('kind',))
    elif kind == 'fixed':
        _keys(value, ('kind', 'seam')); _seam(value['seam'], c)
    elif kind in ('last_current', 'observer'):
        _keys(value, ('kind', 'positive', 'otherwise'))
        _seam(value['positive'], c); _seam(value['otherwise'], c)
    else:
        raise ValueError('selector may read only the last committed current')


def _native_cost(value):
    attempts = sum(command['repeat'] for command in value['commands'])
    repair = sum(command['repeat'] for command in value['commands'] if command['op'] == 'repair')
    capacity = value['memory']['capacity']
    c = value['carrier']
    address_bits = (8+2*c['level'] if c['kind'] == 'federation'
                    else (12*c.get('cells', 1)-1).bit_length())
    record_bits = 32+2*address_bits
    return {'attempts_upper': attempts, 'native_repairs_upper': repair,
            'controller_reads_upper': repair,
            'committed_record_writes_upper': attempts,
            'exported_records_upper': max(0, attempts-capacity) if value['memory']['overflow'] == 'export_oldest' else 0,
            'FIFO_slots': capacity, 'controller_states': 2,
            'controller_state_bits': 1, 'controller_initializations': 1,
            'program_description_bytes': len(canonical({k: v for k, v in value.items()
                                                       if k != 'declared_resources'}).encode()),
            'source_initializations': 1, 'classical_blank_FIFO_initializations': capacity,
            'fresh_coin_draws_upper': repair,
            'selector_draws_upper': repair,
            'law_tag_draws_upper': repair,
            'hold_draws_upper': repair,
            'terminal_status_records': 1,
            'retained_record_capacity_upper': attempts+1,
            'record_bits_per_slot_upper': record_bits,
            'record_storage_bits_upper': (attempts+1)*record_bits,
            'AST_nodes': len(value['commands'])}


def admit(value):
    """Validate and detach a complete literal program; never inspect input state."""
    _keys(value, ('schema', 'carrier', 'parameters', 'memory', 'commands', 'declared_resources'))
    if value['schema'] != SCHEMA or len(canonical(value).encode()) > MAX_DESCRIPTION_BYTES:
        raise ValueError('unknown or oversized native program description')
    c = carrier(value['carrier']); _parameters(value['parameters'], c)
    _keys(value['memory'], ('capacity', 'overflow'))
    _integer(value['memory']['capacity'], 1, 8)
    if value['memory']['overflow'] not in ('export_oldest', 'halt'):
        raise ValueError('complete export or halt policy required')
    commands = value['commands']
    if type(commands) is not list or not 1 <= len(commands) <= MAX_NODES:
        raise ValueError('bounded literal command list required')
    for command in commands:
        if type(command) is not dict:
            raise ValueError('literal command required')
        if command.get('op') == 'repair':
            _keys(command, ('op', 'selector', 'law', 'repeat')+
                  (('cell',) if c['kind'] == 'independent_native_product' else ()))
            if c['kind'] == 'independent_native_product':
                _integer(command['cell'], 0, c['cells']-1)
            _selector(command['selector'], c)
            if command['law'] not in ('C', 'E', 'mixture', 'alternating'):
                raise ValueError('declared native record partition required')
        elif command.get('op') == 'hold':
            _keys(command, ('op', 'repeat'))
        else:
            raise ValueError('undeclared native operation')
        _integer(command['repeat'], 1, MAX_REPEAT)
    expected = _native_cost(value)
    if expected['attempts_upper'] > MAX_ATTEMPTS:
        raise ValueError('program attempt ceiling exceeded')
    if expected['record_storage_bits_upper'] > 2**40:
        raise ValueError('retained record bit ceiling exceeded')
    if type(value['declared_resources']) is not dict or value['declared_resources'] != expected:
        raise ValueError('native resource declaration disagrees with recomputed cost')
    if any(type(x) is not int for x in value['declared_resources'].values()):
        raise ValueError('integer resource declarations required')
    return deepcopy(value)


def program(*, shell='H66', biased=True, commands=None, capacity=1, overflow='export_oldest', federation=None, cells=None):
    """Construct a declaration, without granting admission to supplied callbacks."""
    c = {'kind': 'local', 'shell': shell} if federation is None else {'kind': 'federation', **federation}
    if cells is not None:
        if federation is not None:
            raise ValueError('one carrier grammar per program required')
        c = {'kind': 'independent_native_product', 'shell': shell, 'cells': cells}
    p = local.Parameters(port_weights=((25, 25)+(8,)*10 if biased and federation is None else (1,)*12))
    value = {'schema': SCHEMA, 'carrier': c,
             'parameters': {**p.description(), 'exploration_probability': '1/2'},
             'memory': {'capacity': capacity, 'overflow': overflow},
             'commands': commands if commands is not None else [
                 {'op': 'repair', 'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': 1,
                  **({'cell': 0} if cells is not None else {})}]}
    value['declared_resources'] = _native_cost(value)
    return admit(value)


@dataclass(frozen=True)
class Record:
    kind: str
    selector: str = 'none'
    seam: tuple[int, int] | None = None
    law: str | None = None
    coin: int | None = None
    raw: int | None = None
    label: int | None = None
    current: int | None = None


@dataclass(frozen=True)
class Outcome:
    record: Record
    density: dict
    memory: tuple[Record, ...]
    archive: tuple[Record, ...]
    halted: bool = False


@dataclass(frozen=True)
class History:
    branches: tuple[Outcome, ...]
    attempted_steps: int
    program_digest: str
    status: str
    resources: dict


def _validate_density(value, c):
    if type(value) is not dict:
        raise ValueError('finite sparse source/reference operator required')
    result = {}
    for key, coefficient in value.items():
        if type(key) is not tuple or len(key) != 4:
            raise ValueError('operator key must be source,reference,source,reference')
        for x in key:
            _integer(x)
        for x in (key[0], key[2]):
            if c['kind'] in ('local', 'independent_native_product'):
                if x >= int(c['shell'][1:])**c.get('cells', 1):
                    raise ValueError('source index outside the declared shell')
            elif x.bit_count() != c['total'] or x.bit_length() > (240 << (2*c['level'])):
                raise ValueError('source mask outside the declared global sector')
        # Exact rationals keep zero outcomes distinct from underflow. Complex
        # operators are the complex-linear extension of two exact real calls.
        if type(coefficient) not in (int, F):
            raise ValueError('exact rational coefficients required; apply real and imaginary parts separately')
        if coefficient:
            result[key] = F(coefficient)
    return result


def trace(density):
    return sum((v for (i, a, j, b), v in density.items() if i == j and a == b), F())


def reference_marginal(densities):
    result = defaultdict(F)
    for density in densities:
        for (i, a, j, b), value in density.items():
            if i == j:
                result[a, b] += value
    return {key: value for key, value in result.items() if value}


def _choices(command, memory, p, epsilon, c):
    selector = command['selector']
    if c['kind'] in ('local', 'independent_native_product'):
        seams, weights = local.SEAMS, p.seam_weights
    else:
        seams = federation_source.federation(c['level']).seams
        weights = (1,)*len(seams)
    uniform = tuple(('uniform', edge, F(weight, sum(weights))) for edge, weight in zip(seams, weights, strict=True))
    if selector['kind'] == 'uniform':
        return uniform
    if selector['kind'] == 'fixed':
        return (('fixed', tuple(selector['seam']), F(1)),)
    positive = bool(memory and (memory[-1].current or 0) > 0)
    preferred = tuple(selector['positive' if positive else 'otherwise'])
    if selector['kind'] == 'last_current':
        return (('feedback', preferred, F(1)),)
    return tuple((tag, edge, epsilon*weight) for tag, edge, weight in uniform)+(('feedback', preferred, 1-epsilon),)


@lru_cache(maxsize=256)
def _local_mapping(shell, law, edge, coin, label):
    for c, value, pairs in local.raw_maps(law, edge, shell == 'H66'):
        if (c, value) == (coin, label):
            return dict(pairs)
    return {}


def _apply_native(density, c, law, edge, coin, label, weight, cell=0):
    if not weight or label is None:
        return {}
    if c['kind'] in ('local', 'independent_native_product'):
        base = _local_mapping(c['shell'], law, edge, coin, label)
        dimension = int(c['shell'][1:])
        scale = dimension**cell
        mapping = {i: i+(base[(i//scale) % dimension]-(i//scale) % dimension)*scale
                   for i in {key[k] for key in density for k in (0, 2)}
                   if (i//scale) % dimension in base}
    else:
        mapping = {}
        for mask in {key[i] for key in density for i in (0, 2)}:
            occupied = tuple(i for i in range(mask.bit_length()) if mask & (1 << i))
            image, current = federation_source.binary_step(occupied, edge, coin)
            raw = current if law == 'C' else int(edge[0] in occupied)-int(edge[1] in occupied)
            if raw == label:
                mapping[mask] = sum(1 << i for i in image)
    return {(mapping[i], a, mapping[j], b): weight*value
            for (i, a, j, b), value in density.items() if i in mapping and j in mapping}


def _commit(record, density, memory, archive, capacity):
    exported = memory[:1] if len(memory) == capacity else ()
    retained = memory[1:] if exported else memory
    return Outcome(record, density, (*retained, record), (*archive, *exported))


def _outcomes(value, command, density, memory=(), archive=(), *, repetition=0):
    """Every raw instrument block, including zero maps and zero law weights.

    This kernel defines every finite history even when expanded-tree execution
    hits its resource ceiling. Memory and archive contain committed records;
    neither a selector nor its controller receives ``density`` as an argument.
    """
    c, p = value['carrier'], _parameters(value['parameters'], value['carrier'])[0]
    capacity = value['memory']['capacity']
    if len(memory) == capacity and value['memory']['overflow'] == 'halt':
        return (Outcome(Record('FIFO_TIMEOUT'), dict(density), memory, archive, True),)
    if command['op'] == 'hold':
        return (_commit(Record('hold'), dict(density), memory, archive, capacity),)
    epsilon = _rational(value['parameters']['exploration_probability'], probability=True)
    output = [_commit(Record('hold'), {key: (1-p.repair_probability)*v for key, v in density.items()
                                     if 1-p.repair_probability}, memory, archive, capacity)]
    law = command['law']
    weight_c = p.c_weight if law == 'mixture' else F(law == 'C' or (law == 'alternating' and repetition % 2 == 0))
    for tag, edge, selection in _choices(command, memory, p, epsilon, c):
        for kind, mix in (('C', weight_c), ('E', 1-weight_c)):
            for coin in (0, 1):
                probability = local.coin_probability(p, edge, coin) if c['kind'] != 'federation' else F(1, 2)
                weight = p.repair_probability*selection*mix*probability
                decoding = circuit.record_values(kind, coin)
                for raw in range(1 << circuit.RECORD_BITS[kind]):
                    label = decoding.get(raw)
                    current = (None if label is None else label if kind == 'C'
                               else (label-(2*coin-1)*(abs(label) % 2))//2)
                    record = Record('repair', tag, edge, kind, coin, raw, label, current)
                    output.append(_commit(record, _apply_native(density, c, kind, edge, coin, label, weight,
                                                                  command.get('cell', 0)), memory, archive, capacity))
    return tuple(output)


def execute(value, density, *, max_branches=MAX_EVALUATED_BRANCHES, max_steps=MAX_EVALUATED_STEPS):
    value = admit(value)
    _integer(max_branches, 1, MAX_EVALUATED_BRANCHES); _integer(max_steps, 0, MAX_EVALUATED_STEPS)
    # Stop before dense/sparse source validation could materialize an enormous
    # graph. The supplied operator is returned verbatim and is not called a
    # valid density unless its domain has been evaluated.
    c = value['carrier']
    if c['kind'] == 'federation' and c['level'] > MAX_EVALUATED_LEVEL:
        return History((Outcome(Record('RESOURCE_STOP_UNEVALUATED_CARRIER'), deepcopy(density), (), (), True),),
                       0, digest(value), 'RESOURCE_STOP_UNEVALUATED_CARRIER', {'attempts_executed': 0})
    if c['kind'] == 'independent_native_product' and c['cells'] > 64:
        return History((Outcome(Record('RESOURCE_STOP_UNEVALUATED_CARRIER'), deepcopy(density), (), (), True),),
                       0, digest(value), 'RESOURCE_STOP_UNEVALUATED_CARRIER', {'attempts_executed': 0})
    initial = _validate_density(density, c)
    branches = (Outcome(Record('initial'), initial, (), ()),)
    steps, stopped = 0, False
    for command in value['commands']:
        for repetition in range(command['repeat']):
            # A retained intrinsic terminal outcome ends execution. Later
            # literal commands cannot consume attempts or replace that outcome
            # with a host resource stop.
            if all(branch.halted for branch in branches):
                break
            if steps >= max_steps:
                stopped = True; break
            sizes = []
            for branch in branches:
                if branch.halted or (len(branch.memory) == value['memory']['capacity'] and value['memory']['overflow'] == 'halt'):
                    sizes.append(1)
                elif command['op'] == 'hold':
                    sizes.append(1)
                else:
                    p, epsilon = _parameters(value['parameters'], c)
                    sizes.append(1+24*len(_choices(command, branch.memory, p, epsilon, c)))
            if sum(sizes) > max_branches:
                stopped = True; break
            branches = tuple(result for branch in branches for result in (
                (branch,) if branch.halted else _outcomes(value, command, branch.density, branch.memory, branch.archive, repetition=repetition)))
            steps += 1
        if stopped or all(branch.halted for branch in branches):
            break
    if stopped:
        branches = tuple(branch if branch.halted else Outcome(Record('RESOURCE_STOP'), branch.density,
                         branch.memory, branch.archive, True) for branch in branches)
    status = 'RESOURCE_STOP' if stopped else ('FIFO_TIMEOUT' if all(b.halted for b in branches) else 'COMPLETE')
    return History(branches, steps, digest(value), status,
                   {'attempts_executed': steps, 'complete_expanded_blocks': len(branches),
                    'maximum_FIFO_occupied': max(map(lambda b: len(b.memory), branches)),
                    'maximum_exported_records': max(map(lambda b: len(b.archive), branches)),
                    'terminal_status_records_per_path': 1})


def _word_cost(value):
    repeats, word = value['repeat'], value['word']
    pulses = repeats*len(word)
    return {'native_pulses': pulses, 'proper_code_load_events': pulses,
            'proper_code_export_events': pulses, 'six_level_blank_preparations': pulses,
            'export_complement_outcomes_per_run': repeats*sum(6-(4 if p['kind'] == 'J' else 2) for p in word),
            'coherent_work_qubits': value['qubits'], 'initial_work_qubits_supplied': value['qubits'],
            'terminal_quantum_readout_services': 1,
            'zero_angle_requests': repeats*sum(_rational(p['angle']) == 0 for p in word),
            'literal_word_nodes': len(word), 'word_repetitions': repeats,
            'program_description_bytes': len(canonical({k: v for k, v in value.items()
                                                       if k != 'declared_resources'}).encode()),
            'source_position_counter_bits': max(1, pulses.bit_length()),
            'source_status_bits': 2,
            'angle_payload_bits': sum(max(1, abs(_rational(p['angle']).numerator).bit_length())+
                                      _rational(p['angle']).denominator.bit_length() for p in word)}


def admit_word(value):
    _keys(value, ('schema', 'qubits', 'word', 'repeat', 'declared_resources'))
    if value['schema'] != WORD_SCHEMA or len(canonical(value).encode()) > MAX_DESCRIPTION_BYTES:
        raise ValueError('unknown or oversized coherent service word')
    _integer(value['qubits'], 1, 4096); _integer(value['repeat'], 1, MAX_WORD_REPEATS)
    if type(value['word']) is not list or not 1 <= len(value['word']) <= MAX_WORD:
        raise ValueError('bounded literal source word required')
    for pulse in value['word']:
        _keys(pulse, ('kind', 'targets', 'angle'))
        if pulse['kind'] not in compiler.KINDS:
            raise ValueError('undeclared future source primitive')
        targets = pulse['targets']
        if type(targets) is not list or len(targets) != (2 if pulse['kind'] == 'J' else 1):
            raise ValueError('source primitive arity mismatch')
        for target in targets:
            _integer(target, 0, value['qubits']-1)
        if len(set(targets)) != len(targets):
            raise ValueError('source targets must be distinct')
        angle = _rational(pulse['angle'])
        if (abs(angle) > 2**20 or angle.denominator > 2**1074
                or angle.denominator & (angle.denominator-1) or F(float(angle)) != angle):
            raise ValueError('bounded exact binary64 angle payload required')
    cost = _word_cost(value)
    if cost['native_pulses'] > 2**48:
        raise ValueError('source pulse ceiling exceeded')
    if value['declared_resources'] != cost or any(type(x) is not int for x in value['declared_resources'].values()):
        raise ValueError('source resource declaration differs from literal execution')
    return deepcopy(value)


def source_word(word, qubits, repeat=1):
    value = {'schema': WORD_SCHEMA, 'qubits': qubits, 'word': [
        {'kind': pulse['kind'], 'targets': list(pulse['targets']), 'angle': str(F(pulse['angle']))}
        for pulse in word], 'repeat': repeat}
    value['declared_resources'] = _word_cost(value)
    return admit_word(value)


def _lift(matrix, targets, qubits):
    """Existing compiler packing: one bit or |a,b> mapped to |2a+b>."""
    size = 1 << qubits
    result = np.zeros((size, size), complex)
    for column in range(size):
        local_column = sum(((column >> target) & 1) << (len(targets)-1-i) for i, target in enumerate(targets))
        cleared = column
        for target in targets:
            cleared &= ~(1 << target)
        for local_row in range(1 << len(targets)):
            row = cleared | sum(((local_row >> (len(targets)-1-i)) & 1) << target for i, target in enumerate(targets))
            result[row, column] = matrix[local_row, local_column]
    return result


def execute_word(value, density, *, reference_dimension=1):
    """Complete source-word instrument including first-abort complement labels.

    An export abort resets only its addressed proper code to zero, retains
    every spectator/reference and stops the word. The failed processor row
    is a retained Kraus label. No success normalization or polar repair occurs.
    Floating contractions are diagnostics; exact unitarity and complete
    export certify trace preservation independently of numerical residuals.
    """
    value = admit_word(value)
    _integer(reference_dimension, 1, 256)
    pulses = value['declared_resources']['native_pulses']
    if value['qubits'] > MAX_EVALUATED_QUBITS or pulses > MAX_EVALUATED_WORD_PULSES:
        return {'status': 'RESOURCE_STOP', 'blocks': {('RESOURCE_STOP',): np.array(density, copy=True)},
                'pulses_executed': 0, 'resources': value['declared_resources']}
    size = (1 << value['qubits'])*reference_dimension
    rho = np.asarray(density)
    if rho.shape != (size, size) or rho.dtype.kind not in 'iufc' or not np.isfinite(rho).all():
        raise ValueError('finite work/reference operator of declared dimensions required')
    rho = rho.astype(complex, copy=True)
    blocks, count = {}, 0
    for _ in range(value['repeat']):
        for pulse in value['word']:
            full = compiler.native_pulse(pulse['kind'], float(_rational(pulse['angle'])))
            d = 4 if pulse['kind'] == 'J' else 2
            for complement in range(d, 6):
                failed = np.zeros((d, d), complex)
                failed[0] = full[complement, :d]
                lifted = np.kron(_lift(failed, pulse['targets'], value['qubits']), np.eye(reference_dimension))
                blocks[('EXPORT_ABORT', count, complement)] = lifted@rho@lifted.conj().T
            good = np.kron(_lift(full[:d, :d], pulse['targets'], value['qubits']), np.eye(reference_dimension))
            rho = good@rho@good.conj().T
            count += 1
    blocks[('SUCCESS',)] = rho
    return {'status': 'COMPLETE', 'blocks': blocks, 'pulses_executed': count,
            'resources': value['declared_resources']}


def submit_word(value, density, *, reference_dimension=1):
    """Reject a malformed word before any pulse, retaining the original work.

    The returned rejection block is an identity on the submitted operator;
    rejection does not establish that an unvalidated input is a density.
    """
    try:
        admitted = admit_word(value)
    except (ValueError, TypeError, OverflowError) as error:
        return {'status': 'REJECTED_PROGRAM', 'reason': str(error),
                'blocks': {('REJECTED_PROGRAM',): np.array(density, copy=True)},
                'pulses_executed': 0}
    return execute_word(admitted, density, reference_dimension=reference_dimension)


def occupation_partial_trace(density, ports, keep):
    """Full occupation-factor partial trace; fixed global k may become mixed k.

    This CP restriction is not the static coordinate-site mean R and is not
    asserted to commute with native dynamics across the retained boundary.
    Like execute, this accepts arbitrary exact rational operators, with a
    complex operator represented by separate real and imaginary calls. It
    validates the carrier, not positivity or normalization of a preparation.
    """
    _integer(ports, 1, 4096)
    if type(keep) is not tuple:
        raise ValueError('canonical retained port subset required')
    for port in keep:
        _integer(port, 0, ports-1)
    if tuple(sorted(set(keep))) != keep:
        raise ValueError('canonical retained port subset required')
    if type(density) is not dict:
        raise ValueError('finite sparse source/reference operator required')
    mask = sum(1 << p for p in keep)
    output = defaultdict(F)
    for key, coefficient in density.items():
        if type(key) is not tuple or len(key) != 4:
            raise ValueError('operator key must be source,reference,source,reference')
        for index in key:
            _integer(index)
        i, a, j, b = key
        if max(i, j) >= 1 << ports:
            raise ValueError('occupation state exceeds the finite carrier')
        if type(coefficient) not in (int, F):
            raise ValueError('exact rational coefficients required; apply real and imaginary parts separately')
        if i & ~mask == j & ~mask:
            ii = sum(((i >> p) & 1) << k for k, p in enumerate(keep))
            jj = sum(((j >> p) & 1) << k for k, p in enumerate(keep))
            output[ii, a, jj, b] += F(coefficient)
    return dict(output)
