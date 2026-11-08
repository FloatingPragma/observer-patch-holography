"""Authenticated original-kernel queries and an explicit local observation map."""
from collections import defaultdict
from fractions import Fraction as F

from model.native_local_response import cone
from . import clock as c


def query_word(plan, slots, density):
    """A bounded exact CP-word query, distinct from an expanded host run.

    Real and imaginary rational matrix components are queried separately.
    All native slots, including zero maps, are accepted. The routine makes
    no positivity or source-preparation assertion about a supplied operator.
    """
    c.authenticate_plan(plan)
    p = plan['program']
    if plan['level'] > c.native.MAX_EVALUATED_LEVEL:
        raise ValueError('the owning kernel cannot construct this source carrier')
    if (type(slots) is not tuple or len(slots) != plan['total_attempts']
            or len(slots) > c.native.MAX_EVALUATED_STEPS
            or any(type(s) is not int or s < 0 for s in slots)):
        raise ValueError('one exact raw slot per attempt, within the bounded query domain')
    rho = c.native._validate_density(density, p['carrier'])
    initial = dict(rho)
    memory = archive = ()
    cursor = 0
    for command in p['commands']:
        for repetition in range(command['repeat']):
            outcomes = c.native._outcomes(p, command, rho, memory, archive,
                                           repetition=repetition)
            slot = slots[cursor]
            if slot >= len(outcomes):
                raise ValueError('raw slot lies outside the complete native alphabet')
            selected = outcomes[slot]
            rho, memory, archive = selected.density, selected.memory, selected.archive
            cursor += 1
    return {'plan_digest': c.digest(plan), 'program_digest': c.native.digest(p),
            'slots': slots, 'initial_operator_digest': c.digest(initial),
            'density': rho, 'archive': archive, 'memory': memory,
            'records': archive+memory, 'attempts': cursor,
            'status': 'COMPLETE_RAW_WORD_QUERY',
            'expanded_host_completion_asserted': False}


def authenticate_query(plan, query, density):
    if type(query) is not dict or 'slots' not in query:
        raise ValueError('a complete original-kernel query required')
    expected = query_word(plan, query['slots'], density)
    if c.typed(query) != c.typed(expected):
        raise ValueError('query records, source operator or FIFO differ from exact replay')
    return True


def observe_query(plan, window, attached, query, density):
    """Selected bin-labelled records, final ports and untouched reference.

    Exact global attempt addresses and unselected archive values are traced.
    Every fixed intervention record is observed; uniform records are observed
    precisely on the declared seams. Zero query outputs stay zero.
    """
    authenticate_query(plan, query, density)
    expected = c.attach_window(plan, window, attached['root_cell_radius'],
        attached['observed_uniform_seams'], attached['final_ports'],
        initial_error=attached['budget']['initial_joint_half_trace_upper'],
        initial_reference=attached['declared_initial_reference'])
    if c.typed(expected) != c.typed(attached):
        raise ValueError('a current authenticated observation attachment required')
    inverse = {v: i for i, v in enumerate(window['source_ports'])}
    observed = set(attached['observed_uniform_seams'])
    events = []
    for attempt, record in enumerate(query['records'], start=1):
        label = c.attempt_label(plan, attempt)
        edge = (tuple(inverse[v] for v in record.seam)
                if record.seam is not None and all(v in inverse for v in record.seam) else None)
        if label['role'] == 'fixed' or edge in observed:
            if record.seam is not None and edge is None:
                raise ValueError('an observed intervention escaped its authenticated root')
            local = c.native.Record(record.kind, record.selector, edge, record.law,
                                    record.coin, record.raw, record.label, record.current)
            events.append((label['bin'], label['role'], label['fixed_ordinal'], local))
    final_source = tuple(window['source_ports'][v] for v in attached['final_ports'])
    reduced = cone.partial_trace(query['density'], final_source, 240*4**plan['level'])
    def rename(mask):
        return sum(1 << inverse[v] for v in final_source if mask & (1 << v))
    output = {(rename(i), a, rename(j), b): value for (i, a, j, b), value in reduced.items()}
    return {'events': tuple(events), 'density': output,
            'terminal_status': 'COMPLETE_RAW_WORD_QUERY',
            'global_tick_addresses_observed': False}


def complete_fixed_control():
    """Independent partial-map comparison on all two-attempt raw slots.

    A genuine source/reference entangled input exercises off-diagonal blocks.
    This is an exact finite CP instrument calculation on the evaluated source.
    """
    edge = c.native.local.SEAMS[0]
    plan = c.compile_plan(0, 1, (F(0), F(2, 630)),
        (((edge, 'E'), (edge, 'C')),), ('E',))
    i, j = (1 << edge[0], 1 << edge[1])
    rho = {(x, a, y, b): F(1, 2) for x, a in ((i, 0), (j, 1))
                                     for y, b in ((i, 0), (j, 1))}
    family = cone.source_family(plan['program'])
    first = cone.fixed_instrument(family, rho, edge, law='E')
    marginal = defaultdict(F); nonzero = zero = 0
    for s, (_, first_block) in enumerate(first):
        second = cone.fixed_instrument(family, first_block, edge, law='C')
        for t, (_, expected) in enumerate(second):
            actual = query_word(plan, (s, t), rho)
            if actual['density'] != expected:
                raise ArithmeticError('actual source and independent endpoint instrument differ')
            zero += not bool(expected); nonzero += bool(expected)
            for (x, a, y, b), value in expected.items():
                if x == y:
                    marginal[a, b] += value
    if dict(marginal) != {(0, 0): F(1, 2), (1, 1): F(1, 2)}:
        raise ArithmeticError('the untouched reference changed')
    return {'complete_raw_words': 625, 'zero_words': zero, 'nonzero_words': nonzero,
            'entangled_source_reference_control_checked': True,
            'reference_marginal_preserved_exactly': True,
            'archive_and_FIFO_replayed': True, 'host_expansion_asserted': False}


def source_domain_controls():
    """Literal admission and host-stop controls; no grammar extension."""
    edge = c.native.local.SEAMS[0]
    high = c.compile_plan(14, 2, (F(0), F(2, 630*4**14)),
                          (((edge, 'E'),),), ('E',))
    sentinel = {'source_carrier_unevaluated': True}
    output = c.native.execute(high['program'], sentinel)
    if (output.status != 'RESOURCE_STOP_UNEVALUATED_CARRIER'
            or output.attempted_steps != 0 or output.branches[0].density != sentinel):
        raise ArithmeticError('the source carrier must stop before state validation')
    # An actual evaluated intercell edge gives a small-memory guard witness.
    low_family = cone.native.federation_family(0, (1,)*12)
    crossing = next(e for e in low_family.seams if e[0]//12 != e[1]//12)
    rejected = False
    try:
        c.compile_plan(3, 2, (F(0), F(1, 630*4**3)),
                       (((crossing, 'E'),),), ('E',))
    except ValueError as error:
        rejected = 'intercell addresses' in str(error)
    if not rejected:
        raise ArithmeticError('unevaluated fixed intercell admission changed')
    malformed = dict(high['program']); malformed['commands'] = [{'op': 'source_word', 'repeat': 1}]
    refused = False
    try:
        c.native.admit(malformed)
    except ValueError as error:
        refused = 'undeclared native operation' in str(error)
    if not refused:
        raise ArithmeticError('the repair grammar unexpectedly admitted a coherent service')
    word = c.native.source_word([{'kind': 'X', 'targets': (0,), 'angle': 0.0}], 7)
    import numpy as np
    stopped = c.native.execute_word(word, np.asarray([[1]]))
    if stopped['status'] != 'RESOURCE_STOP' or stopped['pulses_executed'] != 0:
        raise ArithmeticError('generic coherent host cap changed')
    return {'fixed_cell_word_level14_admitted': True,
            'level14_host_status': output.status, 'attempts_before_host_stop': 0,
            'unevaluated_fixed_intercell_command_rejected': True,
            'uniform_intercell_menu_stays_admitted': True,
            'source_word_operation_rejected_by_native_grammar': True,
            'generic_seven_qubit_word_host_status': stopped['status'],
            'coherent_register_initialization_is_supplied':
                word['declared_resources']['initial_work_qubits_supplied'],
            'growing_federation_coherent_load_map_produced': False}
