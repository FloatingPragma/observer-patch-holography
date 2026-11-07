"""Exact FIFO resegmentation for the executable native observer.

This concerns admitted record capacities, with the whole archive retained.
It neither changes the carrier nor identifies FIFO slots with global capacity.
"""
from copy import deepcopy
from fractions import Fraction

from ..native_program import interpreter as source


SOURCE_FILES = (
    'model/particle_identification/native_capacity.py',
    'model/particle_identification/test_native_capacity.py',
    'model/particle_identification/NATIVE_CAPACITY.md',
    'model/native_program/interpreter.py',
    'model/native_program/test_interpreter.py',
)


def with_capacity(program, capacity):
    """Change only an admitted export FIFO size and its counted resources."""
    value = source.admit(program)
    if value['memory']['overflow'] != 'export_oldest':
        raise ValueError('the equivalence requires complete archive export')
    if type(capacity) is not int or not 1 <= capacity <= 8:
        raise ValueError('an admitted integer FIFO capacity from 1 through 8 is required')
    value['memory']['capacity'] = capacity
    value['declared_resources'] = source._native_cost(value)
    return source.admit(value)


def complete_output(history):
    """Forget storage placement, preserving every record and quantum block.

    Program hashes and resource accounts differ between implementations.
    They are not identified with experimental outputs by this map.
    """
    return (history.status, history.attempted_steps, tuple(
        (branch.record, branch.archive + branch.memory,
         branch.halted, branch.density)
        for branch in history.branches))


def resegment(branch, capacity):
    """Reversibly repartition a complete chronological tape."""
    if type(capacity) is not int or not 1 <= capacity <= 8:
        raise ValueError('an admitted integer FIFO capacity from 1 through 8 is required')
    tape = branch.archive + branch.memory
    return source.Outcome(branch.record, deepcopy(branch.density),
                          tape[-capacity:], tape[:-capacity], branch.halted)


def fixture():
    commands = [
        {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [0, 1]},
         'law': 'C', 'repeat': 1},
        {'op': 'hold', 'repeat': 2},
        {'op': 'repair', 'selector': {'kind': 'last_current',
                                    'positive': [0, 5], 'otherwise': [2, 3]},
         'law': 'E', 'repeat': 1},
    ]
    basis = source.local.basis(True)
    i, j = (basis.index(tuple(int(k in pair) for k in range(12)))
            for pair in ((0, 4), (1, 4)))
    # A reference-entangled input tests retained off-diagonal information.
    density = {(x, a, y, b): Fraction(1, 2)
               for x, a in ((i, 0), (j, 1))
               for y, b in ((i, 0), (j, 1))}
    return source.program(commands=commands), density


def certificate():
    program, density = fixture()
    reference = source.execute(program, density)
    output = complete_output(reference)
    for capacity in range(1, 9):
        actual = source.execute(with_capacity(program, capacity), density)
        if complete_output(actual) != output:
            raise ArithmeticError('complete native output depends on export FIFO capacity')
        if tuple(resegment(branch, capacity) for branch in reference.branches) != actual.branches:
            raise ArithmeticError('the full branch maps fail tape resegmentation')
    halt = source.program(commands=[{'op': 'hold', 'repeat': 2}], overflow='halt')
    small = source.execute(halt, density)
    large = source.execute(source.program(commands=[{'op': 'hold', 'repeat': 2}],
                                          capacity=2, overflow='halt'), density)
    if (small.status, large.status) != ('FIFO_TIMEOUT', 'COMPLETE'):
        raise ArithmeticError('the overflow countercontrol failed')
    return {
        'schema': 'oph.native_capacity.resegmentation.v1',
        'admitted_FIFO_capacities': list(range(1, 9)),
        'fixture_attempts': reference.attempted_steps,
        'fixture_complete_raw_blocks': len(reference.branches),
        'fixture_total_probability': str(sum(source.trace(b.density) for b in reference.branches)),
        'complete_record_and_quantum_output_equal': True,
        'theorem_error_after_reversible_record_resegmentation': '0',
        'theorem_domain': 'same admitted native commands, carrier, parameters, input and evaluator limits; export_oldest; complete archive',
        'general_operator_identity_proved_in': 'model/particle_identification/NATIVE_CAPACITY.md',
        'fixtures_are_not_the_general_proof': True,
        'halt_countercontrol': [small.status, large.status],
        'resource_costs_identified': False,
        'specification_reads_of_capacity_identified': False,
        'physical_global_capacity_law_derived': False,
        'source_selected_field_action_derived': False,
    }
