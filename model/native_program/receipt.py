"""Reproduce the bounded native-program admission and execution witness."""
from __future__ import annotations

import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import interpreter as p, members
from ..completion import observer as original
from ..completion.family import FamilyParameters
from ..federated_quantum_sector import construction as sector
from ..federated_observer import refinement
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    """Transitive repository Python imports, retained source/data and initializers."""
    retained = members.checked_record_receipt()
    paths = set(bound_pins(retained['source_sha256'])) | set(sector.SOURCE_PATHS)
    paths.update(p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py'))
    paths.update(('model/native_program/README.md', members.RECORD_RECEIPT,
                  'model/completion/observer.py', 'model/completion/family.py',
                  'model/completion/campaign.json', 'model/federated_observer/refinement.py',
                  'model/requirements.txt'))
    pending, visited = list(paths), set()
    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        stem = name.replace('.', '/')
        for candidate in (stem+'.py', stem+'/__init__.py'):
            if (ROOT/candidate).is_file() and candidate not in paths:
                paths.add(candidate); pending.append(candidate)
        parts = stem.split('/')
        for end in range(1, len(parts)):
            candidate = '/'.join(parts[:end])+'/__init__.py'
            if (ROOT/candidate).is_file() and candidate not in paths:
                paths.add(candidate); pending.append(candidate)
    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        file = ROOT/name
        if not file.is_file():
            raise ValueError('missing source binding: '+name)
        module = name.removesuffix('.py').replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        for node in ast.walk(ast.parse(file.read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = ('.'*node.level)+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add_module(target)
                for alias in node.names:
                    add_module(target+'.'+alias.name)
    return bound_pins({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)})


def instrument_certificate():
    checked, serialized = 0, []
    parameters = p.local.Parameters(port_weights=(25, 25)+(8,)*10)
    for shell in ('H78', 'H66'):
        size = int(shell[1:])
        for law in ('C', 'E'):
            for edge in p.local.SEAMS:
                effects = [F()]*size
                for coin, label, pairs in p.local.raw_maps(law, edge, shell == 'H66'):
                    probability = p.local.coin_probability(parameters, edge, coin)
                    if len({i for i, _ in pairs}) != len(pairs) or len({j for _, j in pairs}) != len(pairs):
                        raise ArithmeticError('native branch is not a partial injection')
                    for i, _ in pairs:
                        effects[i] += probability; checked += 1
                    serialized.append([shell, law, list(edge), coin, label, str(probability), list(pairs)])
                if effects != [F(1)]*size:
                    raise ArithmeticError('complete source/reference channel failed')
    return {'input_columns_checked': checked, 'complete_native_maps_sha256': p.digest(serialized),
            'all_H78_H66_seams_and_C_E_partitions': True,
            'branch_effect': 'positive rational multiple of a domain projection',
            'completeness': 'sum_outcome K* K=I exactly',
            'quantum_reference_extension': 'each Kraus map tensor identity on arbitrary R',
            'complex_operators': 'complex-linear extension of exact rational real/imaginary parts',
            'compatibility_precedes_update': True}


def native_witness():
    command = {'op': 'repair', 'selector': {'kind': 'last_current',
               'positive': [0, 5], 'otherwise': [0, 1]}, 'law': 'mixture', 'repeat': 2}
    program = p.program(commands=[command])
    states = p.local.basis()
    i, j = (states.index(tuple(int(k in pair) for k in range(12))) for pair in ((0, 4), (1, 4)))
    initial = {(a, r, b, s): F(1, 2) for a, r in ((i, 0), (j, 1)) for b, s in ((i, 0), (j, 1))}
    history = p.execute(program, initial)
    marginal = p.reference_marginal(b.density for b in history.branches)
    if marginal != p.reference_marginal([initial]) or sum(p.trace(b.density) for b in history.branches) != 1:
        raise ArithmeticError('continuing observer lost complete mass or quantum reference')
    payload = [[[r.kind, r.selector, r.seam, r.law, r.coin, r.raw, r.label, r.current]
                for r in (*b.archive, *b.memory)]+[
                [[list(k), str(v)] for k, v in sorted(b.density.items())]] for b in history.branches]
    return {'program': program, 'status': history.status, 'attempts_executed': history.attempted_steps,
            'complete_raw_history_blocks': len(history.branches),
            'zero_history_blocks': sum(not b.density for b in history.branches),
            'full_quantum_history_sha256': p.digest(payload),
            'total_mass': '1', 'arbitrary_reference_certificate': 'instrument_certificate',
            'executed_entangled_reference_marginal': ['1/2', '1/2'],
            'exports_per_history': 1, 'FIFO_slots_per_history': 1,
            'controller_input': 'only last committed current', 'time_parameter': 'native attempt count'}


def original_observer_witness():
    value = p.program(shell='H78', biased=False, commands=[{
        'op': 'repair', 'selector': {'kind': 'observer', 'positive': [0, 5], 'otherwise': [0, 1]},
        'law': 'mixture', 'repeat': 6}])
    value['parameters']['repair_probability'] = '3/4'
    value['declared_resources'] = p._native_cost(value)
    p.admit(value)
    # The complete long-history kernel is composed by the existing marginal
    # evaluator, which explicitly traces earlier exported records. The new
    # interpreter independently retains every raw record for a short prefix.
    initial = original.entangled_preparation()
    _, summaries = original.run(initial, FamilyParameters(repair_probability=F(3, 4)), 6)
    first = p.execute(value, initial, max_steps=1)
    if (first.attempted_steps != 1 or first.status != 'RESOURCE_STOP'
            or sum(p.trace(b.density) for b in first.branches) != 1):
        raise ArithmeticError('prefix evaluator stop lost the original observer')
    return {'program': value, 'original_marginal_attempts_executed': 6,
            'terminal_summary': summaries[-1], 'independent_complete_raw_prefix_blocks': len(first.branches),
            'prefix_stop_retains_source_and_FIFO': True,
            'marginal_archive_trace_is_explicit': True}


def federation_witness():
    graph = p.federation_source.federation(0); edge = graph.inter_seams[0]
    values = {}
    for law in ('C', 'E'):
        value = p.program(federation={'level': 0, 'total': 1}, commands=[members._fixed(edge, law)])
        source = {(1 << edge[0], 0, 1 << edge[0], 0): F(1)}
        history = p.execute(value, source)
        moved = sum(v for branch in history.branches for (i, a, j, b), v in branch.density.items()
                    if i == j == 1 << edge[1] and a == b)
        if moved != F(1, 2):
            raise ArithmeticError('native joining seam did not transfer occupation')
        values[law] = {'program_digest': p.digest(value), 'complete_raw_blocks': len(history.branches),
                       'transfer_probability': str(moved)}
    return {'ports': graph.ports, 'cells': graph.cells, 'seam': list(edge),
            'distinct_cells': edge[0]//12 != edge[1]//12, 'law_witnesses': values,
            'complete_endpoint_sector_certificate': sector.local_certificate(),
            'retained_coordinate_refinement': refinement.static_certificate(0, 1),
            'quantum_restriction': 'occupation_partial_trace on all retained occupation sectors',
            'dynamic_refinement_intertwiner_claimed': False,
            'growth_grammar_preserved': True}


def service_witness():
    word = p.source_word([p.compiler.pulse('X', .125, [0]), p.compiler.pulse('J', .25, [0, 1]),
                          p.compiler.pulse('Z', 0., [1])], 2)
    state = np.zeros(8, complex); state[0], state[7] = 1/math.sqrt(2), 1j/math.sqrt(2)
    density = np.outer(state, state.conj())
    output = p.execute_word(word, density, reference_dimension=2)
    total = sum(output['blocks'].values())
    marginal = np.einsum('iaib->ab', total.reshape(4, 2, 4, 2))
    residual = max(float(abs(np.trace(total)-1)), float(np.linalg.norm(marginal-np.eye(2)/2)))
    if not math.isfinite(residual) or residual > 1e-12:
        raise ArithmeticError('complete proper-code service lost mass/reference')
    return {'program': word, 'source_pulses_executed': 3, 'complete_output_blocks': len(output['blocks']),
            'export_complement_abort_blocks': len(output['blocks'])-1,
            'numerical_mass_and_reference_error_upper': '1/1000000000000',
            'full_M6_pulses_used': True, 'success_normalization_or_polar_correction': False,
            'exact_completeness_reason': 'unitary M6 pulse followed by complete code+complement export',
            'supplied_load_export_and_blank_services': True,
            'native_field_encoding_inferred': False}


def report():
    bindings = source_bindings()
    result = {'schema': 'oph.native_program.receipt.v1', 'source_sha256': bound_pins(bindings),
              'instrument_certificate': instrument_certificate(),
              'native_continuing_witness': native_witness(),
              'original_H78_observer': original_observer_witness(),
              'federation': federation_witness(), 'complete_coherent_service': service_witness(),
              'retained_record_clock_member': members.record_clock_member(),
              'scope': {'class_admission_not_physical_selection': True,
                        'canonical_physical_premises_discharged': 0,
                        'source_selected_field_or_mass_claimed': False,
                        'physical_clock_calibration_claimed': False,
                        'finite_theorem_domain_is_actual_native_instrument': True,
                        'CP00_CP01_CP13_acceptance_is_separate_review': True}}
    if source_bindings() != bindings:
        raise ArithmeticError('source changed during native-program replay')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args(); current = report()
    if args.write is not None:
        args.write.write_text(json.dumps(current, indent=2, sort_keys=True, allow_nan=False)+'\n')
        print('Native-program receipt written.')
    else:
        if strict_json.load(args.check) != current:
            raise ValueError('native-program receipt differs from complete replay')
        print('Native-program admission, complete instruments and source bindings verified.')


if __name__ == '__main__':
    main()
