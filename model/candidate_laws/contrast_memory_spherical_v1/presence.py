"""A reversible pair-presence instrument and its two-reading observer record law."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import carrier
from . import reader, sectors

HERE = Path(__file__).resolve().parent
SOURCE_FILES = ('__init__.py', '../../executable_model/carrier.py', 'common.py', 'law.py', 'reader.py', 'sectors.py',
                'presence.py', 'presence_protocol.json')


def check(condition, message):
    if not condition:
        raise ValueError(message)


def load_json(path):
    def unique(items):
        result = dict(items)
        check(len(result) == len(items), 'duplicate JSON key')
        return result
    def nonfinite(token):
        raise ValueError(f'nonfinite JSON constant: {token}')
    return json.loads(path.read_text(encoding='utf8'), object_pairs_hook=unique,
                      parse_constant=nonfinite)


def pins():
    return {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in SOURCE_FILES}


def tree_layers(leaves):
    """Successive disjoint (left, right, target) reversible OR operations."""
    if type(leaves) is not int or leaves < 1:
        raise ValueError('at least one leaf is required')
    active = list(range(leaves))
    layers = []
    target = leaves
    while len(active) > 1:
        layer = []
        parents = []
        for left, right in zip(active[::2], active[1::2]):
            layer.append((left, right, target))
            parents.append(target)
            target += 1
        if len(active) % 2:
            parents.append(active[-1])
        layers.append(layer)
        active = parents
    return layers


def apply_layer(work, layer):
    for left, right, target in layer:
        work[:, target] ^= work[:, left] | work[:, right]


def status_bits(public):
    """OR flags: any yes, any mismatch, any unavailable, any unread/malformed."""
    return np.select([public == 29, public == 1, public == 2, public == 3],
                     [1, 0, 2, 4], default=8)


def compile_records(flags, controls, keep_seams):
    """Execute the apparatus permutation on fresh workspace for every basis input.

    Each column of flags is the diagonal predicate controlling one native reader.
    The register copied to the archive only depends on the readers' public bits.
    """
    flags = np.asarray(flags)
    check(flags.ndim == 2 and flags.dtype == bool and flags.shape[1] == len(controls),
          'boolean source predicates and one reader per seam required')
    check(type(keep_seams) is bool, 'address switch must be boolean')
    check(all(control in reader.CONTROLS for control in controls), 'unknown apparatus control')
    size, leaves = flags.shape
    layers = tree_layers(leaves)
    tables = {control: [reader.stages(flag, **reader.CONTROLS[control]) for flag in (0, 1)]
              for control in set(controls)}
    domain = np.arange(reader.APP)
    for settings in tables.values():
        for stages in settings:
            check(len(stages) == 4, 'four native reader stages required')
            for stage in stages:
                check(stage.shape == domain.shape and np.issubdtype(stage.dtype, np.integer)
                      and np.all((stage >= 0) & (stage < reader.APP))
                      and np.array_equal(stage[stage], domain), 'reader stage fails involution')
    words = np.zeros((size, leaves), dtype=int)

    def read_stage(stage):
        for edge, control in enumerate(controls):
            off, on = tables[control]
            words[:, edge] = np.where(flags[:, edge], on[stage][words[:, edge]],
                                      off[stage][words[:, edge]])

    for stage in range(4):
        read_stage(stage)
    codes = words >> 6
    work = np.zeros((size, 2*leaves-1), dtype=int)
    work[:, :leaves] ^= status_bits(codes)
    for layer in layers:
        apply_layer(work, layer)
    # XOR-copy into fresh archive registers; both switches allocate the same space.
    aggregate = work[:, -1].copy()
    addresses = codes.copy() if keep_seams else np.zeros_like(codes)
    for layer in reversed(layers):
        apply_layer(work, layer)
    work[:, :leaves] ^= status_bits(words >> 6)
    for stage in reversed(range(4)):
        read_stage(stage)
    check(not np.any(words) and not np.any(work), 'private workspace retained source information')
    return [(int(label), tuple(map(int, address))) for label, address in zip(aggregate, addresses)]


def partition(records):
    """Each group is one private-free Kraus projector onto its basis indices."""
    groups = defaultdict(list)
    for index, label in enumerate(records):
        groups[label].append(index)
    return dict(sorted(groups.items()))


def flags_with_fuel(sector):
    loads = np.asarray([(0,)*sector.car.n, *sector.states])
    edges = np.asarray(sector.car.seams)
    return (loads[:, edges[:, 0]] == 1) & (loads[:, edges[:, 1]] == 1)


def paired_histories(source, prepared, diagonal, first, last, delay, window):
    """Both samplings include the diagonal source phase through the read window."""
    # Each column is an unnormalized first branch, with its full residual source.
    columns = np.zeros((len(prepared), len(first)), dtype=complex)
    for index, indices in enumerate(first.values()):
        columns[indices, index] = prepared[indices]
    columns *= np.exp(-1j*window*diagonal)[:, None]
    columns = expm_multiply(-1j*delay*source, columns)
    # The last diagonal phase commutes through every output projector and
    # cancels from the probabilities. Keep every possible ordered label pair.
    return np.stack([np.sum(abs(columns[indices])**2, axis=0) for indices in last.values()], axis=1)


def controls_for(name, count):
    names = ['nominal'] * count
    if name == 'all common flip':
        names = ['wrong_agreement'] * count
    elif name in ('first-seam mismatch', 'first-seam unavailable'):
        names[0] = name.split()[-1]
    else:
        check(name == 'nominal', 'unknown fault placement')
    return names


def joint_statistics(aggregate, resolved, codes, labels):
    """Hide only the first seam addresses, keeping each final record separate."""
    hidden = np.zeros_like(aggregate)
    for probabilities, (code, _) in zip(resolved, labels, strict=True):
        hidden[codes.index(code)] += probabilities
    margin_error = float(np.max(abs(aggregate.sum(axis=1) - hidden.sum(axis=1))))
    return margin_error, float(abs(aggregate - hidden).sum()/2)


def case(sector, last, flags, kappa, fault, protocol):
    matrix = sparse.csr_matrix(sector.hamiltonian(kappa))
    seed = [0] * sector.car.n
    for port in carrier.EMITTER_SEAM['cell']:
        seed[port] = 1
    pulse, fuel_energy = sectors.emit(matrix, sector.index[tuple(seed)])
    # One common energy origin, including the residual fuel. The emitter
    # returns the state in this same rotating frame.
    source = sparse.block_diag(([0.], matrix - fuel_energy*sparse.eye(matrix.shape[0])),
                               format='csr')
    diagonal = np.r_[0., kappa*sector.C - fuel_energy]
    before, after = protocol['free_actions']
    prepared = expm_multiply(-1j*before*source, pulse)
    controls = controls_for(fault, flags.shape[1])
    first = [partition(compile_records(flags, controls, keep)) for keep in (False, True)]
    joint = [paired_histories(source, prepared, diagonal, groups, last, after, protocol['window'])
             for groups in first]
    aggregate_codes = [code for code, _ in first[0]]
    first_error, separation = joint_statistics(*joint, aggregate_codes, first[1])
    return {'kappa': kappa, 'first_reader': fault, 'residual_fuel': float(abs(pulse[0])**2),
            'aggregate_codes': aggregate_codes,
            'resolved_first_labels': list(first[1]),
            'aggregate_joint': joint[0].tolist(), 'resolved_joint': joint[1].tolist(),
            'first_margin_error': first_error, 'public_joint_tv': separation}


def campaign():
    protocol = load_json(HERE/'presence_protocol.json')
    source_pins = pins()
    check(source_pins == load_json(HERE/'presence_freeze.json'), 'frozen source mismatch')
    sector = sectors.MultisetSector(carrier.cell(), 2)
    flags = flags_with_fuel(sector)
    nominal = ['nominal'] * flags.shape[1]
    last = partition(compile_records(flags, nominal, True))
    cases = [case(sector, last, flags, kappa, fault, protocol)
             for kappa in protocol['laws'] for fault in protocol['faults']
             if kappa == 8 or fault == 'nominal']
    return {'construction': protocol['construction'], 'source_pins': source_pins,
            'final_labels': list(last), 'cases': cases}


def scientific_gates(result):
    protocol = load_json(HERE/'presence_protocol.json')
    expected = [(k, f) for k in protocol['laws'] for f in protocol['faults']
                if k == 8 or f == 'nominal']
    check([(r['kappa'], r['first_reader']) for r in result['cases']] == expected,
          'incomplete campaign or reordered cases')
    for row in result['cases']:
        aggregate = np.asarray(row['aggregate_joint'])
        resolved = np.asarray(row['resolved_joint'])
        codes, labels = row['aggregate_codes'], row['resolved_first_labels']
        for law, rows in ((aggregate, len(codes)), (resolved, len(labels))):
            check(law.shape == (rows, len(result['final_labels']))
                  and np.all(np.isfinite(law)) and np.all(law >= 0)
                  and abs(law.sum() - 1) < 1e-10, 'incomplete joint law')
        margin, separation = joint_statistics(aggregate, resolved, codes, labels)
        check(margin < 1e-10 and row['first_margin_error'] < 1e-10, 'unequal first public margins')
        # Recompute derived fields from the very same saved joint probabilities.
        # The absolute allowance covers the final sums at binary64 precision.
        for key, value in [('first_margin_error', margin), ('public_joint_tv', separation)]:
            reported = row[key]
            check(type(reported) is float and 0 <= reported and abs(reported-value) <= 5e-15,
                  f'inconsistent {key}')
    main = next(row for row in result['cases'] if row['kappa'] == 8 and row['first_reader'] == 'nominal')
    check(main['public_joint_tv'] > 48*protocol['window'],
          'public separation does not exceed the two-experiment pulse budget')


def compare(fresh, stored):
    if isinstance(fresh, (list, tuple)):
        check(type(stored) is list and len(fresh) == len(stored), 'receipt list mismatch')
        for expected, actual in zip(fresh, stored):
            compare(expected, actual)
    elif type(fresh) is dict:
        check(type(stored) is dict and fresh.keys() == stored.keys(), 'receipt keys mismatch')
        for key in fresh:
            compare(fresh[key], stored[key])
    elif type(fresh) is float:
        check(type(stored) is float and math.isfinite(stored) and stored >= 0, 'invalid receipt number')
        # Only the explicitly reported numerical equality residual is allowed an
        # absolute roundoff allowance (handled by verify); joint weights use relative error.
        check(math.isclose(fresh, stored, rel_tol=1e-8, abs_tol=0.), 'receipt weight mismatch')
    else:
        check(type(stored) is type(fresh) and stored == fresh, 'receipt value mismatch')


def verify(fresh, stored):
    # Compare schema and reconstructed labels before doing arithmetic on input.
    # The equality residual can differ by machine roundoff across platforms.
    reference = dict(fresh, cases=[dict(row, first_margin_error=0.) for row in fresh['cases']])
    candidate = dict(stored, cases=[dict(row, first_margin_error=0.) for row in stored['cases']])
    compare(reference, candidate)
    scientific_gates(stored)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--record', type=Path)
    group.add_argument('--verify', type=Path, default=HERE/'presence_results.json')
    args = parser.parse_args(argv)
    result = campaign()
    if args.record:
        with args.record.open('x', encoding='utf8', newline='\n') as stream:
            json.dump(result, stream, separators=(',', ':'), allow_nan=False)
            stream.write('\n')
    else:
        verify(result, load_json(args.verify))
    scientific_gates(result)
    print('Verified presence instrument and six complete two-reading laws.')


if __name__ == '__main__':
    main()
