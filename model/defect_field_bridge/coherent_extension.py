"""Complete CAR repair with a coherent, unflagged nonbinary hold.

This is an explicitly chosen extension of the binary C/E law. The flagged
extension remains available in full_domain.py as an independent control.
The two extensions agree on every old hard-core two-particle branch.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ..common_source import native_routing, record_code
from ..defect_field_interactions import interaction
from ..exchange_sector import statistics_service as source
from ..observer.source_scalar import compiler
from . import full_domain

RAW_SLOTS = full_domain.RAW_SLOTS
SOURCE_BOUND = full_domain.SOURCE_BOUND


def branch_action(mask, law, coin, *, sites=12, flavors=90, seam=(0, 5)):
    """Keep binary columns; coherently append nonbinary identity to raw zero."""
    raw, output, sign, current = full_domain.branch_action(
        mask, law, coin, sites=sites, flavors=flavors, seam=seam)
    if raw == full_domain.COLLISION_RAW:
        raw = 0
    return raw, output, sign, current


def finite_branches(*, sites, flavors, law, seam):
    modes = full_domain._domain(sites, flavors, seam)
    if modes > 12:
        raise ValueError('explicit complete-Fock oracle is bounded to twelve modes')
    rows = {(coin, raw): [] for coin in (0, 1) for raw in range(RAW_SLOTS)}
    for coin in (0, 1):
        for mask in range(1 << modes):
            raw, output, sign, _ = branch_action(mask, law, coin, sites=sites,
                                                flavors=flavors, seam=seam)
            rows[coin, raw].append((mask, output, sign))
    for row in rows.values():
        if len({item[1] for item in row}) != len(row):
            raise ArithmeticError('coherent direct-sum branch lost output injectivity')
    return {key: tuple(value) for key, value in rows.items()}


def channel(density, *, sites, flavors, law, seam):
    rows = finite_branches(sites=sites, flavors=flavors, law=law, seam=seam)
    dimension = 1 << (sites*flavors)
    value = np.asarray(density, dtype=complex)
    if value.shape != (dimension, dimension):
        raise ValueError('complete density matrix on the specified CAR space required')
    output = np.zeros_like(value)
    for row in rows.values():
        inputs = np.array([a for a, _, _ in row], dtype=int)
        outputs = np.array([b for _, b, _ in row], dtype=int)
        signs = np.array([sign for _, _, sign in row], dtype=int)
        output[np.ix_(outputs, outputs)] += value[np.ix_(inputs, inputs)]*np.outer(signs, signs)/2
    return output


@lru_cache(maxsize=4)
def source_isometry(sites, flavors, law):
    """Literal complete source dilation with the same512-address service."""
    if (sites, flavors) not in ((2, 2), (3, 1)):
        raise ValueError('source witnesses are the collision and crossing cases')
    modes, seam = sites*flavors, (0, sites-1)
    dimension = 1 << modes
    assignments, quarters = [], [0]*512
    target = np.zeros((512, dimension), complex)
    for coin in (0, 1):
        for mask in range(dimension):
            raw, output, sign, _ = branch_action(mask, law, coin, sites=sites,
                                                flavors=flavors, seam=seam)
            origin = mask+(coin << modes)
            destination = output+(coin << modes)+(raw << (modes+1))
            assignments.append((origin, destination))
            quarters[origin] = 2 if sign == -1 else 0
            target[destination, mask] = sign/math.sqrt(2)
    permutation = native_routing.complete_permutation(tuple(assignments))
    word = [compiler.pulse('Y', -math.pi/4, [modes])]+source.diagonal_word(tuple(quarters))
    actual = compiler.execute(word, 9, np.eye(512, dtype=complex)[:, :dimension])
    actual = source._permutation(permutation, actual)
    if np.linalg.norm(actual-target) > SOURCE_BOUND:
        raise ArithmeticError('coherent source dilation exceeded its full-column admission')
    counts = source._add(compiler.resources(word), record_code._permutation_resources(permutation))
    return actual, target, source._failure_account(counts)


def even_endpoint_identity(sites, flavors, law, seam):
    """Exact identity-on-subspace certificate, including arbitrary references."""
    modes = full_domain._domain(sites, flavors, seam)
    if modes > 12:
        raise ValueError('bounded explicit identity witness required')
    local_mask = (1 << flavors)-1
    admitted = []
    for mask in range(1 << modes):
        if any(((mask >> (site*flavors)) & local_mask).bit_count() % 2 for site in seam):
            continue
        for coin in (0, 1):
            if branch_action(mask, law, coin, sites=sites, flavors=flavors, seam=seam) != (0, mask, 1, 0):
                raise ArithmeticError('an even-endpoint column no longer shares identity raw zero')
        admitted.append(mask)
    return tuple(admitted)


def source_echo(energy=.102759757, duration=1.):
    gate, _ = interaction.pair_gate(energy, math.pi/2, duration)
    initial = np.zeros(16, complex)
    initial[:4] = gate[:, 0]
    inverse = np.kron(np.eye(4), gate.conj().T)
    result = {}
    for law in ('C', 'E'):
        isometry, _, _ = source_isometry(2, 2, law)
        output = (inverse@(isometry@initial).reshape(32, 16).T).T
        returned = float(np.sum(abs(output[:, 0])**2))
        if abs(returned-1) > 4*SOURCE_BOUND:
            raise ArithmeticError('unflagged source hold failed to preserve the inverse-drive echo')
        result[law] = {'return_probability': returned,
                       'all_record_field_amplitudes': output}
    return result


def certificate():
    audits, executions, arrays = [], [], {}
    for sites, flavors in ((2, 2), (3, 1), (4, 2)):
        for law in ('C', 'E'):
            seam = (0, sites-1)
            rows = finite_branches(sites=sites, flavors=flavors, law=law, seam=seam)
            dimension = 1 << (sites*flavors)
            in_counts, out_counts = [0]*dimension, [0]*dimension
            for row in rows.values():
                for origin, output, sign in row:
                    in_counts[origin] += 1
                    out_counts[output] += 1
                    if origin.bit_count() != output.bit_count() or sign not in (-1, 1):
                        raise ArithmeticError('coherent repair changed number or CAR amplitudes')
            if set(in_counts) != {2} or set(out_counts) != {2}:
                raise ArithmeticError('complete fair-coin channel is not trace preserving and unital')
            support = even_endpoint_identity(sites, flavors, law, seam)
            audits.append({'sites': sites, 'flavors': flavors, 'law': law,
                           'Fock_dimension': dimension, 'all_number_sectors': sites*flavors+1,
                           'all_coin_raw_slots': len(rows), 'zero_slots': sum(not row for row in rows.values()),
                           'even_endpoint_identity_subspace_dimension': len(support)})
            if (sites, flavors) != (4, 2):
                actual, target, resources = source_isometry(sites, flavors, law)
                arrays[f'{sites}_{flavors}_{law}_isometry'] = actual
                executions.append({'sites': sites, 'flavors': flavors, 'law': law,
                                   'complete_source_input_columns': target.shape[1],
                                   'all_source_output_addresses': target.shape[0],
                                   'operator_Frobenius_error_below': SOURCE_BOUND,
                                   'resources': resources})
    echoes = source_echo()
    echo_rows = {}
    for law, row in echoes.items():
        arrays[f'{law}_complete_echo'] = row.pop('all_record_field_amplitudes')
        echo_rows[law] = row
    flagged = interaction.full_repair_echo()
    paths = ['model/defect_field_bridge/coherent_extension.py',
             'model/defect_field_bridge/test_coherent_extension.py',
             'model/defect_field_bridge/COHERENT_EXTENSION.md',
             'model/defect_field_bridge/full_domain.py',
             'model/defect_field_interactions/interaction.py',
             'model/common_source/native_routing.py', 'model/common_source/record_code.py',
             'model/exchange_sector/statistics_service.py', 'model/observer/source_scalar/compiler.py']
    root = Path(__file__).resolve().parents[2]
    result = {'schema': 'oph.coherent_collision_hold_extension.v1',
              'design': 'K_coin,raw0 = K_binary_coin,raw0 direct-sum I_nonbinary/sqrt(2); other raw maps zero on nonbinary endpoints',
              'C_E_record_scope': 'Existing binary labels unchanged; assigning raw0 to a nonbinary hold is an explicit new extension, including for E.',
              'original_hard_core_colored_restriction_exact': True,
              'full_spatial_domain': 'all2^1080 states of the12-port90-flavor CAR space, with arbitrary reference',
              'completeness_proof': 'Binary outputs remain binary; nonbinary identity outputs are disjoint. Each enlarged branch remains injective; each input appears once per fair coin.',
              'identity_subspace_proof': 'Even endpoint occupancies are either both zero or include a nonbinary endpoint. Every admitted column belongs to raw0 with sign+1 for both coins, so the complete channel is identity on all coherences of this subspace.',
              'coherent_hold_removes_all_repair_backaction': False,
              'all_number_sectors_retained': True, 'unital_fair_coin_channel': True,
              'bounded_exact_audits': audits, 'bounded_source_dilations': executions,
              'source_echo': echo_rows,
              'full_spatial_drive_control': {'phase': math.pi/2, 'duration_GeV_inverse': .002,
                  'current': 0, 'endpoint_vacuum_probability': flagged['two_endpoint_vacuum_probability'],
                  'flagged_hold_return_probability': flagged['inverse_drive_return_probability'],
                  'coherent_hold_return_probability_exact': 1,
                  'postmeasurement_Gaussian_assumption': False},
              'full_1080_mode_source_isometry_executed': False,
              'source_scope': 'Complete16- and8-input source dilations; full-domain extension proved algebraically. Four local source isometries and all echo outputs are retained.',
              'physical_particle_identification_complete': False,
              'source_sha256': {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}}
    return result, arrays


def main():
    import argparse
    from .dynamics_packet import compare_arrays, compare_json

    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result, arrays = certificate()
    output = Path(__file__).with_name('coherent_extension_results.json')
    archive = output.with_suffix('.npz')
    if args.write:
        output.write_text(json.dumps(result, indent=2)+'\n')
        np.savez_compressed(archive, **arrays)
    elif args.check:
        compare_json(result, json.loads(output.read_text()))
        with np.load(archive, allow_pickle=False) as saved:
            compare_arrays(arrays, {key: saved[key] for key in saved.files})
    print('complete coherent collision-hold extension verified')


if __name__ == '__main__':
    main()
