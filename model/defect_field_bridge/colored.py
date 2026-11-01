"""Actual two-fermion exchange on the native hard-core port configurations.

Each occupied port carries its own internal V, rather than assigning one
collective V label to a pair. The field's same-port two-occupancy complement
is counted explicitly and is outside this restricted repair law.
"""
from __future__ import annotations

from math import comb
from pathlib import Path
import hashlib
import json
import sympy as s

from ..defect_protection import protection
from ..exchange_sector import statistics as native


def dimensions(dimension):
    if type(dimension) is not int or dimension < 1:
        raise ValueError('positive exact per-particle internal dimension required')
    return {'internal_per_particle': dimension, 'ports': 12,
            'one_particle_modes': 12*dimension,
            'native_hard_core_two_particle_dimension': 66*dimension**2,
            'full_CAR_two_particle_dimension': comb(12*dimension, 2),
            'same_port_collision_complement_dimension': 12*comb(dimension, 2)}


def address(pair_index, a, b, dimension):
    return (pair_index*dimension+a)*dimension+b


def ordered_cover_embedding(dimension):
    """W=sqrt(2)E=[I;-SWAP] for the actual native orientation-bit lift."""
    dimensions(dimension)
    if dimension > 4:
        raise ValueError('explicit internal matrix witness is bounded to d<=4')
    n = dimension**2
    exchange = s.ImmutableSparseMatrix(n, n, {(b*dimension+a, a*dimension+b): -1
                                              for a in range(dimension) for b in range(dimension)})
    W = s.ImmutableSparseMatrix.vstack(s.ImmutableSparseMatrix(s.eye(n)), exchange)
    flip = s.kronecker_product(s.Matrix([[0, 1], [1, 0]]), s.eye(n))
    if exchange*exchange != s.eye(n) or W.H*W != 2*s.eye(n) or flip*W != W*exchange:
        raise ArithmeticError('actual ordered-cover/colored-field code failed')
    return W, exchange


def colored_pairs(row, dimension):
    """Integer map of sqrt(2)K: the internal flip is -SWAP, not parity."""
    dimensions(dimension)
    return tuple((address(i, a, b, dimension),
                  address(j, b, a, dimension) if sign == -1 else address(j, a, b, dimension),
                  sign)
                 for i, j, sign in zip(row.inputs, row.outputs, row.signs, strict=True)
                 for a in range(dimension) for b in range(dimension))


def _annihilate(mask, mode):
    if not mask >> mode & 1:
        return None
    return mask ^ (1 << mode), (-1)**(mask & ((1 << mode)-1)).bit_count()


def _create(mask, mode):
    if mask >> mode & 1:
        return None
    return mask | (1 << mode), (-1)**(mask & ((1 << mode)-1)).bit_count()


def car_oracle(row, dimension):
    """Independent bitstring c_target^dagger c_origin on all admissible states.

    Native row.inputs/outputs/signs are not consulted. Endpoint balancing
    determines the recorded branch and the fermionic CAR determines sign.
    """
    dimensions(dimension)
    if dimension > 3:
        raise ValueError('explicit CAR oracle is bounded to d<=3; the formula is general')
    left, right = row.seam
    output = []
    for index, pair in enumerate(native.PAIRS):
        n = int(left in pair)+int(right in pair)
        current = int(left in pair)-(n//2+row.coin*(n % 2))
        raw = {0: 0, 1: 1, -1: 2}[current] + (4*(n % 2) if row.law == 'E' else 0)
        if raw != row.raw:
            continue
        for a in range(dimension):
            for b in range(dimension):
                modes = (pair[0]*dimension+a, pair[1]*dimension+b)
                mask, sign = sum(1 << mode for mode in modes), 1
                if current:
                    origin, destination = (left, right) if current == 1 else (right, left)
                    internal = a if pair[0] == origin else b
                    first = _annihilate(mask, origin*dimension+internal)
                    second = _create(first[0], destination*dimension+internal)
                    mask, sign = second[0], first[1]*second[1]
                result_modes = [mode for mode in range(12*dimension) if mask >> mode & 1]
                first, second = [divmod(mode, dimension) for mode in result_modes]
                if first[0] == second[0]:
                    raise ArithmeticError('native hard-core transport entered the collision complement')
                target = native.INDEX[(first[0], second[0])]
                output.append((address(index, a, b, dimension), address(target, first[1], second[1], dimension), sign))
    return tuple(output)


def complete_native_certificate(dimension=3):
    """Exhaustive full-branch CAR witness at bounded d; symbolic lift all d."""
    size = dimensions(dimension)['native_hard_core_two_particle_dimension']
    ordered_cover_embedding(dimension)
    results = {}
    for law in ('C', 'E'):
        rows = protection.branches(law)
        by_seam = {}
        columns = 0
        for row in rows:
            actual, wanted = colored_pairs(row, dimension), car_oracle(row, dimension)
            if actual != wanted:
                raise ArithmeticError('colored repair does not equal the independently applied CAR hopping')
            if len({i for i, _, _ in actual}) != len(actual) or len({j for _, j, _ in actual}) != len(actual):
                raise ArithmeticError('colored native branch is not a partial isometry')
            counts = by_seam.setdefault(row.seam, [0]*size)
            for i, _, sign in actual:
                if sign not in (-1, 1):
                    raise ArithmeticError('nonunitary internal swap coefficient')
                counts[i] += 1
            columns += len(actual)
        if len(by_seam) != 30 or any(set(counts) != {2} for counts in by_seam.values()):
            raise ArithmeticError('complete native colored instrument lost probability')
        loop = protection.exchange_loop(law)
        maps = {(row.law, row.seam, row.coin, row.raw): dict((i, (j, sign)) for i, j, sign in colored_pairs(row, dimension)) for row in rows}
        pair_index = native.INDEX[tuple(loop['initial_pair'])]
        for a in range(dimension):
            for b in range(dimension):
                at, sign = address(pair_index, a, b, dimension), 1
                for label in loop['history']:
                    at, step = maps[(label[0], tuple(label[1]), label[2], label[3])][at]
                    sign *= step
                if at != address(pair_index, b, a, dimension) or sign != -1:
                    raise ArithmeticError('actual odd exchange loop is not minus internal SWAP')
        results[law] = {'all_seams': 30, 'raw_slots': len(rows),
                        'abort_slots': sum(row.abort for row in rows),
                        'zero_slots': sum(not row.inputs for row in rows),
                        'independently_checked_active_columns': columns,
                        'all_column_effect_sums': '2 before unchanged Kraus factor1/sqrt(2)',
                        'closed_odd_exchange_internal_action': '-SWAP'}
    return {'oracle_dimension': dimension, 'native': results,
            'target_dimensions': dimensions(90),
            'general_law': 'K_h=(1/sqrt(2)) sum_x |f_h(x)><x| tensor (-SWAP)^(p_h(x))',
            'CAR_embedding': '|{x<y};a,b> -> c^dagger_(x,a) c^dagger_(y,b)|vacuum>',
            'general_sign_proof': 'A hop across the occupied spectator crosses one creation operator, producing minus and reversing the two internal slots; a noncrossing hop produces neither.',
            'gauge_compatibility': 'SWAP (rho(g) tensor rho(g))=(rho(g) tensor rho(g)) SWAP for every common/global single-particle internal action.',
            'identical_internal_states_exchange': '-1',
            'ordered_cover_code': {
                'isometry': 'E=[I;-SWAP]/sqrt(2)', 'identity': 'X_orientation E = E (-SWAP)',
                'native_branch_join': 'The actual ordered-cover K_h tensor I compressed by E is exactly K_h^colored, on every source column.',
                'natural_added_fiber_dimension': 'd^2',
                'negative_exchange_multiplicity': 'dim Sym^2(V)=d(d+1)/2',
                'positive_exchange_multiplicity': 'dim Lambda^2(V)=d(d-1)/2',
                'target_negative_dimension': 4095, 'target_positive_dimension': 4005,
                'full_orientation_plus_fiber_at_d90': 16200,
                'interpretation': 'Both native characters occur for identical fermions with internal states. Fixed chi=- is the symmetric-internal channel, not the whole two-particle module.'},
            'global_Fock_parity_on_this_N2_domain': '+1; distinct from the true exchange operation',
            'domain': 'Hard-core per-port two-particle sector only. Same-port collision states are counted, not silently deleted from a claimed full Fock law.',
            'not_claimed': ['law on collision states or other particle-number sectors',
                            'site-local gauge covariance without the required link transports',
                            'interacting measured particle spectrum', 'spacetime continuum or Lorentz spin-statistics theorem']}


def certificate():
    result = complete_native_certificate(3)
    root = Path(__file__).resolve().parents[2]
    names = ['model/defect_field_bridge/colored.py', 'model/defect_protection/protection.py',
             'model/exchange_sector/statistics.py']
    result['schema'] = 'oph.defect_field_bridge.colored_exchange.v1'
    result['source_sha256'] = {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in names}
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = Path(__file__).with_name('colored_results.json')
    if args.check:
        if result != json.loads(path.read_text()):
            raise SystemExit('colored exchange receipt drift')
        print('colored exchange receipt OK')
    elif args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(f'wrote {path}')
    else:
        print(json.dumps(result, indent=2))
