"""Complete local repair on every particle-number sector of a spatial CAR space.

The binary endpoint rule is unchanged. An endpoint with more than one
quasiparticle causes a recorded hold. This is an explicit extension, not
an assertion that the original binary rule already specified collisions.
The basis is the static per-site quasiparticle basis, not bare occupation.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ..common_source import native_routing, record_code
from ..exchange_sector import statistics_service as source
from ..observer.source_scalar import compiler


RAW_SLOTS = 16
COLLISION_RAW = 8
SOURCE_BOUND = 2e-7


def _domain(sites, flavors, seam):
    if type(sites) is not int or sites < 2 or type(flavors) is not int or flavors < 1:
        raise ValueError('at least two sites and one internal state required')
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(x) is not int for x in seam)
            or not 0 <= seam[0] < seam[1] < sites):
        raise ValueError('two increasing valid endpoint indices required')
    return sites*flavors


def branch_action(mask, law, coin, *, sites=12, flavors=90, seam=(0, 5)):
    """One exact basis column: (raw record, output mask, CAR sign, current).

    Occupation order is site outer, internal state inner. Nonbinary local
    occupancies are held and recorded; nothing is removed or renormalized.
    """
    modes = _domain(sites, flavors, seam)
    if type(mask) is not int or not 0 <= mask < (1 << modes):
        raise ValueError('an exact CAR occupation mask in the full domain required')
    if law not in ('C', 'E') or type(coin) is not int or coin not in (0, 1):
        raise ValueError('C/E law and an exact binary coin required')
    local_mask = (1 << flavors)-1
    a, b = seam
    left, right = ((mask >> (site*flavors)) & local_mask for site in seam)
    na, nb = left.bit_count(), right.bit_count()
    if na > 1 or nb > 1:
        return COLLISION_RAW, mask, 1, 0
    total = na+nb
    current = na-(total//2+coin*(total % 2))
    raw = {0: 0, 1: 1, -1: 2}[current]
    if law == 'E':
        raw += 4*(total % 2)
    if not current:
        return raw, mask, 1, current
    origin, target, internal = (a, b, left) if current == 1 else (b, a, right)
    flavor = internal.bit_length()-1
    i, j = origin*flavors+flavor, target*flavors+flavor
    low, high = sorted((i, j))
    between = ((1 << high)-1) ^ ((1 << (low+1))-1)
    sign = (-1)**(mask & between).bit_count()
    return raw, mask ^ (1 << i) ^ (1 << j), sign, current


def finite_branches(*, sites, flavors, law, seam):
    """All exact columns on a bounded oracle, with every padded raw slot."""
    modes = _domain(sites, flavors, seam)
    if modes > 12:
        raise ValueError('explicit full-Fock oracle is bounded to twelve modes')
    rows = {(coin, raw): [] for coin in (0, 1) for raw in range(RAW_SLOTS)}
    for coin in (0, 1):
        for mask in range(1 << modes):
            raw, output, sign, _ = branch_action(mask, law, coin, sites=sites,
                                                 flavors=flavors, seam=seam)
            rows[coin, raw].append((mask, output, sign))
    for row in rows.values():
        if len({x[1] for x in row}) != len(row):
            raise ArithmeticError('recorded full-domain repair is not a partial injection')
    return {key: tuple(value) for key, value in rows.items()}


def channel(density, *, sites, flavors, law, seam):
    """Complete small oracle channel, retaining the fair coin normalization."""
    rows = finite_branches(sites=sites, flavors=flavors, law=law, seam=seam)
    dimension = 1 << (sites*flavors)
    value = np.asarray(density, dtype=complex)
    if value.shape != (dimension, dimension):
        raise ValueError('a full-Fock density matrix of the specified size required')
    result = np.zeros_like(value)
    for row in rows.values():
        inputs = np.array([a for a, _, _ in row], dtype=int)
        outputs = np.array([b for _, b, _ in row], dtype=int)
        signs = np.array([s for _, _, s in row], dtype=int)
        result[np.ix_(outputs, outputs)] += value[np.ix_(inputs, inputs)]*np.outer(signs, signs)/2
    return result


@lru_cache(maxsize=4)
def source_isometry(sites, flavors, law):
    """Execute a complete small source dilation, not merely its target map.

    Every field input column is retained. The fair coin and four record
    bits start blank. Remaining addresses are padding of the nine-bit
    source register. The permutation uses the executed native Gray macro.
    """
    seam = (0, sites-1)
    modes = _domain(sites, flavors, seam)
    if (sites, flavors) not in ((2, 2), (3, 1)):
        raise ValueError('source witnesses are the collision and crossing cases')
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
    images = native_routing.complete_permutation(tuple(assignments))
    coin_word = [compiler.pulse('Y', -math.pi/4, [modes])]
    sign_word = source.diagonal_word(tuple(quarters))
    initial = np.eye(512, dtype=complex)[:, :dimension]
    actual = compiler.execute(coin_word+sign_word, 9, initial)
    actual = source._permutation(images, actual)
    if np.linalg.norm(actual-target) > SOURCE_BOUND:
        raise ArithmeticError('full-domain source dilation failed its complete input-column bound')
    counts = source._add(compiler.resources(coin_word+sign_word),
                         record_code._permutation_resources(images))
    return actual, target, source._failure_account(counts)


def certificate():
    audits = []
    for sites, flavors in ((2, 2), (3, 1), (4, 2)):
        for law in ('C', 'E'):
            seam = (0, sites-1)
            rows = finite_branches(sites=sites, flavors=flavors, law=law, seam=seam)
            dimension = 1 << (sites*flavors)
            counts = [0]*dimension
            output_counts = [0]*dimension
            for row in rows.values():
                for origin, output, sign in row:
                    if origin.bit_count() != output.bit_count() or sign not in (-1, 1):
                        raise ArithmeticError('repair changed total number or violated unit amplitude')
                    counts[origin] += 1
                    output_counts[output] += 1
            if set(counts) != {2} or set(output_counts) != {2}:
                raise ArithmeticError('complete effects or output projectors do not sum to identity')
            audits.append({'sites': sites, 'internal_states': flavors, 'law': law,
                           'all_particle_number_sectors': sites*flavors+1,
                           'full_Fock_dimension': dimension,
                           'complete_coin_and_raw_slots': 2*RAW_SLOTS,
                           'zero_slots': sum(not row for row in rows.values()),
                           'collision_columns': sum(len(rows[coin, COLLISION_RAW]) for coin in (0, 1))})
    executions = []
    for sites, flavors in ((2, 2), (3, 1)):
        for law in ('C', 'E'):
            _, target, counts = source_isometry(sites, flavors, law)
            executions.append({'sites': sites, 'internal_states': flavors, 'law': law,
                               'all_input_columns': target.shape[1], 'padded_output_addresses': 512,
                               'initialized_coin_bits': 1, 'initialized_record_bits': 4,
                               'isometry_Frobenius_error_below': SOURCE_BOUND, 'resources': counts})
    # A pair rotation can create a same-site collision from the vacuum.
    # The extra record explicitly destroys their off-diagonal coherence.
    psi = np.zeros(16, complex)
    psi[0] = psi[3] = 1/math.sqrt(2)
    before = np.outer(psi, psi.conj())
    after = channel(before, sites=2, flavors=2, law='C', seam=(0, 1))
    if abs(after[0, 3]) > 1e-15 or abs(np.trace(after@after).real-.5) > 1e-15:
        raise ArithmeticError('collision-record disturbance control changed')
    root = Path(__file__).resolve().parents[2]
    paths = ['model/defect_field_bridge/full_domain.py',
             'model/common_source/native_routing.py', 'model/common_source/record_code.py',
             'model/exchange_sector/statistics_service.py',
             'model/observer/source_scalar/compiler.py']
    return {'schema': 'oph.full_spatial_CAR_repair.v1',
            'design': 'binary endpoint C/E balancing, otherwise a recorded collision hold',
            'basis': 'site outer, static per-site quasiparticle internal state inner',
            'collision_raw_record': COLLISION_RAW, 'raw_slots_per_coin': RAW_SLOTS,
            'general_domain': 'all 2^(sites*internal_states) CAR basis states, with arbitrary reference',
            'general_proof': 'Within each recorded branch the map is injective: holds are identities; a signed hop is inverted by the recorded seam/current. Distinct inputs give orthogonal outputs. Every input occurs once for each fair coin. Thus sum K_h^dagger K_h=I on all sectors.',
            'unital_repair': 'Every output has two preimages across the recorded branches: two coin holds on local counts zero, two or nonbinary; two possible initial endpoints for local count one. Thus sum K_h K_h^dagger=I and the finite normalized identity is a faithful stationary state of each fixed-seam repair channel.',
            'unital_scope': 'This does not assert stationarity under arbitrary record-conditioned field operations or feedback schedules. The tracial reference has trivial modular flow.',
            'original_restriction': 'On the two-particle hard-core subspace the collision flag never occurs and every existing C/E colored branch, CAR sign and weight is unchanged.',
            'history_composition': 'Any complete unitary propagation or particle-changing operation on the same full CAR space can precede or follow this instrument. Every resulting record is retained; no sector is projected out of the dynamics.',
            'collision_disturbance': {'vacuum_pair_superposition_output_purity': .5,
                                    'vacuum_pair_coherence_after_unread_collision_record': 0,
                                    'hold_outcome_is_not_a_no_disturbance_claim': True},
            'exact_bounded_audits': audits, 'source_executions': executions,
            'full_1080_mode_source_execution_claimed': False,
            'source_scope': 'All input columns of the four stated small witnesses. General extension is exact finite algebra and a specified algorithm, not an executed 1080-mode dilation.',
            'floating_execution_is_interval_or_hardware_certificate': False,
            'physical_particle_identification_complete': False,
            'source_sha256': {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}}


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = Path(__file__).with_name('full_domain_results.json')
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    elif args.check:
        if json.loads(path.read_text()) != result:
            raise SystemExit('full-domain repair receipt drift')
        print('full-domain repair receipt OK')
    else:
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
