"""Native coin-prefix consensus and a nonconstant record Gibbs family."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
from itertools import product
import json
from pathlib import Path

from tools import strict_json
from ..native_born_bridge import certificate as born
from ..native_thermodynamic_bridge import certificate as entropy

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PREFIXES = ((0,), (1, 0), (1, 1))
COSTS = (1, 2, 2)
WEIGHTS = (Q(1, 2), Q(1, 4), Q(1, 4))
require = born.require


def integer(value, minimum, maximum, name):
    require(type(value) is int and minimum <= value <= maximum, f'{name}: original integer required')
    return value


def parameters(beta, cap):
    integer(beta, 1, 128, 'replica count')
    integer(cap, 1, 4096, 'proposal cap')


def source_state(value):
    require(type(value) is tuple and all(type(x) is int for x in value)
            and value in born.INDEX, 'complete total-two native configuration required')
    return value


def reference(beta):
    integer(beta, 1, 128, 'replica count')
    masses = tuple(w**beta for w in WEIGHTS)
    z = sum(masses)
    return tuple(a/z for a in masses), z


def kernel(beta, cap=1):
    parameters(beta, cap)
    pi, z = reference(beta)
    fail = (1-z)**cap
    return tuple(tuple((fail if i == j else 0)+(1-fail)*pi[j] for j in range(3)) for i in range(3))


def complete_register_law(initial, beta, cap=1):
    """Public agreement-label/timeout marginal, before any conditioning."""
    integer(initial, 0, 2, 'incoming register')
    parameters(beta, cap)
    pi, z = reference(beta)
    fail = (1-z)**cap
    return {('accept', h): (1-fail)*pi[h] for h in range(3)} | {('timeout', initial): fail}


@dataclass(frozen=True)
class Prefix:
    label: int
    records: tuple
    output: tuple


def prefix(initial, coins, seam=born.native.E):
    """The source is hidden from the stopping policy, which reads only coins."""
    state = source_state(initial)
    born.seam(seam)
    stream = iter(coins)
    records = []
    bits = []
    for _ in range(2):
        try:
            coin = next(stream)
        except StopIteration as exc:
            raise ValueError('native prefix coin stream ended') from exc
        integer(coin, 0, 1, 'native coin')
        state, record = born.endpoint(state, seam, coin)
        records.append(record)
        bits.append(coin)
        if bits == [0] or len(bits) == 2:
            return Prefix(PREFIXES.index(tuple(bits)), tuple(records), state)
    raise ArithmeticError('prefix controller failed to halt')


@dataclass(frozen=True)
class Execution:
    flag: str
    register: int
    batches: tuple
    output: tuple
    attempts: int


def execute(initial, coins, register, *, beta=2, cap=3, seam=born.native.E):
    parameters(beta, cap)
    integer(register, 0, 2, 'incoming register')
    state = source_state(initial)
    stream = iter(coins)
    batches = []
    attempts = 0
    for _ in range(cap):
        batch = []
        for _ in range(beta):
            item = prefix(state, stream, seam)
            state = item.output
            attempts += len(item.records)
            batch.append(item)
        batches.append(tuple(batch))
        if len({item.label for item in batch}) == 1:
            return Execution('accept', batch[0].label, tuple(batches), state, attempts)
    return Execution('timeout', register, tuple(batches), state, attempts)


def source_law():
    """Actual stopped records on all 78 initial source configurations."""
    table = []
    for state in born.BASIS:
        row = [Q()]*3
        for word in PREFIXES:
            item = prefix(state, word)
            require(item.label == PREFIXES.index(word) and len(item.records) == len(word),
                    'native prefix did not implement its record stopping rule')
            row[item.label] += Q(1, 2**len(item.records))
        require(tuple(row) == WEIGHTS, 'source input changed the prefix distribution')
        table.append(row)
    return tuple(map(tuple, table))


def independent_proposal(beta):
    """Exact finite prefix-history enumeration, retaining rejection mass and cost."""
    integer(beta, 1, 6, 'enumerated replica count')
    accepted = [Q()]*3
    reject, cost = Q(), Q()
    for labels in product(range(3), repeat=beta):
        probability = Q(1)
        for h in labels:
            probability *= WEIGHTS[h]
        cost += probability*sum(COSTS[h] for h in labels)
        if len(set(labels)) == 1:
            accepted[labels[0]] += probability
        else:
            reject += probability
    pi, z = reference(beta)
    require(sum(accepted)+reject == 1 and reject == 1-z, 'proposal lost disagreement histories')
    require(tuple(a/z for a in accepted) == pi and cost == Q(3*beta, 2),
            'native prefix consensus failed Gibbs or cost identity')
    return tuple(accepted), reject, cost


def enumerated_controller(initial, register, *, beta=2, cap=2):
    """Independent finite source executor versus the closed register formula.

    Unused suffix labels pad early agreement paths; their probabilities sum
    to one and are not counted as executed attempts or retained records.
    """
    parameters(beta, cap)
    require(beta*cap <= 6, 'at most six prefix labels in the finite path oracle')
    out = {key: Q() for key in complete_register_law(register, beta, cap)}
    cost = Q()
    for labels in product(range(3), repeat=beta*cap):
        probability = Q(1)
        coins = []
        for h in labels:
            probability *= WEIGHTS[h]
            coins.extend(PREFIXES[h])
        result = execute(initial, tuple(coins), register, beta=beta, cap=cap)
        out[result.flag, result.register] += probability
        cost += probability*result.attempts
    require(out == complete_register_law(register, beta, cap),
            'actual bounded comparison executor disagrees with register law')
    require(cost == expected_attempts(beta, cap), 'actual rejected-history resource cost disagrees')
    return out


def expected_attempts(beta, cap=None):
    integer(beta, 1, 128, 'replica count')
    pi, z = reference(beta)
    if cap is None:
        return Q(3*beta, 2)/z
    integer(cap, 1, 4096, 'proposal cap')
    return Q(3*beta, 2)*(1-(1-z)**cap)/z


def energy(p):
    """Coefficient of log(2) in the operational expected prefix cost."""
    p = entropy.law(p, size=3)
    return sum(p[i]*COSTS[i] for i in range(3))


def check_kernel(beta, cap):
    pi, z = reference(beta)
    k = kernel(beta, cap)
    require(all(sum(row) == 1 for row in k), 'register kernel not stochastic')
    require(entropy.push(pi, k) == pi, 'counted reference not stationary')
    require(all(pi[i]*k[i][j] == pi[j]*k[j][i] for i in range(3) for j in range(3)),
            'register detailed balance failed')
    square = entropy.multiply(k, k)
    require((square == k) == (beta == 1), 'finite rejection kernel incorrectly identified as projector')
    p = (Q(1, 2), Q(1, 3), Q(1, 6))
    q = entropy.push(p, k)
    # Exact free-energy decrease: entropy split and data-processing residual.
    lhs = entropy.combine((1, entropy.entropy(q)), (-1, entropy.entropy(p)),
                          (-beta*(energy(q)-energy(p)), ((2, Q(1)),)))
    rhs = entropy.combine((1, entropy.relative_entropy(p, pi)), (-1, entropy.relative_entropy(q, pi)))
    require(lhs == rhs, 'dimensionless Clausius split failed')
    # The reference-preserving fluctuation exponent is a ratio of probabilities.
    fluct = sum((p[i]*k[i][j]*(q[j]*pi[i])/(p[i]*pi[j]) for i in range(3) for j in range(3)), Q())
    require(fluct == 1, 'positive-support register fluctuation identity failed')
    if beta == 1:
        require(lhs == entropy.relative_entropy(p, q), 'exact heat-bath dissipation failed')
    return {'beta': beta, 'cap': cap, 'reference': pi, 'agreement_probability_per_proposal': z,
            'timeout_probability': (1-z)**cap, 'projector': square == k,
            'mean_attempts_with_cap': expected_attempts(beta, cap),
            'mean_attempts_until_agreement': expected_attempts(beta),
            'native_attempts_max': 2*beta*cap,
            'positive_support_fluctuation_average': fluct,
            'excited_mass': pi[1]+pi[2], 'excited_mass_upper': Q(2, 2**beta)}


def encode(value):
    if isinstance(value, Q):return str(value)
    if isinstance(value, dict):return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):return [encode(v) for v in value]
    return value


def report():
    inherited = strict_json.load(ROOT/'model/native_born_bridge/results.json')
    for path, expected in inherited['source_sha256'].items():
        require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected,
                f'upstream native record source binding changed: {path}')
    table = source_law()
    for beta in (1, 2, 3, 4):
        independent_proposal(beta)
    for state in (born.BASIS[0], born.BASIS[-1]):
        for register in range(3):
            enumerated_controller(state, register)
    rows = [check_kernel(beta, cap) for beta in (1, 2, 3, 8) for cap in (1, 3)]
    pi2, _ = reference(2)
    require(pi2 == (Q(2, 3), Q(1, 6), Q(1, 6)) and pi2 != WEIGHTS,
            'nonindependent-copy control lost its distinction')
    require(pi2[0].denominator & (pi2[0].denominator-1) != 0,
            'bounded fair-coin nonrealizability control became dyadic')
    paths = ['model/native_record_thermodynamics/certificate.py', 'model/native_record_thermodynamics/PROOF.md',
             'model/native_record_thermodynamics/test_certificate.py', 'model/native_born_bridge/certificate.py',
             'model/native_born_bridge/results.json', 'model/native_thermodynamic_bridge/certificate.py',
             'model/observer/native.py', 'Lean/Thermodynamics/FourLawAdequacySurface.lean']
    return encode({'schema': 'oph.native_record_thermodynamics.v1', 'all_pass': True,
        'native_prefix': {'source_inputs': len(table), 'prefixes': PREFIXES, 'weights': WEIGHTS,
                          'energy_coefficients_in_log2_units': COSTS,
                          'native_stopped_paths_checked': len(table)*3,
                          'bounded_controller_padded_histories_checked': 2*3*3**4,
                          'conditional_prefix_law_independent_of_every_source_input': True,
                          'independent_carrier_replicas_or_reset_required': False},
        'register_kernels': rows,
        'controls': {'copy_one_label_beta2_accepted_law': WEIGHTS, 'independent_beta2_accepted_law': pi2,
                     'beta2_exact_sampler_uniformly_bounded_fair_coin_tree': False,
                     'finite_beta_greater_than_one_kernel_idempotent': False},
        'scope': {'domain': 'declared three-state observer record register',
                  'energy_defined_from_native_prefix_attempt_cost_before_reference': True,
                  'beta_selected_by_integer_replica_comparison_controller': True,
                  'beta_one_exact_heat_bath_max_native_attempts': 2,
                  'beta_greater_than_one_exact_heat_bath_requires_unbounded_retry': True,
                  'finite_caps_keep_every_rejection_and_timeout': True,
                  'full_source_and_growing_archive_stationary': False,
                  'physical_temperature_energy_work_or_erasure_cost_identified': False,
                  'same_H_as_occupation_or_transfer_reconstruction': False,
                  'cofinal_spacetime_low_temperature_control': False},
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true');parser.add_argument('--check', action='store_true')
    args = parser.parse_args(); result = report(); path = HERE/'results.json'
    if args.write:path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check:require(strict_json.load(path) == result, 'native record thermodynamics receipt changed')
    print(json.dumps({'all_pass': True, 'record_states': 3, 'source_inputs': len(born.BASIS),
                      'physical_temperature_derived': False}))


if __name__ == '__main__':main()
