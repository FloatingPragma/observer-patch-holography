"""Current-only census composition for native closed-history Gibbs records.

The finite outcome law is derived before conditioning. The transfer Hamiltonian
comes from the native incident kernel; the loop length then selects internal
inverse temperature. No physical temperature or coherent state preparation is
inferred from this reconstructed history ensemble.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..native_clock import configuration as census
from ..observer_time_fields import core, thermal
from ..record_laws import burn_in

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASIS = census.native.states(binary=True)
N = len(BASIS)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, lower, upper, name):
    require(type(value) is int and lower <= value <= upper, f'{name}: exact integer required')
    return value


def rational(value):
    require(type(value) in (int, Q), 'exact rational required; no float or Boolean')
    return Q(value)


def population(p=None):
    if p is None:
        return tuple(Q(1, N) for _ in BASIS)
    require(type(p) is tuple and len(p) == N, 'complete 66-column population required')
    p = tuple(rational(x) for x in p)
    require(all(x >= 0 for x in p) and sum(p) == 1, 'normalized nonnegative population required')
    return p


def parameters(ports, cap, seed_scans, retries):
    thermal._protocol(ports, cap, False)
    require(cap is not None, 'the executed controller requires a finite waiting cap')
    integer(seed_scans, 1, 128, 'seed scans')
    integer(retries, 1, 128, 'retries')


def success_tables(seed_scans=4, retries=32):
    r = census.success_probabilities(seed_scans, retries, mode='restore')
    m = census.success_probabilities(seed_scans, retries, mode='measure')
    require(len(r) == len(m) == N and all(0 < x <= 1 for x in r+m),
            'invalid current-census success table')
    return r, m


def _compose(body, p, restore, measure):
    """Marginalize complete fine census records to their explicit public flags."""
    p = population(p)
    require(len(restore) == len(measure) == N, 'complete census columns required')
    restore, measure = tuple(map(rational, restore)), tuple(map(rational, measure))
    require(all(0 <= x <= 1 for x in restore+measure), 'census success outside [0,1]')
    require(body.denominator > 0 and body.denominator % N == 0, 'invalid body denominator')
    # The inherited body stores uniform-input joint endpoint counts. Checking
    # every column permits this exact reuse for arbitrary initial populations.
    require(all(sum(matrix[y][x] for matrix in body.counts for y in range(N))
                == body.denominator//N for x in range(N)), 'body lost an input column')
    out = defaultdict(Q)
    out['initial_census_failure', ()] = sum((p[x]*(1-restore[x]) for x in range(N)), Q())
    for (kind, word), matrix in zip(body.records, body.counts, strict=True):
        if kind == 'timeout':
            out['waiting_timeout', word] += sum((Q(N, body.denominator)*p[x]*restore[x]
                *sum(matrix[y][x] for y in range(N)) for x in range(N)), Q())
        else:
            require(kind == 'completed', 'unrecognized body record')
            for flag in ('terminal_census_failure', 'return', 'different'):
                out[flag, word] += Q()
            for y, row in enumerate(matrix):
                for x, count in enumerate(row):
                    if not count or not p[x]:
                        continue
                    mass = Q(N*count, body.denominator)*p[x]*restore[x]
                    out['terminal_census_failure', word] += mass*(1-measure[y])
                    out['return' if x == y else 'different', word] += mass*measure[y]
    require(all(x >= 0 for x in out.values()) and sum(out.values()) == 1,
            'complete census/body composition lost probability')
    return dict(sorted(out.items()))


def record_law(p=None, *, ports=(0, 1), cap=2, seed_scans=4, retries=32):
    parameters(ports, cap, seed_scans, retries)
    body = thermal.loop_law(ports, cap=cap)
    return _compose(body, population(p), *success_tables(seed_scans, retries))


def ideal_census_law(p=None, *, ports=(0, 1), cap=2):
    parameters(ports, cap, 1, 1)
    return _compose(thermal.loop_law(ports, cap=cap), population(p), (Q(1),)*N, (Q(1),)*N)


@lru_cache(maxsize=1)
def attempt_rows():
    """All 3960 actual fair seam/coin branches, in input-row convention."""
    positions = {x: i for i, x in enumerate(BASIS)}
    rows = []
    for state in BASIS:
        counts = defaultdict(int)
        for edge in census.quantum.seams():
            for coin in (0, 1):
                out, _ = census.native.native(state, edge, coin)
                counts[positions[out]] += 1
        require(sum(counts.values()) == 60, 'native attempt branch lost mass')
        rows.append(tuple(sorted(counts.items())))
    return tuple(rows)


def attempt_record_law(p=None, *, length=2, seed_scans=4, retries=32, ideal=False):
    """Fixed global attempts with full configuration words from current records.

    This is a separate clock from the incident body. Zero-probability full
    configuration words are included for the executed one/two-attempt scope.
    """
    integer(length, 1, 2, 'executed full-word length')
    require(type(ideal) is bool, 'explicit Boolean ideal-census control required')
    p = population(p)
    restore, measure = ((Q(1),)*N, (Q(1),)*N) if ideal else success_tables(seed_scans, retries)
    active = {(x, x, ()): 1 for x in range(N)}
    for _ in range(length):
        following = defaultdict(int)
        for (origin, state, word), count in active.items():
            for final, weight in attempt_rows()[state]:
                following[origin, final, (*word, final)] += count*weight
        active = following
    require(sum(active.values()) == N*60**length, 'native full-word law lost a branch')
    from itertools import product
    out = {(flag, word): Q() for word in product(range(N), repeat=length)
           for flag in ('terminal_census_failure', 'return', 'different')}
    out['initial_census_failure', ()] = sum((p[x]*(1-restore[x]) for x in range(N)), Q())
    for (origin, final, word), count in active.items():
        mass = p[origin]*restore[origin]*Q(count, 60**length)
        out['terminal_census_failure', word] += mass*(1-measure[final])
        out['return' if origin == final else 'different', word] += mass*measure[final]
    require(sum(out.values()) == 1, 'global-attempt census composition lost probability')
    return dict(sorted(out.items()))


def tv(p, q):
    return sum((abs(p.get(k, Q())-q.get(k, Q())) for k in set(p)|set(q)), Q())/2


def conditioned(law, flag='return'):
    selected = {word: value for (kind, word), value in law.items() if kind == flag}
    total = sum(selected.values(), Q())
    require(total > 0, 'conditioning event has zero probability')
    return {key: value/total for key, value in selected.items()}


def bounds(*, length=2, cap=32, seed_scans=4, retries=32, warmup=None, absorb=False,
           clock='incident'):
    integer(length, 1, 4096, 'loop length')
    integer(cap, 1, 4096, 'waiting cap')
    integer(seed_scans, 1, 128, 'seed scans'); integer(retries, 1, 128, 'retries')
    require(type(absorb) is bool, 'explicit Boolean covering-sweep choice required')
    require(type(clock) is str and clock in ('incident', 'attempt'), 'incident or attempt clock required')
    if warmup is not None:
        integer(warmup, 0, 4096, 'warmup attempts')
    r = census.budget(seed_scans, retries, mode='restore')
    m = census.budget(seed_scans, retries, mode='measure')
    ec = min(Q(1), r['failure_upper']+m['failure_upper'])
    ew = Q() if warmup is None else burn_in.binary_bound(warmup)
    wait = 1-(1-Q(7, 10)**cap)**length if clock == 'incident' else Q()
    delta = min(Q(1), ec+ew+wait)
    p = Q(1, N)
    conditional = min(Q(1), 2*delta/(p-delta)) if delta < p else Q(1)
    return {'clock': clock, 'census_error_upper': ec, 'waiting_error_upper': wait,
            'preparation_error_upper': ew, 'complete_law_error_upper': delta,
            'return_probability_lower': p,
            'conditional_return_error_upper': conditional,
            'exact_uniform_conditional_error_upper': min(Q(1), (ec+wait)/p),
            'native_attempts_max': (30 if absorb else 0)+(warmup or 0)
                +r['native_attempts_max']+length*(cap if clock == 'incident' else 1)+m['native_attempts_max'],
            'exact_uniform_preparation_assumed': warmup is None,
            'binary_input_required_unless_covering_sweep': not absorb}


@dataclass(frozen=True)
class Execution:
    flag: str
    word: tuple
    initial_census: census.Execution | None
    terminal_census: census.Execution | None
    preparation_records: tuple
    body_records: tuple
    output: tuple
    attempts: int
    clock: str = 'incident'


def execute(initial, coins, seams, *, ports=(0, 1), cap=2, seed_scans=4,
            retries=32, warmup=0, absorb=False, law='C', clock='incident', length=2):
    """Run the hidden source with policies receiving only declared records.

    `initial` initializes the source executor, never the current-census policy.
    The probability theorem requires fresh fair coins and IID uniform body and
    warmup seam choices. Deterministic streams are accepted for path replay.
    """
    parameters(ports, cap, seed_scans, retries)
    require(type(clock) is str and clock in ('incident', 'attempt'), 'incident or attempt clock required')
    integer(length, 1, 4096, 'loop length')
    integer(warmup, 0, 4096, 'warmup attempts')
    require(type(absorb) is bool, 'explicit Boolean covering-sweep choice required')
    require(type(law) is str and law in ('C', 'E'), 'C or E fine record law required')
    state = census._configuration(initial, binary=not absorb)
    coin_stream, seam_stream = iter(coins), iter(seams)
    edges = tuple(census.quantum.seams())
    prep, body = [], []

    def coin():
        try:
            value = next(coin_stream)
        except StopIteration as exc:
            raise ValueError('coin stream ended before complete outcome') from exc
        return integer(value, 0, 1, 'placement coin')

    def edge():
        try:
            value = next(seam_stream)
        except StopIteration as exc:
            raise ValueError('seam stream ended before complete outcome') from exc
        return edges[integer(value, 0, 29, 'seam index')]

    def step(e, records):
        nonlocal state
        bit = coin()
        before = state
        state, (_, current) = census.native.native(state, e, bit)
        label = current if law == 'C' else before[e[0]]-before[e[1]]
        records.append((e, bit, current, label))
        return current

    if absorb:
        for e in edges:
            step(e, prep)
        require(max(state) <= 1, 'covering sweep failed to remove collisions')
    for _ in range(warmup):
        step(edge(), prep)
    first = census.execute(state, coin_stream, seed_scans, retries, mode='restore', law=law)
    state = first.output
    if not first.outcome.success:
        return Execution('initial_census_failure', (), first, None, tuple(prep), (), state,
                         len(prep)+first.outcome.attempts, clock)
    # These are records produced by the first controller, not hidden source loads.
    original = first.outcome.initial
    field = tuple(original[p] for p in core.source().region)
    tracked = list(original)
    word = []
    for port in (ports if clock == 'incident' else range(length)):
        for _ in range(cap if clock == 'incident' else 1):
            e = edge()
            current = step(e, body)
            if clock == 'attempt':
                tracked[e[0]] -= current
                tracked[e[1]] += current
                word.append(BASIS.index(tuple(tracked)))
                break
            incident = bool(set(e).intersection(core.source().region))
            if incident:
                field = core.update_field(field, e, current, region=core.source().region)
                word.append(field[core.source().region.index(port)])
                break
        else:
            return Execution('waiting_timeout', tuple(word), first, None, tuple(prep), tuple(body),
                             state, len(prep)+first.outcome.attempts+len(body), clock)
    last = census.execute(state, coin_stream, seed_scans, retries, mode='measure', law=law)
    flag = ('return' if last.outcome.initial == original else 'different') if last.outcome.success else 'terminal_census_failure'
    return Execution(flag, tuple(word), first, last, tuple(prep), tuple(body), last.output,
                     len(prep)+first.outcome.attempts+len(body)+last.outcome.attempts, clock)


def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    return value


def digest(value):
    return hashlib.sha256(json.dumps(encode(value), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def inherited_pins():
    results = {}
    for receipt in ('model/native_clock/configuration_results.json',
                    'model/observer_time_fields/thermal_results.json',
                    'model/record_laws/burn_in_results.json'):
        stored = strict_json.load(ROOT/receipt)
        for path, expected in stored['source_sha256'].items():
            require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected,
                    f'inherited source binding changed: {path}')
        results[receipt] = len(stored['source_sha256'])
    return results


def report():
    pins = inherited_pins()
    stored = strict_json.load(ROOT/'model/native_clock/configuration_results.json')
    for mode, key in (('restore', 'binary_restoration'), ('measure', 'binary_measurement')):
        require(census.certificate(mode=mode) == stored[key], 'native census certificate changed')
    require(census.source_audit() == stored['source_kraus_audit'], 'source current-inverse audit changed')
    body = thermal.loop_law(cap=2)
    require(body == thermal.native_oracle(cap=2), 'direct native history oracle disagrees')
    actual, ideal = record_law(), ideal_census_law()
    first, last = success_tables()
    ec = 1-min(first)+1-min(last)
    require(tv(actual, ideal) <= ec, 'complete census error exceeds uniform bound')
    require(any(key[0] == 'waiting_timeout' and value for key, value in actual.items()), 'waiting failures lost')
    require(any(key[0] == 'terminal_census_failure' and value for key, value in actual.items()), 'terminal failures lost')
    require(actual['initial_census_failure', ()] > 0, 'initial failures lost')
    # At finite resources the successful census can reweight the return word.
    coarse = record_law(seed_scans=1, retries=2)
    bias = tv(conditioned(coarse), conditioned(ideal))
    require(bias > 0, 'census conditioning-bias control disappeared')
    spectral = core.spectral_certificate()
    require(spectral['effective_transfer_lower'] == Q(3, 38), 'source positivity enclosure changed')
    closure = thermal.thermal_cycle(2)
    boundary = thermal.boundary_certificate(2)
    local = thermal.local_closure_control()
    budget = bounds(warmup=512, absorb=True)
    attempts, ideal_attempts = attempt_record_law(), attempt_record_law(ideal=True)
    require(tv(attempts, ideal_attempts) <= ec, 'global-attempt census error exceeds bound')
    P = [[dict(attempt_rows()[i]).get(j, 0) for j in range(N)] for i in range(N)]
    trace2 = sum(P[i][j]*P[j][i] for i in range(N) for j in range(N))
    p_return = sum(value for (flag, _), value in ideal_attempts.items() if flag == 'return')
    require(p_return == Q(trace2, N*60**2), 'native global return law differs from transfer trace')
    paths = [str(path.relative_to(ROOT)) for path in (HERE/'certificate.py', HERE/'PROOF.md', HERE/'test_certificate.py')]
    paths += ['model/native_clock/configuration.py', 'model/native_clock/configuration_results.json',
              'model/observer_time_fields/core.py', 'model/observer_time_fields/thermal.py',
              'model/observer_time_fields/thermal_results.json', 'model/observer/native.py',
              'model/record_laws/burn_in.py', 'model/record_laws/burn_in_results.json']
    return encode({'schema': 'oph.native_thermal_return.v1', 'all_pass': True,
        'source_pins_rechecked': pins,
        'native_census': {'binary_inputs': N, 'successful_seed_cases_per_mode': 600,
                         'initial_mode': 'restore', 'terminal_mode': 'measure',
                         'first_success_table_sha256': digest(first), 'terminal_success_table_sha256': digest(last),
                         'exact_combined_failure_upper': ec, 'simple_combined_failure_upper': Q(17, 2**29)},
        'exact_composition': {'body_intervals': 2, 'body_attempt_cap_per_interval': 2,
                              'public_outcomes_including_zero': len(actual), 'complete_probability': sum(actual.values()),
                              'complete_output_sha256': digest(actual), 'complete_census_distance': tv(actual, ideal),
                              'finite_success_conditioning_changes_return_words': True,
                              'conditioning_bias_control_sha256': digest(bias)},
        'source_transfer': spectral, 'ideal_closed_history': closure,
        'exact_KMS_boundary': boundary, 'local_return_negative_control': local,
        'finite_preparation_and_readout_budget': budget,
        'complete_attempt_variant': {'body_length': 2, 'public_outcomes_including_zero': len(attempts),
             'output_sha256': digest(attempts), 'complete_census_distance': tv(attempts, ideal_attempts),
             'ideal_return_probability': p_return,
             'budget': bounds(clock='attempt', warmup=512, absorb=True),
             'full_configuration_words_from_initial_census_and_global_signed_currents': True,
             'incident_transfer_or_incident_clock_used': False},
        'scope': {'initial_configuration_oracle_removed_to_certified_finite_error': True,
                  'initial_regional_bits_derived_by_first_census': True,
                  'two_censuses_suffice_for_all_admitted_regional_field_reads': True,
                  'all_failures_and_private_tapes_retained': True,
                  'law_is_declared_marginal_ignoring_private_census_and_preparation_tapes': True,
                  'uniform_seam_and_fair_coin_sampling_supplied': True,
                  'H_reconstructed_from_native_transfer_before_loop_length': True,
                  'internal_inverse_temperature_is_supplied_loop_length': True,
                  'native_quantum_Gibbs_state_prepared': False,
                  'full_complex_effect_quotient_derived': False,
                  'physical_energy_temperature_or_time_units_derived': False},
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = report()
    target = HERE/'results.json'
    if args.write:
        target.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check:
        require(strict_json.load(target) == result, 'native thermal return receipt changed')
    print(json.dumps({'all_pass': True, 'binary_inputs': N,
                      'public_outcomes': result['exact_composition']['public_outcomes_including_zero'],
                      'physical_temperature_derived': False}))


if __name__ == '__main__':
    main()
