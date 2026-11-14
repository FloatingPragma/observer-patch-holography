"""A fixed, exact decision protocol; no data are used to select its settings.

The finite C/E hypotheses differ in their retained quantum instruments.
Fresh independent trials, a common preparation/readout identification and an
independently established apparatus-error budget are premises. Count files
cannot authenticate these premises or their own experimental provenance.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import random

from . import observer_instrument as instrument

ROOT = Path(__file__).resolve().parents[2]
TRIALS = 8192
HALF_WIDTH = F(1, 32)
SERVICE_ERROR = F(1, 64)
DELTA = F(1, 1000)
SETTINGS = ('phase_0', 'phase_pi', 'phase_pi_over_2', 'diagonal_control')
BINS = tuple((coin, current, outcome) for coin in (0, 1)
             for current in (-1, 0, 1) for outcome in (-1, 1))
SOURCE_FILES = ('model/observer/native.py',
                'model/predictions/observer_instrument.py',
                'model/predictions/protocol.py')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def key(coin, current, outcome):
    return f'{coin},{current},{outcome}'


def concentration_certificate():
    """Rational enclosure of 2*m*exp(-2*n*b²), with no rounded logarithm.

    Hoeffding applies to the indicator of every bin, abort, and aggregate
    plus outcome. A union bound makes their intervals simultaneous. Testing
    several hypotheses against these same intervals adds no multiplicity.
    """
    m = len(SETTINGS)*(len(BINS)+2)
    exponent = 2*TRIALS*HALF_WIDTH**2
    term = total = F(1)
    for order in range(1, 33):
        term *= exponent/order
        total += term
    bound = 2*m/total  # exp(exponent) >= positive Taylor partial sum
    if bound >= DELTA:
        raise ArithmeticError('fixed sample plan does not meet its error bound')
    return {'simultaneous_statistics': m, 'exponential_argument': str(exponent),
            'exp_lower_partial_sum': str(total), 'taylor_degree': 32,
            'familywise_error_upper': str(bound), 'required_upper': str(DELTA)}


def build_protocol():
    predictions = {}
    for candidate in instrument.CANDIDATES:
        predictions[candidate] = {}
        for setting in SETTINGS:
            law = instrument.public_distribution(candidate, setting)
            if sum(law.values()) != 1:
                raise ArithmeticError('instrument probabilities do not sum to one')
            predictions[candidate][setting] = {
                'plus': str(instrument.prediction(candidate, setting)),
                'abort': '0',
                'joint': {key(*event): str(law.get(event, F(0))) for event in BINS}}
    return {
        'schema': 'oph.observer.discrimination.protocol.v1',
        'scope': 'finite quantum observer instrument C versus E; conditional apparatus test',
        'settings': list(SETTINGS), 'trials_per_setting': TRIALS,
        'sampling_half_width': str(HALF_WIDTH),
        'apparatus_probability_error_budget': str(SERVICE_ERROR),
        'concentration': concentration_certificate(), 'predictions': predictions,
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in SOURCE_FILES},
        'trial_contract': [
            'fresh independent trials at the fixed sample size for every setting; no optional stopping',
            'common phase preparation and source-projector readout fixed independently of validation data',
            'coin and current are committed records; E contrast is retained internally and only summed for the shared interface',
            'all plus, minus and abort outcomes count; no selection on a record, outcome or successful service',
            'each joint/plus/abort probability differs from its ideal candidate by at most the independently established apparatus budget',
            'unconditional joint frequencies are scored; no inference of apparatus calibration or provenance from a passing count file'],
        'meaning': 'fixed basis-state preparation, seam (0,1), original current record, binary source-projector measurement',
        'physical_verdict': 'no independently authenticated experimental dataset is supplied',
    }


def checked_protocol(path=None):
    path = path or Path(__file__).with_name('protocol.json')
    retained = json.loads(Path(path).read_text(encoding='utf-8'))
    if retained != build_protocol():
        raise ValueError('protocol differs from the retained source-bound specification')
    return retained


def data_template(protocol):
    return {
        'schema': 'oph.observer.discrimination.counts.v1',
        'protocol_sha256': digest(protocol),
        'data_kind': 'external_unverified',
        'provenance': 'Replace with acquisition identifier and independent calibration references.',
        'settings': [{'setting': setting, 'aborts': None,
                      'events': [{'coin': c, 'current': j, 'outcome': y, 'count': None}
                                 for c, j, y in BINS]} for setting in SETTINGS],
    }


def _count(value):
    if type(value) is not int or value < 0:
        raise ValueError('every outcome requires an explicit nonnegative integer count')
    return value


def _validate_data(data, protocol):
    if not isinstance(data, dict) or set(data) != {
            'schema', 'protocol_sha256', 'data_kind', 'provenance', 'settings'}:
        raise ValueError('incorrect counts schema or undeclared fields')
    if data['schema'] != 'oph.observer.discrimination.counts.v1':
        raise ValueError('unsupported counts schema')
    if data['protocol_sha256'] != digest(protocol):
        raise ValueError('counts belong to a different protocol')
    if data['data_kind'] not in ('synthetic', 'external_unverified'):
        raise ValueError('data must be explicitly synthetic or externally supplied and unverified')
    if not isinstance(data['provenance'], str) or not data['provenance'].strip():
        raise ValueError('data provenance description required')
    if not isinstance(data['settings'], list) or len(data['settings']) != len(SETTINGS):
        raise ValueError('all four fixed settings are required')
    parsed = {}
    for row in data['settings']:
        if not isinstance(row, dict) or set(row) != {'setting', 'aborts', 'events'}:
            raise ValueError('incorrect setting schema')
        setting = row['setting']
        if setting not in SETTINGS or setting in parsed:
            raise ValueError('unknown or duplicated setting')
        aborts = _count(row['aborts'])
        if not isinstance(row['events'], list) or len(row['events']) != len(BINS):
            raise ValueError('all committed-record and outcome bins must be present')
        counts = {}
        for event in row['events']:
            if not isinstance(event, dict) or set(event) != {'coin', 'current', 'outcome', 'count'}:
                raise ValueError('incorrect event schema')
            values = tuple(event[name] for name in ('coin', 'current', 'outcome'))
            if any(type(v) is not int for v in values) or values not in BINS or values in counts:
                raise ValueError('unknown or duplicated committed-record/outcome bin')
            counts[values] = _count(event['count'])
        if sum(counts.values())+aborts != TRIALS:
            raise ValueError('fixed trial count not met; aborts and every outcome must be retained')
        parsed[setting] = (counts, aborts)
    return parsed


def _interval(center, radius):
    return max(F(0), center-radius), min(F(1), center+radius)


def evaluate(data=None, protocol=None):
    protocol = protocol if protocol is not None else checked_protocol()
    if protocol != build_protocol():
        raise ValueError('only the fixed source-bound protocol can be evaluated')
    result = {'schema': 'oph.observer.discrimination.verdict.v1',
              'protocol_sha256': digest(protocol),
              'physical_scope': 'conditional finite-instrument comparison; no OPH-wide verdict',
              'provenance_verified': False, 'empirical_confirmation': False}
    if data is None:
        return {**result, 'status': 'no_data', 'candidate_verdicts': {}}
    parsed = _validate_data(data, protocol)
    observations = {}
    budget_exceeded = []
    for setting, (counts, aborts) in parsed.items():
        means = {key(*event): F(count, TRIALS) for event, count in counts.items()}
        means['plus'] = F(sum(count for (_, _, y), count in counts.items() if y == 1), TRIALS)
        means['abort'] = F(aborts, TRIALS)
        observations[setting] = {name: list(map(str, _interval(mean, HALF_WIDTH)))
                                 for name, mean in means.items()}
        if means['abort']-HALF_WIDTH > SERVICE_ERROR:
            budget_exceeded.append(setting)
    verdicts = {}
    for candidate, settings in protocol['predictions'].items():
        failures = []
        for setting, prediction in settings.items():
            means = {**prediction['joint'], 'plus': prediction['plus'], 'abort': '0'}
            for statistic, center in means.items():
                lo, hi = _interval(F(center), SERVICE_ERROR)
                observed_lo, observed_hi = map(F, observations[setting][statistic])
                if observed_hi < lo or observed_lo > hi:
                    failures.append({'setting': setting, 'statistic': statistic,
                                     'model_interval': [str(lo), str(hi)],
                                     'observed_interval': [str(observed_lo), str(observed_hi)]})
        verdicts[candidate] = {'verdict': 'rejected_under_contract' if failures else 'not_rejected',
                               'separating_statistics': failures}
    return {**result, 'status': 'apparatus_budget_exceeded' if budget_exceeded else 'comparison_complete',
            'data_kind': data['data_kind'], 'provenance': data['provenance'],
            'apparatus_budget_failures': budget_exceeded,
            'candidate_verdicts': verdicts, 'simultaneous_observation_intervals': observations,
            'interpretation': 'An exceeded apparatus budget prevents attributing rejection to the ideal rule; otherwise verdicts require every declared preparation, readout, independence and service premise. Not rejected does not mean established.'}


def synthetic_counts(candidate, seed, protocol=None):
    """Labeled software negative controls; never experimental evidence."""
    if candidate not in instrument.CANDIDATES or type(seed) is not int:
        raise ValueError('declared candidate and integer seed required')
    protocol = protocol if protocol is not None else checked_protocol()
    data = data_template(protocol)
    data['data_kind'] = 'synthetic'
    data['provenance'] = f'Software draws from instrument {candidate}; seed {seed}; no laboratory data.'
    rng = random.Random(seed)
    for row in data['settings']:
        probabilities = instrument.public_distribution(candidate, row['setting'])
        denominator = 1
        from math import lcm
        for chance in probabilities.values():
            denominator = lcm(denominator, chance.denominator)
        bins = []
        for event in BINS:
            bins.extend([event]*int(probabilities.get(event, F(0))*denominator))
        if len(bins) != denominator:
            raise ArithmeticError('incomplete synthetic event law')
        counts = dict.fromkeys(BINS, 0)
        for _ in range(TRIALS):
            counts[bins[rng.randrange(denominator)]] += 1
        row['aborts'] = 0
        for event in row['events']:
            event['count'] = counts[event['coin'], event['current'], event['outcome']]
    return data


def run():
    protocol = checked_protocol()
    controls = {}
    for candidate, seed in (('C', 2026101201), ('E', 2026101202)):
        data = synthetic_counts(candidate, seed, protocol)
        verdict = evaluate(data, protocol)
        for hypothesis in instrument.CANDIDATES:
            expected = 'not_rejected' if hypothesis == candidate else 'rejected_under_contract'
            if verdict['candidate_verdicts'][hypothesis]['verdict'] != expected:
                raise AssertionError('synthetic discrimination control failed')
        controls[candidate] = {'data': data, 'verdict': verdict}
    from .endpoint_audit import audit
    return {'schema': 'oph.observer.prediction_checks.v1', 'protocol_sha256': digest(protocol),
            'no_data_verdict': evaluate(None, protocol),
            'synthetic_controls': controls, 'retrospective_endpoint_audit': audit(),
            'software_checks_pass': True, 'empirical_confirmation': False}
