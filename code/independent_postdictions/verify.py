"""Independent direct sums, quadrature, algebraic tau and exact rule checks.

No imports from this packet's root solver, emitter, decision code or either
historical producer. A custody match alone never accepts scientific evidence.
"""
from decimal import Decimal
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
from urllib.parse import parse_qs, urlparse

from mpmath import MPContext
import pin_scope


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CLAIM = 'OPH-INDEPENDENT-ALPHA-TAU-REPRODUCTION'
CERT = 'code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json'
TRUNK = 'code/P_derivation/runtime/p_closure_trunk_current.json'
TAU = 'code/particles/runs/leptons/koide_balance_comparison.json'
FROZEN = 'evidence/custody/falsification/frozen_targets/fz10_2026-07-28/'
OWN = ('__init__.py', 'replay.py', 'comparison.py', 'admission.py', 'fetch_catalogue.py',
       'build.py', 'verify.py', 'inner_check.py', 'test_reproduction.py', 'inputs.json', 'protocol.json',
       'catalogue_discovery.json', 'README.md', 'CONTRACT.md')
SOURCES = ['code/independent_postdictions/'+p for p in OWN]+[
    CERT, TRUNK, TAU, 'code/P_derivation/paper_math.py',
    'code/P_derivation/codata_2022_alpha_fixture.json',
    'code/P_derivation/runtime/selection_accounting.json', 'tracking/null_model_scorecard.md',
    'code/particles/leptons/koide_balance_comparison_certificate.py',
    'paper/deriving_the_particle_zoo_from_observer_consistency.tex',
    FROZEN+'frozen_target_koide_conditional_tau_2026-07-28.md',
    FROZEN+'koide_balance_comparison_frozen_2026-07-28.json',
    'docs/research/INDEPENDENT_POSTDICTION_COMPARISON.md',
    '.github/workflows/independent-postdictions.yml', 'requirements.txt', '.gitattributes']
SOURCES = pin_scope.bound(SOURCES)

# Normative, reviewed contracts: rebuilding a receipt cannot loosen these.
# A policy change must separately update this review pin and its explanation.
POLICY_PINS = {
    'protocol.json': '7802a3dec6293d6979f51903aa93ae319d2f60f801f66755e859c4e938be8316',
    'inputs.json': '0cd8c7944a48c790fb9f499cf874471d30a73ed1aabaadbc471b2c2746be61cb',
    # Human classification is of this exact historical metadata, not a
    # keyword classifier that can safely label a later catalogue response.
    'catalogue_discovery.json': '84a021209fa5e7d9cb4971ac3df86dda4b8e655e2e0e955932a0aec168806059',
}
HISTORICAL_PINS = {
    CERT: '20ab213f836fbe95149d4b40ba57aa144038580a0acb899322f4b72d963c2ad6',
    FROZEN+'koide_balance_comparison_frozen_2026-07-28.json':
        '09efbad8de813790d10671b425daef69429831c58daa212f4be131c61f653898',
}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def keys(value, expected):
    need(type(value) is dict and value.keys() == set(expected.split()), 'exact object fields')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('ascii')


def same(actual, expected):
    need(type(actual) is type(expected), 'exact JSON type')
    if type(expected) is dict:
        need(actual.keys() == expected.keys(), 'exact fields')
        for k in expected:
            same(actual[k], expected[k])
    elif type(expected) is list:
        need(len(actual) == len(expected), 'complete census')
        for a, b in zip(actual, expected):
            same(a, b)
    else:
        need(actual == expected, 'exact semantic value')


def load(path):
    need(Path(path).stat().st_size < 200_000, 'bounded primary input')
    raw = Path(path).read_bytes()
    need(len(raw) < 200_000, 'bounded primary input')
    def unique(pairs):
        out = {}
        for k, v in pairs:
            need(k not in out, 'duplicate JSON key')
            out[k] = v
        return out
    def forbidden(x):
        raise ValueError('nonfinite JSON token')
    def integer(x):
        need(len(x) < 50, 'bounded JSON integer')
        return int(x)
    def floating(x):
        need(len(x) < 50 and math.isfinite(float(x)), 'bounded finite JSON float')
        return float(x)
    try:
        result = json.loads(raw.decode('utf-8'), object_pairs_hook=unique,
                            parse_constant=forbidden, parse_int=integer, parse_float=floating)
    except (RecursionError, UnicodeError) as exc:
        raise ValueError('invalid JSON encoding or nesting') from exc
    pending, count = [(result, 0)], 0
    while pending:
        value, depth = pending.pop()
        count += 1
        need(depth <= 32 and count <= 10_000, 'bounded JSON structure')
        children = value.values() if type(value) is dict else value if type(value) is list else ()
        pending.extend((child, depth+1) for child in children)
    return result


def pins():
    result = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(SOURCES)}
    registry = json.loads((ROOT/'claims/claim_registry.yaml').read_text(encoding='utf-8'))
    claims = [x for x in registry['claims'] if x['claim_id'] == CLAIM]
    need(len(claims) == 1, 'one reviewed claim')
    definitions = (ROOT/'claims/assumption_dictionary.md').read_text(encoding='utf-8').splitlines()
    rows = []
    for name in claims[0]['assumptions']:
        matches = [x for x in definitions if x.startswith('| `'+name+'` |')]
        need(len(matches) == 1, 'one consumed assumption definition')
        rows.extend(matches)
    targets = json.loads((ROOT/'claims/frozen_prediction_register.json').read_text())['rows']
    targets = [r for r in targets if r['id'] in ('FZ-10', 'FZ-11')]
    need(len(targets) == 2, 'canonical target projection')
    result['projection:claim_and_assumptions'] = hashlib.sha256(canonical([claims[0], rows])).hexdigest()
    result['projection:canonical_FZ10_FZ11'] = hashlib.sha256(canonical(targets)).hexdigest()
    return result


def numeric(c, value):
    need(type(value) is str and len(value) < 160 and
         re.fullmatch(r'-?[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]{1,3})?', value) is not None,
         'bounded decimal evidence')
    out = c.mpf(value)
    need(c.isfinite(out) and abs(out) < c.mpf('1e10'), 'finite bounded number')
    return out


def verify_policy():
    for name, digest in HISTORICAL_PINS.items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest,
             'immutable historical evidence')
    for name, digest in POLICY_PINS.items():
        need(hashlib.sha256(canonical(load(HERE/name))).hexdigest() == digest,
             'reviewed policy or reference-data contract changed')
    snapshot = load(HERE/'catalogue_discovery.json')
    keys(snapshot, 'limit outcomes_requested queries retrieved_utc scope')
    same(snapshot['limit'], 25)
    same(snapshot['outcomes_requested'], False)
    need(type(snapshot['queries']) is list and len(snapshot['queries']) == 4, 'bounded discovery menu')
    expected_queries = [
        'date > 2026-07-28 and title tau and title mass',
        'date > 2026-07-28 and (title tau or title \u03c4 or title lepton) and title mass',
        'date > 2023-01-01 and date < 2024-01-01 and title tau and title mass',
        'date > 2023-01-01 and date < 2024-01-01 and title lepton and title mass']
    fields = 'titles,arxiv_eprints,dois,preprint_date,document_type,collaborations'
    for row, query, count in zip(snapshot['queries'], expected_queries, (0, 4, 0, 33)):
        keys(row, 'endpoint query rows total')
        same(row['query'], query)
        same(row['total'], count)
        url = urlparse(row['endpoint'])
        same((url.scheme, url.netloc, url.path), ('https', 'inspirehep.net', '/api/literature'))
        same(parse_qs(url.query), dict(q=[query], fields=[fields], size=['25'], sort=['mostrecent']))
        need(type(row['rows']) is list and len(row['rows']) == min(count, 25), 'retained metadata count')
        ids = set()
        for record in row['rows']:
            keys(record, 'id created '+fields.replace(',', ' '))
            need(type(record['id']) is str and record['id'].isdigit() and record['id'] not in ids, 'unique catalogue identifier')
            ids.add(record['id'])
    # These are a bounded human classification of the retained metadata,
    # not a theorem that the world's unindexed data are absent.
    future = {r['id']: r for r in snapshot['queries'][1]['rows']}
    same(set(future), {'3200778', '3193587', '3096034', '3193022'})
    for identity, title in [('3200778', 'vector-like top partner'), ('3193587', 'Koide')]:
        need(any(title in t['title'] for t in future[identity]['titles']), 'retained nonmeasurement classification')
    for identity in ('3096034', '3193022'):
        need(future[identity]['preprint_date'] < '2026-07-28', 'first-public date overrides later index date')
    control = [r for r in snapshot['queries'][3]['rows'] if r['id'] == '2663717']
    need(len(control) == 1 and any(a['value'] == '2305.19116' for a in control[0]['arxiv_eprints']),
         'historical measurement positive control')


def close(c, actual, expected, tolerance='1e-35'):
    need(abs(numeric(c, actual)-expected) < c.mpf(tolerance), 'independent numerical replay')


def verify_alpha(data):
    keys(data, 'rows electron_menu masses_over_v')
    c = MPContext()
    c.dps = 60
    roots = sorted([1+c.sqrt(2)*c.cos(c.mpf(2)/9+2*k*c.pi/3) for k in range(3)])
    need(type(data['electron_menu']) is list and len(data['electron_menu']) == 7, 'complete exponent menu')
    scores = []
    for row, n in zip(data['electron_menu'], range(5, 12)):
        keys(row, 'exponent score')
        same(row['exponent'], n)
        powers = (n, 4, 3)
        score = max(abs(c.log((roots[i]/roots[j])**2/(c.mpf(1)/6)**(powers[i]-powers[j])))
                    for i, j in itertools.combinations(range(3), 2))
        close(c, row['score'], score)
        scores.append((score, n))
    need(min(scores)[1] == 7, 'electron selector')
    factor = c.exp(c.log(2)/6-c.fsum(c.log(r*r*c.sqrt(2)*6**n)
                     for r, n in zip(roots, (7, 4, 3)))/3)
    masses = dict(zip(('e', 'mu', 'tau'), [factor*r*r for r in roots]))
    masses.update({name: (c.mpf(1)/6)**n/c.sqrt(2) for name, n in
                   zip(('u', 'c', 't', 'd', 's', 'b'), (6, 3, 0, 6, 4, 2))})
    need(type(data['masses_over_v']) is dict and data['masses_over_v'].keys() == masses.keys(), 'dimensionless mass census')
    for k, value in masses.items():
        close(c, data['masses_over_v'][k], value)
    cert = load(ROOT/CERT)['modes']
    need(type(data['rows']) is list and len(data['rows']) == 3, 'all three declared maps')
    for row, mode in zip(data['rows'], ('structured', 'gauge_width', 'asymptotic')):
        keys(row, 'mode P alpha_inverse alpha_U mZ_over_v couplings anchor_inverse lepton_transport unscreened_quark_transport printed_residual_tolerance')
        same(row['mode'], mode)
        same(row['printed_residual_tolerance'], '1e-38')
        P, g, z = [numeric(c, row[k]) for k in ('P', 'alpha_U', 'mZ_over_v')]
        need(c.mpf('1.62') < P < c.mpf('1.64') and c.mpf('.03') < g < c.mpf('.05') and c.mpf('.3') < z < c.mpf('.5'), 'physical branch box')
        # Reconstruct dimensionful ratios then run downward, instead of the
        # producer's eliminated log-ratio expression and simultaneous solve.
        v = P**(-c.mpf('.5'))*c.exp(-c.pi/(2*g))
        mu = c.exp(-2*c.pi)*P**(c.mpf(1)/6)
        alphas = [1/(1/g+b*c.log(mu/(z*v))/(2*c.pi)) for b in (c.mpf(33)/5, 1, -3)]
        need(type(row['couplings']) is list and len(row['couplings']) == 3, 'three running couplings')
        for x, y in zip(row['couplings'], alphas):
            close(c, x, y)
        need(abs(z*z-c.pi*(alphas[1]+c.mpf(3)/5*alphas[0])) < c.mpf('1e-38'), 'weak-scale equation')
        sums = []
        for group, t in ((2, 4*c.pi**2*alphas[1]), (3, 4*c.pi**2*alphas[2])):
            weights, numerators = [], []
            labels = ((n,) for n in range(121)) if group == 2 else itertools.product(range(91), repeat=2)
            for label in labels:
                if group == 2:
                    n, = label
                    dimension, casimir = n+1, c.mpf(n*(n+2))/4
                else:
                    p, q = label
                    dimension = (p+1)*(q+1)*(p+q+2)//2
                    casimir = c.mpf(p*p+q*q+p*q+3*p+3*q)/3
                w = dimension*c.exp(-t*casimir)
                weights.append(w)
                numerators.append(w*c.log(dimension))
            sums.append(c.fsum(numerators)/c.fsum(weights))
        need(abs(sum(sums)-P/4) < c.mpf('1e-38'), 'pixel equation')
        def kernel(mass):
            if mode == 'asymptotic':
                return (c.log(z*z/(mass*mass))-c.mpf(5)/3)/(3*c.pi)
            # Direct integral; no closed-form kernel imported or duplicated.
            f = lambda x: x*(1-x)*c.log1p((z/mass)**2*x*(1-x))
            return 4/c.pi*c.quad(f, [0, c.mpf('.000001'), c.mpf('.01'), c.mpf('.5')])
        lepton = c.fsum(kernel(masses[k]) for k in ('e', 'mu', 'tau'))
        quark = c.fsum(kernel(masses[k])*w for k, w in
                       (('u', c.mpf(4)/3), ('c', c.mpf(4)/3), ('d', c.mpf(1)/3), ('s', c.mpf(1)/3), ('b', c.mpf(1)/3)))
        anchor = 1/alphas[1]+5/(3*alphas[0])
        inverse = anchor+lepton+(1-3*alphas[2]/c.pi)*quark+(g if mode == 'gauge_width' else 0)
        for key, value in (('anchor_inverse', anchor), ('lepton_transport', lepton),
                           ('unscreened_quark_transport', quark), ('alpha_inverse', inverse)):
            close(c, row[key], value, '1e-38')
        need(abs((P-(1+c.sqrt(5))/2)*inverse-c.sqrt(c.pi)) < c.mpf('1e-38'), 'outer equation')
        if mode != 'asymptotic':
            key = 'thomson_structured_running'+('_plus_gauge_width' if mode == 'gauge_width' else '')
            for field, value in (('P', P), ('alpha_inv', inverse)):
                bounds = cert[key]['certified_enclosure'][field]
                need(c.mpf(bounds['lo']) <= value <= c.mpf(bounds['hi']), 'agreement with retained certified enclosure')
    return c


def verify_tau(row, inputs):
    keys(row, 'central_mev corners_mev excluded_small_root_mev measured_Q measured_Q_corner_range outward_mev jacobian first_order_sigma_max_mev residual_in_sigma_with_arbitrary_correlations covariance_convention')
    c = MPContext()
    c.dps = 65
    e, m, se, sm = [c.mpf(inputs['codata'][k][field]) for k, field in
                    (('electron', 'value'), ('muon', 'value'), ('electron', 'sigma'), ('muon', 'sigma'))]
    def root(a, b, sign=1):
        s = c.sqrt(a)+c.sqrt(b)
        return (2*s+sign*c.sqrt(6*s*s-3*(a+b)))**2
    prediction = root(e, m)
    close(c, row['central_mev'], prediction)
    close(c, row['excluded_small_root_mev'], root(e, m, -1))
    need(root(e, m, -1) < m < prediction, 'mass-ordering selector')
    corners = [root(e+i*se, m+j*sm) for i, j in itertools.product((-1, 1), repeat=2)]
    need(type(row['corners_mev']) is list and len(row['corners_mev']) == 4, 'complete corner census')
    for a, b in zip(row['corners_mev'], corners):
        close(c, a, b)
    lo, hi = c.floor(min(corners)*10**6)/10**6, c.ceil(max(corners)*10**6)/10**6
    need(type(row['outward_mev']) is list and len(row['outward_mev']) == 2, 'two outward endpoints')
    close(c, row['outward_mev'][0], lo)
    close(c, row['outward_mev'][1], hi)
    legacy = load(ROOT/TAU)['conditional_tau']
    for value, old in zip((lo, hi), legacy['tau_enclosure_mev_outward']):
        need(value == c.mpf(old), 'historical enclosure agreement')
    a, b = c.sqrt(e), c.sqrt(m)
    discriminant = c.sqrt(3*a*a+12*a*b+3*b*b)
    derivatives = [c.sqrt(prediction)/a*(2+(3*a+6*b)/discriminant),
                   c.sqrt(prediction)/b*(2+(6*a+3*b)/discriminant)]
    need(type(row['jacobian']) is list and len(row['jacobian']) == 2, 'both input sensitivities')
    for value, exact in zip(row['jacobian'], derivatives):
        close(c, value, exact, '1e-30')
    bound = derivatives[0]*se+derivatives[1]*sm
    close(c, row['first_order_sigma_max_mev'], bound)
    measured, sigma = c.mpf(inputs['tau']['value']), c.mpf(inputs['tau']['sigma'])
    expected = [abs(prediction-measured)/(sigma+bound), abs(prediction-measured)/(sigma-bound)]
    need(type(row['residual_in_sigma_with_arbitrary_correlations']) is list and len(row['residual_in_sigma_with_arbitrary_correlations']) == 2, 'covariance bound endpoints')
    for value, exact in zip(row['residual_in_sigma_with_arbitrary_correlations'], expected):
        close(c, value, exact)
    Q = lambda a, b, d: (a+b+d)/(c.sqrt(a)+c.sqrt(b)+c.sqrt(d))**2
    close(c, row['measured_Q'], Q(e, m, measured))
    corners_q = [Q(e+i*se, m+j*sm, measured+k*sigma) for i, j, k in itertools.product((-1, 1), repeat=3)]
    need(type(row['measured_Q_corner_range']) is list and len(row['measured_Q_corner_range']) == 2, 'Q range endpoints')
    for value, exact in zip(row['measured_Q_corner_range'], (min(corners_q), max(corners_q))):
        close(c, value, exact)
    same(row['covariance_convention'], 'first-order Cauchy bounds; no independence, Gaussian confidence level or theory-error distribution assumed')
    return prediction


def verify_calculations(data):
    keys(data, 'alpha comparison_diagnostics measured_pixel_diagnostic approximate_trunk tau historical_tau normal_error_reference_diagnostic conclusions')
    c = verify_alpha(data['alpha'])
    inputs = load(HERE/'inputs.json')
    historical = load(ROOT/TAU)['measured_imports']['values_mev']
    for name in ('electron', 'muon'):
        same([inputs['codata'][name]['value'], str(Decimal(inputs['codata'][name]['sigma']))],
             [historical[name][0], str(Decimal(historical[name][1]))])
        same(inputs['codata'][name]['unit'], 'MeV')
    same([inputs['tau'][k] for k in ('value', 'sigma', 'unit')], ['1776.93', '0.09', 'MeV'])
    alpha_ref = load(ROOT/'code/P_derivation/codata_2022_alpha_fixture.json')['inverse_fine_structure_constant']
    same(inputs['codata']['alpha_inverse']['value'], alpha_ref['value'])
    same(inputs['codata']['alpha_inverse']['sigma'], alpha_ref['standard_uncertainty'])
    ref, sigma = c.mpf(alpha_ref['value']), c.mpf(alpha_ref['standard_uncertainty'])
    need(type(data['comparison_diagnostics']) is list and len(data['comparison_diagnostics']) == 3, 'all residuals retained')
    for row, solved in zip(data['comparison_diagnostics'], data['alpha']['rows']):
        keys(row, 'mode signed_reference_minus_output absolute_ppm absolute_in_reference_sigma_only minimum_symmetric_remainder_for_two_sigma_overlap')
        same(row['mode'], solved['mode'])
        d = ref-c.mpf(solved['alpha_inverse'])
        for key, value in (('signed_reference_minus_output', d), ('absolute_ppm', abs(d)/ref*1000000),
                           ('absolute_in_reference_sigma_only', abs(d)/sigma),
                           ('minimum_symmetric_remainder_for_two_sigma_overlap', max(abs(d)-2*sigma, 0))):
            close(c, row[key], value, '1e-30')
    mixed = data['measured_pixel_diagnostic']
    keys(mixed, 'P_measured alpha_U_measured_pixel mixed_inverse')
    pc = (1+c.sqrt(5))/2+c.sqrt(c.pi)/ref
    close(c, mixed['P_measured'], pc)
    # The measured-pixel alpha_U also satisfies the *inner* equations. Reuse
    # only the direct heat formula via a two-variable independent solve.
    from .inner_check import measured_gauge
    g = measured_gauge(pc)
    close(c, mixed['alpha_U_measured_pixel'], g)
    close(c, mixed['mixed_inverse'], c.mpf(data['alpha']['rows'][0]['alpha_inverse'])+g)
    trunk = load(ROOT/TRUNK)['fixed_point_candidate']
    row = data['approximate_trunk']
    keys(row, 'printed_inverse relative_pair_defect printed_minus_converged_asymptotic')
    same(row['printed_inverse'], trunk['alpha_inv'])
    defect = (c.mpf(trunk['P'])-(1+c.sqrt(5))/2)*c.mpf(trunk['alpha_inv'])/c.sqrt(c.pi)-1
    close(c, row['relative_pair_defect'], defect)
    close(c, row['printed_minus_converged_asymptotic'], c.mpf(trunk['alpha_inv'])-c.mpf(data['alpha']['rows'][2]['alpha_inverse']))
    tau = verify_tau(data['tau'], inputs)
    same(data['historical_tau']['frozen_center_verdict'], 'INCONCLUSIVE')
    same(data['historical_tau']['window_robust_verdict'], 'INCONCLUSIVE')
    same(data['historical_tau']['interpretation'], 'conditional balanced-mass relation; no OPH-versus-Koide discrimination')
    keys(data['historical_tau'], 'frozen_center_verdict distance_mev distance_over_reported_sigma window_robust_verdict interpretation')
    close(c, data['historical_tau']['distance_mev'], c.mpf('1776.969027')-c.mpf('1776.93'))
    close(c, data['historical_tau']['distance_over_reported_sigma'], (c.mpf('1776.969027')-c.mpf('1776.93'))/c.mpf('.09'))
    keys(data['normal_error_reference_diagnostic'], 'oph_to_koide twice_log_free_mass_to_balanced')
    same(data['normal_error_reference_diagnostic']['oph_to_koide'], '0')
    close(c, data['normal_error_reference_diagnostic']['twice_log_free_mass_to_balanced'], ((tau-c.mpf('1776.93'))/c.mpf('.09'))**2)
    same(data['conclusions'], dict(reproduction='AGREEMENT_UNDER_DECLARED_CONVENTIONS', natural_comparison='NOT_READY',
         selected_target='FZ-10', new_natural_outcomes=0, global_chance_probability=None,
         corners_are_joint_68_percent_interval=False, tau_independently_tests_P=False,
         oph_versus_koide_log_likelihood_ratio='0'))


def verify(path=HERE/'receipt.json'):
    packet = load(path)
    keys(packet, 'schema sources calculations')
    same(packet['schema'], 'oph-independent-postdictions-v1')
    same(pin_scope.bound_pins(packet['sources']), pins())
    verify_policy()
    verify_calculations(packet['calculations'])
    return packet


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified independent alpha/tau arithmetic, historical comparisons and NOT_READY boundary')
