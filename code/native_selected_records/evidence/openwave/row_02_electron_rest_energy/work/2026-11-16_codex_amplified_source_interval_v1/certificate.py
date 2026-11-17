"""Independent finite-source interval certificate; explicit first-record mode."""
import argparse
import hashlib
import json
from pathlib import Path
import time
from fractions import Fraction as F
from evaluator import *
from service_error import budget, envelope, CONTRACT

HERE = Path(__file__).resolve().parent


def check_counts(source):
    """Bind the weighted source mean to the actual measured eighteen columns."""
    columns = source['native_columns']
    table = {(r['input_position'], tuple(r['edge']), r['coin']): r for r in columns}
    require(len(columns) == len(table) == 18, 'complete measured source columns')
    counts = [[F(0) for _ in range(3)] for _ in range(3)]
    for position in range(3):
        for edge in ((0, 1), (0, 1), (0, 2), (0, 2), (1, 2)):
            for coin in (0, 1):
                row = table[position, edge, coin]
                vector = [source['sign']*int(j == position) for j in range(3)]
                vector[edge[0]] -= row['current']
                vector[edge[1]] += row['current']
                require(vector == [source['sign']*int(j == row['decoded_output']) for j in range(3)],
                        'calibrated column decoder')
                counts[row['decoded_output']][position] += F(1, 10)
    expected = [[F(3, 5), F(1, 5), F(1, 5)], [F(1, 5), F(7, 10), F(1, 10)],
                [F(1, 5), F(1, 10), F(7, 10)]]
    require(counts == expected, 'measured native mean changed')
    prep = source['preparation']
    cool, contrast = prep['cool'], prep['contrast']
    selector = cool['selector_phase']
    re, im = F(selector['real']), F(selector['imag'])
    require(re*re+im*im == 1, 'selector phase is not unitary')
    a = 5*((5*(re*re-im*im)+6*re-3)**2+(10*re*im+6*im)**2)/512
    require(0 < a <= 1, 'selector active weight')
    finite = [[a*counts[i][j]+(1-a)*int(i == j) for j in range(3)] for i in range(3)]
    require([[F(x) for x in row] for row in cool['counts']] == finite, 'inactive selector outcome omitted')
    require([[F(x) for x in row] for row in contrast['counts']] ==
            [[(F(int(i == j))-finite[i][j])/2 for j in range(3)] for i in range(3)],
            'contrast source block differs')


def exact_bounds(record):
    def point(raw):
        sign, mantissa, exponent = raw
        return (-1 if sign else 1)*F(int(mantissa))*F(2)**exponent
    return point(record['lo_dyadic']), point(record['hi_dyadic'])


def refinement_gate(rows):
    indexed = {(r['k'], r['charge_sign'], r['turns'], r['direction'], r['filter_stride']): r
               for r in rows if r['law'] == 'native'}
    for key, fine in indexed.items():
        if key[-1] != 1:
            continue
        coarse = indexed[(*key[:-1], 2)]
        flo, fhi = exact_bounds(fine['visibility'])
        clo, chi = exact_bounds(coarse['visibility'])
        require(fhi < clo if fine['turns'] == 1 else flo > chi,
                'finer source frames did not improve return magnitude')


def run(ks):
    started = time.perf_counter()
    source_path = HERE/'SOURCE_INPUTS.json'
    data = json.loads(source_path.read_text())
    require(data['root_join_bound'] is True, 'fresh common pair not bound')
    require(set(data['sources']) == {'-1', '1'}, 'signed source domain')
    for sign, source in data['sources'].items():
        require(source['sign'] == int(sign), 'source sign binding')
        require(set(source['catalogs']) == {'native', 'swap', 'swap_fixed_native_calibration'}, 'complete law menu')
        check_counts(source)
    cached_arms = {}
    for sign, source in data['sources'].items():
        for law, catalog in source['catalogs'].items():
            require(len(catalog) == 96, 'complete M32 frame catalog')
            require([record['slot'] for record in catalog] == list(range(96)), 'ordered frame slot coverage')
            for stride, h in ((1, 2), (2, 4)):
                for turns in (1, 2):
                    for direction in (-1, 1):
                        cached_arms[sign, law, stride, turns, direction] = arms(
                            catalog, 32, h, turns, direction, stride, law != 'native')
    rows, preparations = [], []
    for k in ks:
        gram, z = counter(k)
        for sign, source in data['sources'].items():
            position, probability, weights = preparation(source['preparation'], gram)
            require(weights[0] > real('.999999'), 'native prepared mode weight')
            preparations.append({'k': k, 'charge_sign': int(sign), 'success': packed(probability),
                'optical_weights': [packed(x) for x in weights],
                'counter_decode_success_given_preparation': packed(hermitian_trace(position)/probability)})
            for law in source['catalogs']:
                if law == 'swap_fixed_native_calibration' and k not in (0, 6):
                    continue
                for stride, h in ((1, 2), (2, 4)):
                    for turns in (1, 2):
                        for direction in (-1, 1):
                            setting = {'k': k, 'charge_sign': int(sign), 'law': law,
                                'M': 32, 'h': h, 'filter_stride': stride, 'turns': turns, 'direction': direction}
                            result = case(position, gram, z, cached_arms[sign, law, stride, turns, direction],
                                          setting, data['classical_timeout_union_bound_each'])
                            rows.append(result)
        print(json.dumps({'progress_k': k, 'cases': len(rows)}), flush=True)
    refinement_gate(rows)
    require(len(rows) == sum(48 if k in (0, 6) else 32 for k in ks), 'complete signed source case menu')
    return {'schema': 'finite-amplified-common-source-interval-v1', 'passed': True,
        'source_inputs_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
        'source_bindings': data['bindings'], 'k_values': ks, 'case_count': len(rows),
        'preparations': preparations, 'cases': rows, 'seconds': time.perf_counter()-started,
        'interval_dps': iv.dps, 'reader_impurity': '0.01', 'reader_probability_allowance': '0.00002',
        'complete_service_budget': budget(), 'service_error_contract': CONTRACT,
        'one_source_error_envelope': {key: str(value) for key, value in envelope(40, 1).items()},
        'pair_error_envelope': {key: str(value) for key, value in envelope(52, 2).items()},
        'finite_implementation': 'Exact rational matching phases, all inactive weights and conditioning retained; outward arithmetic encloses actual circuit.',
        'pair_scope': 'Full independent product instrument on separately admitted signed sources, bound by the fresh global-column join.',
        'external_G6': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, nargs='+', default=[2])
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--record-to', type=Path)
    group.add_argument('--verify', type=Path)
    args = parser.parse_args()
    require(len(set(args.k)) == len(args.k) and all(0 <= k <= 6 for k in args.k), 'certificate k menu')
    if args.record_to:
        require(not args.record_to.exists(), 'refusing to overwrite evidence')
    result = run(args.k)
    if args.record_to:
        with args.record_to.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    else:
        saved = json.loads((args.verify or HERE/'HEADLINE.json').read_text())
        observed = dict(result)
        saved.pop('seconds'); observed.pop('seconds')
        require(saved == observed, 'saved interval evidence differs')
    print(json.dumps({key: value for key, value in result.items() if key not in ('cases', 'preparations', 'source_bindings')}))


if __name__ == '__main__':
    main()
