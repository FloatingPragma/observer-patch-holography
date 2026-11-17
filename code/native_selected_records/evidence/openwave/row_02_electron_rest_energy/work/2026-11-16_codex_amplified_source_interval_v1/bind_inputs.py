"""Bind immutable compiler headers to two actually continued source backgrounds."""
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OPENWAVE = HERE.parents[2]
SPIN = OPENWAVE/'row_07_double_cover/work/2026-11-16_codex_stationary_readout_repair_v2'
PAIR = HERE.parent/'2026-11-16_codex_amplified_common_pair_v1'


def require(value, why):
    if not value:
        raise ValueError(why)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with gzip.open(path, 'rt') as stream:
        return json.load(stream)


def descriptor(record):
    keys = ('counts', 'phase_record', 'preboost', 'selector_phase', 'selector_active_weight',
            'slot', 'M', 'coordinate_system', 'charge_sign', 'source_kind')
    return {key: record[key] for key in keys if key in record}


def build():
    freeze = json.loads((SPIN/'PREREG_V5.json').read_text())
    for name, expected in freeze['files'].items():
        require(sha(SPIN/name) == expected, 'amplified source changed after freeze')
    prep_path = SPIN/'CONFIRMATION/PREPARATION.json.gz'
    cat_path = SPIN/'CONFIRMATION/CATALOG.json.gz'
    cal_path = SPIN/'CONFIRMATION/CALIBRATION_INPUT.json.gz'
    require(sha(prep_path) == 'ef449037d53703c0d8f183edca2b72fd7ad227c203aaaeb1a1ceaa2187c35e6c', 'preparation header differs')
    require(sha(cat_path) == '813b8c2aa45b1f08b2010f53e7ad7f705635059560ba9eb70cfd5a50dfc002ab', 'frame catalog differs')
    prep, catalog, calibration = read(prep_path), read(cat_path), read(cal_path)
    index_path = SPIN/'CONFIRMATION/INDEX.json'
    index = json.loads(index_path.read_text())
    require(len({r['path'] for r in index['files']}) == len(index['files']), 'duplicate final source artifact')
    for row in index['files']:
        require(Path(row['path']).name == row['path'] and sha(index_path.parent/row['path']) == row['sha256'],
                'final source artifact digest')
    summary_path = SPIN/'CONFIRMATION/SUMMARY.json.gz'
    summary = read(summary_path)
    require(summary['gates_pass'] is True and len(summary['rows']) == 256, 'full source confirmation absent')
    pair = read(PAIR/'SOURCE_JOIN.json.gz')
    require(pair['passed'] is True and pair['pair_count'] == 16 and len(pair['rows']) == 16,
            'actual common source join absent')
    tables = {(r['sign'], r['law']): r['columns'] for r in calibration['entries']}
    source_count = 0
    for row in pair['rows']:
        require(len(row['joint_columns']) == 2178, 'common joint columns absent')
        for side, sign in enumerate(row['signs']):
            actual = row['observer_records'][side]
            require(actual['outcome'] == 'ready', 'common source not ready')
            require(actual['calibration']['columns'] == tables[sign, 'native'],
                    'second actual source differs from frozen compiler input')
            require(actual['calibration']['return_position'] == 1 and
                    actual['calibration']['sign'] == sign, 'actual coherent preparation input differs')
            require(actual['blank']['outcome'] == 'blank_checked', 'actual source blank unchecked')
            source_count += 1
    require(source_count == 32, 'actual source coverage')
    paths = [SPIN/'PREREG_V5.json', prep_path, cat_path, cal_path, index_path, summary_path, PAIR/'SOURCE_JOIN.json.gz',
             PAIR/'PREREG.json', PAIR/'AMENDMENT_PREREG.json', PAIR/'TIMING_CLARIFICATION.json']
    bindings = {str(path.relative_to(OPENWAVE)): sha(path) for path in paths}
    return {'candidate': 'A-native-amplified-common-pair-v1', 'root_join_bound': True,
        'root_join_binding_scope': 'Complete primitive tensor map and equality of actual local compiler inputs; complete quantum instruments inherit the frozen deterministic compiler. Numerical full-case comparison is separate.',
        'source_count': source_count, 'bindings': bindings,
        'nominal_source_cases': summary['rows'],
        'classical_timeout_union_bound_each': pair['classical_timeout_union_bound_each'],
        'sources': {str(sign): {'sign': sign, 'native_columns': tables[sign, 'native'],
            'preparation': {name: descriptor(prep[str(sign)]['native'][name]) for name in ('contrast', 'cool')},
            'catalogs': {law: [descriptor(rec) for rec in frames] for law, frames in catalog[str(sign)].items()}}
            for sign in (-1, 1)},
        'pair_hypotheses': 'Separate initial source/ancilla registers; local complete instruments on disjoint ports; common product clock tick; untouched connecting edge; every joint source column checked in both orders.',
        'unregistered_pilot_retained': 'Common classical join first run; explicitly preserved and followed by registered read-only verification.',
        'external_G6': False}


if __name__ == '__main__':
    result = build()
    path = HERE/'SOURCE_INPUTS.json'
    with path.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'bound_sources': result['source_count'], 'source_inputs_sha256': sha(path)}))
