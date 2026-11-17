"""Exact scope correction; no source regeneration or floating comparison."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'2026-11-16_codex_amplified_source_interval_v1'


def require(value, why):
    if not value:
        raise ValueError(why)


def upper(record):
    sign, mantissa, exponent = record['hi_dyadic']
    return (-1 if sign else 1)*F(int(mantissa))*F(2)**exponent


def run():
    path = SOURCE/'CERTIFICATE.json'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    require(digest == '96e52814f637bfdbeb596c01dd345cefd0153af2073d39ac2e5bd982c22c8872', 'source certificate changed')
    data = json.loads(path.read_text())
    require(data['passed'] is True and len(data['cases']) == 256, 'source case coverage')
    native = [r for r in data['cases'] if r['law'] == 'native']
    control = [r for r in data['cases'] if r['law'] != 'native']
    require(len(native) == 112 and len(control) == 144, 'law coverage')
    n = max(upper(r['pair_off_band']) for r in native)
    c = max(upper(r['pair_off_band']) for r in control)
    require(n < F('0.005486130010'), 'native impurity summary bound')
    require(max(n, c) < F('0.005486165872'), 'all-case impurity summary bound')
    require(c > F('0.005486130010'), 'reported scope distinction not present')
    require(max(n, c) < F('0.01'), 'registered pair gate changed')
    return {'passed': True, 'source_certificate_sha256': digest,
        'native_cases': len(native), 'control_cases': len(control),
        'maximum_native_upper_exact': str(n), 'maximum_control_upper_exact': str(c),
        'native_display_upper': '0.005486130010', 'all_cases_display_upper': '0.005486165872',
        'registered_gate': '0.01', 'source_data_or_gates_changed': False,
        'correction': 'The tight 0.005486130010 bound in the sealed PROOF and FINAL_RECEIPT summary is native-only. The valid all-case bound is 0.005486165872.',
        'external_G6': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    result = run()
    if args.record_to:
        with args.record_to.open('x') as stream:
            stream.write(json.dumps(result, indent=2)+'\n')
    else:
        saved = HERE/'RESULT.json'
        require(saved.is_file(), 'missing saved scope result')
        require(result == json.loads(saved.read_text()), 'scope correction differs')
    print(json.dumps(result, indent=2))
