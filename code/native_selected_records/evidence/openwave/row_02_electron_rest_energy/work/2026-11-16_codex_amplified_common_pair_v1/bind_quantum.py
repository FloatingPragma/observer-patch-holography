"""Read-only binding of actual pair entry to complete finite quantum records.

This checks the common operator parameters and reviewed state/cost evidence;
it does not count those matrix calculations as additional source histories.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OPENWAVE = HERE.parents[2]
INTERVAL = HERE.parent/'2026-11-16_codex_amplified_source_interval_v1'
REVIEW = OPENWAVE/'row_07_double_cover/reviews/2026-11-16_codex_amplified_source_adversarial'


def require(value, why):
    if not value:
        raise ValueError(why)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(folder):
    lines = (folder/'MANIFEST.sha256').read_text().splitlines()
    for line in lines:
        expected, name = line.split('  ', 1)
        path = (folder/name).resolve()
        require(path.is_relative_to(folder.resolve()) and sha(path) == expected, 'changed bound source leaf')
    return {'manifest_sha256': sha(folder/'MANIFEST.sha256'), 'leaves': len(lines)}


def run():
    folders = [INTERVAL, REVIEW, HERE.parent/'2026-11-16_codex_record_calibration_v1',
        OPENWAVE/'row_02_electron_rest_energy/reviews/2026-11-16_codex_amplified_interval_adversarial']
    providers = {str(folder.relative_to(OPENWAVE)): verify_manifest(folder) for folder in folders}
    inputs = json.loads((INTERVAL/'SOURCE_INPUTS.json').read_text())
    require(inputs['bindings'][str((HERE/'SOURCE_JOIN.json.gz').relative_to(OPENWAVE))] == sha(HERE/'SOURCE_JOIN.json.gz'),
            'interval did not consume this actual source join')
    with gzip.open(HERE/'SOURCE_JOIN.json.gz', 'rt') as stream:
        source = json.load(stream)
    count = 0
    for row in source['rows']:
        for side, sign in enumerate(row['signs']):
            observed = row['observer_records'][side]
            require(observed['calibration']['columns'] == inputs['sources'][str(sign)]['native_columns'],
                    'actual continued source and compiler table differ')
            require(observed['outcome'] == 'ready' and observed['calibration']['return_position'] == 1,
                    'actual finite preparation entry differs')
            count += 1
    require(count == 32, 'two-source input coverage')
    left, right = inputs['sources']['-1'], inputs['sources']['1']
    records = [(left['preparation'][key], right['preparation'][key]) for key in ('contrast', 'cool')]
    records += [(a, b) for law in left['catalogs'] for a, b in zip(left['catalogs'][law], right['catalogs'][law])]
    for a, b in records:
        require(a['counts'] == b['counts'] and a['preboost'] == b['preboost'] and
                a['phase_record']['real'] == b['phase_record']['real'] and
                a['phase_record']['imag'] == b['phase_record']['imag'], 'signed optical instruments differ')
    audit = json.loads((REVIEW/'WORKSPACE_BRIDGE.json').read_text())
    bridge = audit['signed_bridge']
    require(audit['passed'] is True and bridge['clock_refinements'] == 7 and
            bridge['corresponding_case_pairs'] == 128 and bridge['returned_density_pairs'] == 256,
            'whole signed-clock covariance review incomplete')
    interval = json.loads((INTERVAL/'CERTIFICATE.json').read_text())
    require(interval['passed'] is True and interval['case_count'] == 256, 'joint finite interval evidence incomplete')
    return {'passed': True, 'actual_continued_sources': count,
        'signed_optical_parameter_pairs': len(records), 'all_four_charge_sign_products': True,
        'whole_signed_clock_bridge': bridge, 'providers': providers,
        'interval_certificate_sha256': sha(INTERVAL/'CERTIFICATE.json'),
        'actual_source_sha256': sha(HERE/'SOURCE_JOIN.json.gz'),
        'scope': 'Actual common entry, complete primitive tensor operator, identical full local compiler inputs, reviewed signed state covariance and fresh finite pair interval.',
        'physical_claim_proof': 'PROOF.md', 'external_G6': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    if args.record_to:
        require(not args.record_to.exists(), 'refusing to overwrite quantum binding')
    result = run()
    if args.record_to:
        with args.record_to.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    else:
        require(result == json.loads((HERE/'QUANTUM_BINDING.json').read_text()), 'saved quantum binding differs')
    print(json.dumps(result))
