"""Bind the signed-source charge, rest-energy and frame results to their evidence."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_signed_source'
OUTPUT = ROOT / 'claims/openwave_signed_source.json'

ROWS = [
    {'row': 1, 'criterion': 'Charge quantization',
     'result': ('A primitive integer source on the twelve-port cell: the integer face certificate '
                'ker B = im D excludes two sharp half-charge daughters of one signed pair, and the '
                'law-dependent stopped flow E[F] = B^T L^+ q is measured: 264 native trials absorb, '
                '48 changed-law controls do not, a planted half-current is detected.'),
     'scope': ('A probabilistic alternative between two integer outcomes is distinguished from a '
               'fractional charge. The charge dictionary reads the primitive boundary source as '
               'electric charge.')},
    {'row': 2, 'criterion': 'Electron rest energy (mass)',
     'result': ('A compact four-port charged ray has positive energy E6 = -1260 log(209/210) = '
                '6.014331228348673 from the full 630-seam transfer operator of the same source; '
                'the disjoint pair has E12 = 12.057508280349891, which is 1024.4508487 keV in the '
                'unit fixed by 511 keV.'),
     'scope': ('511 keV is the permitted calibration of the energy unit, not a predicted number. '
               'Lower nonconstant modes and the complementary energy-pointer outcomes are retained.')},
    {'row': 7, 'criterion': 'Spin-1/2 statistics (720 degree double cover)',
     'result': ('Independently measured first-current mean m = ((phi-3)/15) v_p and covariance '
                'C_parallel = (3-phi)/9, C_transverse = (phi+2)/15 over 720 branches: the executed '
                'proper half-turn returns the covariance and reverses the nonzero mean; the second '
                'half-turn returns both. Fixed-seam and off-axis controls discriminate.'),
     'scope': ('The executed rotation is a finite symmetry of the addressed cell, distinct from a '
               'continuous interpolation; no exchange-statistics claim is made.')},
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    spec = importlib.util.spec_from_file_location('signed_source_verify', PACKAGE / 'verify.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    verification = module.verify()
    return {
        'schema': 'oph.openwave.signed_source.v1',
        'construction': ('One signed source on the twelve-port carrier under the native repair law, '
                         'with a complete quantum instrument and current pointer; the same work '
                         'register is prepared, evolved and read.'),
        'declared_architecture': ['positive-amplitude memory lift of the repair branches',
                                  'finite carrier graph (240 ports, 630 seams; twelve-port cell)',
                                  'preparation and pointer services',
                                  'physical dictionary: boundary source as charge, transfer '
                                  'spectrum as energy, first-current moments as frame and field'],
        'rows': [dict(r, criterion_met_on_construction=True, official_status='untested')
                 for r in ROWS],
        'evidence': {'proof': 'code/openwave_signed_source/PROOF.md',
                     'recompute': 'python3 code/openwave_signed_source/verify.py',
                     'full_replay': verification['full_replay'],
                     'published_package': 'code/openwave_static_signed_source',
                     'published_evidence_files_verified': verification['published_evidence_files_verified'],
                     'sealed_native_calls': verification['sealed_native_calls'],
                     'recomputed': verification['recomputed']},
        'official_OpenWave_awards_added': 0,
        'source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in
                          [Path(__file__), PACKAGE / 'verify.py', PACKAGE / 'recompute.py', PACKAGE / 'PROOF.md',
                           ROOT / 'code/openwave_static_signed_source/MANIFEST.json',
                           ROOT / 'code/openwave_static_signed_source/replay.py']},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = json.dumps(generate(), indent=2, sort_keys=True) + '\n'
    if args.write:
        OUTPUT.write_text(text)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != text:
        raise SystemExit('Stale signed-source claim ledger')
    print('PASS_SIGNED_SOURCE_CLAIM')


if __name__ == '__main__':
    main()
