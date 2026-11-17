"""Generate the official executable-model claim projection."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools.build_ccf_lt2_claims import generate as official_claims


def generate():
    official = official_claims()
    return {
        'schema': 'oph.native_execution_claims.v2',
        'official_model': 'CCF-LT-2',
        'specification': 'model/SPEC.md',
        'runtime': 'model/ccf_lt2',
        'mission': 'One specified executable construction reproduces compatible observer records for physical observables.',
        'architecture_uniqueness_required': False,
        'finite_construction': official,
        'physical_premises_closed': [],
        'independent_physical_masses_GeV_added': [],
        'complete_physical_mission_proved': False,
        'acceptance': 'Calculate complete propagated master-source histories and instantiate physical spectrum, common geometry, clock and matched Einstein premises on that source.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    payload=json.dumps(generate(),sort_keys=True,indent=2)+'\n'
    out=ROOT/'claims/native_execution.json'
    if args.write:
        out.write_text(payload)
    elif out.read_text()!=payload:
        raise SystemExit('Official native execution register differs')
    print('PASS_CCF_LT2_NATIVE_EXECUTION_CLAIMS')


if __name__ == '__main__':
    main()
