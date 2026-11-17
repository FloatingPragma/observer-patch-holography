from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from . import waves, electrostatic, mass_shell, orbitals

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = (
    'model/native_wilson_maxwell/field.py',
    'model/maxwell_local_selection/PROOF.md',
    'model/native_mass_channels/algebra.py',
    'model/source_mass_shell/PROOF.md',
)


def produce():
    jobs = [waves.run(), electrostatic.run(), mass_shell.run(), orbitals.run()]
    for job in jobs:
        for group in ('checks', 'mutations'):
            if not job[group] or any(v is not True for v in job[group].values()):
                raise ArithmeticError('missing passed check or broken-source control')
    paths = [*SOURCES, *(str(p.relative_to(ROOT)) for p in sorted(HERE.glob('*.py'))),
             'model/openwave_field_observables/PROOF.md']
    return {'schema': 'oph.finite_field_observables.v1',
            'source_basis_commit': 'd88d4854744037accdee6a24cae21b1d7623daf1',
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
            'jobs': jobs, 'verdict': 'PASS_SCOPED_FIELD_OBSERVABLES'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = produce()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    elif json.loads(path.read_text()) != value:
        raise SystemExit('stored scientific output or source custody differs')
    print(value['verdict'])


if __name__ == '__main__':
    main()
