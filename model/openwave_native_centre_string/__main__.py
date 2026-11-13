"""Produce and cold-replay the canonical quotient-Gauss string."""
import argparse
import hashlib
import json
from pathlib import Path
from model.openwave_strong_step_scaling import source as parent_source
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .centre import certificate

HERE = Path(__file__).resolve().parent
ROOT = Path(parent_source.__file__).resolve().parents[2]
PREFIX = 'model/openwave_native_centre_string'


def pins():
    assert HERE == ROOT / PREFIX, 'Canonical admission requires the declared repository source path.'
    files = set()
    for name in ('openwave_strong_step_scaling', 'openwave_quark_flux_records'):
        parent = f'model/{name}/results.json'
        prior = json.loads((ROOT / parent).read_text())
        for p, digest in prior['source_sha256'].items():
            assert hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == digest, p
        files |= set(prior['source_sha256']) | {parent}
    files |= {str(p.relative_to(ROOT)) for p in HERE.iterdir()
              if p.suffix in ('.py', '.md', '.tex')}
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
            for p in sorted(import_closure(files))}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    before = pins(); result = certificate(); assert before == pins()
    result['source_sha256'] = before
    path = HERE / 'results.json'
    if args.write:
        path.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    else:
        assert equivalent(json.loads(path.read_text()), result)
    print(result['verdict'])


if __name__ == '__main__':
    main()
