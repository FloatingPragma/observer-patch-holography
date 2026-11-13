"""Produce and cold-replay the canonical native collision source."""
import argparse
import hashlib
import json
from pathlib import Path
from model.openwave_passive_clock_records.clock import ROOT
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .clock import certificate

HERE = Path(__file__).resolve().parent
PREFIX = 'model/openwave_autonomous_clock_channel'


def pins():
    assert HERE == ROOT / PREFIX, 'Canonical admission requires the declared repository source path.'
    parent = 'model/openwave_passive_clock_records/results.json'
    prior = json.loads((ROOT / parent).read_text())
    files = set(prior['source_sha256']) | {parent}
    for p, digest in prior['source_sha256'].items():
        assert hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == digest, p
    files |= {str(p.relative_to(ROOT)) for p in HERE.iterdir()
              if p.suffix in ('.py', '.md', '.tex') or p.name == 'witnesses.json'}
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
            for p in sorted(import_closure(files))}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    before = pins()
    result = certificate()
    assert before == pins()
    result['source_sha256'] = before
    path = HERE / 'results.json'
    if args.write:
        path.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    else:
        assert equivalent(json.loads(path.read_text()), result)
    print(result['verdict'])


if __name__ == '__main__':
    main()
