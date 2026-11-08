"""Check four retained ordered-spin receipts against their complete producers."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'results.json'
PROVIDERS = (
    'spin_particle_scaling',
    'spin_canonical_field',
    'spin_canonical_energy',
    'spin_current_spectrum',
)


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n'


def source_bindings():
    paths = {'model/__init__.py', 'tools/strict_json.py'}
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.suffix in ('.py', '.md'))
    for owner in PROVIDERS:
        provider = importlib.import_module(f'model.{owner}.packet')
        pins = (provider.source_bindings() if hasattr(provider, 'source_bindings')
                else provider.bindings())
        paths.update(pins)
        paths.add(f'model/{owner}/results.json')
    # Package initialization is executed when the declared providers and their
    # local dependencies are imported. Bind those files even when a provider
    # lists only its mathematical inputs.
    for package in ('model', 'tools', 'model/spin_goldstone_sector',
                    *(f'model/{owner}' for owner in PROVIDERS)):
        init = ROOT / package / '__init__.py'
        if init.is_file():
            paths.add(str(init.relative_to(ROOT)))
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def build():
    before = source_bindings()
    providers = {}
    for owner in PROVIDERS:
        module = importlib.import_module(f'model.{owner}.packet')
        path = ROOT / 'model' / owner / 'results.json'
        retained = strict_json.load(path)
        # validate compares the entire retained JSON to a fresh complete
        # producer result, including types, scientific flags and source pins.
        module.validate(retained)
        providers[owner] = {
            'receipt': str(path.relative_to(ROOT)),
            'schema': retained['schema'],
            'retained_equals_complete_replay': True,
            'theorem': retained['theorem'],
            'inherited_custody': retained['inherited_custody'],
            'provider_status_flags': {
                key: retained[key] for key in
                ('new_checks_pass', 'dependency_custody_pass', 'all_pass')
                if key in retained
            },
        }
    result = {
        'schema': 'oph.spin_source_chain.v1',
        'source_bindings': before,
        'providers': providers,
        'verification_checks_pass': True,
        'scope': {
            'kind': 'Retained-receipt replay adapter; no additional physical theorem',
            'complete_replays': list(PROVIDERS),
            'all_ancestral_producers_replayed': False,
            'finite_arithmetic_alone_proves_continuum_theorems': False,
            'source_geometry_or_selected_state_removed_as_inputs': False,
            'physical_particle_species_or_masses_identified': False,
        },
    }
    if source_bindings() != before:
        raise RuntimeError('source changed during complete retained-receipt replay')
    return result


def validate(value):
    if type(value) is not dict or value.get('schema') != 'oph.spin_source_chain.v1':
        raise ValueError('invalid spin source chain receipt')
    if canonical(value) != canonical(build()):
        raise ValueError('retained receipt, source, provider result or scope changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    mode.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    if args.check is not None:
        validate(strict_json.load(args.check))
    else:
        args.write.write_text(canonical(build()))
    print('spin_source_chain: retained receipts reproduce exactly')


if __name__ == '__main__':
    main()
