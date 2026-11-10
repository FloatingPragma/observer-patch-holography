#!/usr/bin/env python3
"""Replay five owners against their retained a036 typed receipt pointers.

The predecessor owner authority and a036 receipt bytes remain historical inputs. This
separate authority admits only the explicitly listed receipt hashes; all scientific values, source-map keys, rendered
claims and producer implementations remain unchanged. No generic hash masking
or new source-scope admission is provided.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import strict_json  # noqa: E402
from tools.renew_canonical_live_receipts import (canonical, carry_metadata,
    metadata_view, verify_executable_source_pins)  # noqa: E402
from tools.replay_model_receipts import bindings  # noqa: E402

AUTHORITY = 'tools/field_current_pointer_owner_authority.json'
AUTHORITY_SHA256 = '9239e2d87ce417906cea9014ecc438c9fee7e6f38671b2c3001bd7090043dbf3'
ORDER = ('particle_comparison', 'standard_model_source_map', 'charged_lepton_response', 'native_local_energy_clock', 'native_fair_ground_readout')
MODULES = ('model.particle_comparison.compare', 'model.standard_model_source_map.packet', 'model.charged_lepton_response.packet', 'model.completion.field_current_measurement', 'model.completion.field_current_measurement')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authority():
    if digest(ROOT / AUTHORITY) != AUTHORITY_SHA256:
        raise ValueError('atlas custody authority changed')
    value = strict_json.load(ROOT / AUTHORITY)
    parent=value['predecessor_authority']
    if digest(ROOT/parent['path']) != parent['sha256']:
        raise ValueError('retained atlas owner authority changed')
    if tuple(value['ordered_owners']) != ORDER or set(value['owners']) != set(ORDER):
        raise ValueError('atlas owner inventory changed')
    if tuple(value['owners'][name]['module'] for name in ORDER) != MODULES:
        raise ValueError('atlas producer inventory changed')
    verify_executable_source_pins(value['reviewed_executable_source_sha256'], ROOT)
    return value


def conserve(spec, retained, fresh, *, allow_retained=False, root=ROOT):
    pointers = spec['reviewed_nested_receipt_pointers']
    projected = carry_metadata(spec['receipt'], retained, fresh, pointers=pointers,
                               allow_retained_pointers=allow_retained, scope_append=spec.get('reviewed_source_scope_append'), root=root)
    if canonical(metadata_view(spec['receipt'], projected, pointers, root)) != canonical(
            metadata_view(spec['receipt'], fresh, pointers, root)):
        raise ValueError('atlas scientific payload changed')
    if canonical(projected) != canonical(fresh):
        raise ValueError('atlas replay differs outside reviewed custody')
    return projected


def verify_live_sources(receipt, fresh):
    for _, name, expected, path in bindings(ROOT, receipt, fresh):
        if digest(path) != expected:
            raise ValueError('atlas replay source is not current: ' + name)


def unchanged_renderings(owner, producer, fresh):
    if owner == 'standard_model_source_map':
        expected = {producer.MARKDOWN: producer.render(fresh),
                    producer.MASS_MARKDOWN: producer.render_masses(fresh)}
    elif owner == 'charged_lepton_response':
        expected = {producer.MARKDOWN: producer.render(fresh),
                    producer.CANDIDATES: producer.canonical(fresh['certificates']['candidate_particle_map'])}
    elif owner == 'candidate_correspondence':
        expected = {producer.HERE / 'CANDIDATES.md': producer.markdown(fresh)}
    else:
        expected = {}
    for path, text in expected.items():
        if path.read_text() != text:
            raise ValueError('atlas rendered scientific output changed: ' + str(path))
    return {str(path.relative_to(ROOT)): digest(path) for path in expected}


def execute(owner, *, renew=False):
    config = authority()
    if owner not in ORDER:
        raise ValueError('unknown atlas owner')
    spec = config['owners'][owner]
    retained_path = ROOT / spec['retained_path']
    if digest(retained_path) != spec['retained_sha256']:
        raise ValueError('retained atlas receipt changed')
    if digest(ROOT / spec['producer_source']) != spec['producer_sha256']:
        raise ValueError('atlas scientific producer changed')
    retained = strict_json.load(retained_path)
    output = ROOT / spec['receipt']
    before = output.read_bytes()
    current = strict_json.loads(before.decode())
    checkpoints = spec.get('executed_intermediate_receipts', [])
    for checkpoint in checkpoints:
        if digest(ROOT / checkpoint['path']) != checkpoint['sha256']:
            raise ValueError('executed intermediate atlas receipt changed')
    if hashlib.sha256(before).hexdigest() in {row['sha256'] for row in checkpoints}:
        # Exact bytes of an already executed owner output can retain a superseded
        # live pointer. This admission is enumerated, not an arbitrary hash alias.
        pointers = spec['reviewed_nested_receipt_pointers']
        if canonical(metadata_view(spec['receipt'], current, pointers)) != canonical(
                metadata_view(spec['receipt'], retained, pointers)):
            raise ValueError('executed intermediate atlas science changed')
    else:
        conserve(spec, retained, current, allow_retained=True)
    producer = importlib.import_module(spec['module'])
    fresh = getattr(producer, spec['entrypoint'])()  # Actual unchanged finite producer.
    if hasattr(producer, 'validate'):
        producer.validate(fresh)
    admitted = conserve(spec, retained, fresh)
    verify_live_sources(spec['receipt'], admitted)
    rendered = unchanged_renderings(owner, producer, fresh)
    authority()  # No executable may change during the replay.
    if output.read_bytes() != before:
        raise ValueError('canonical atlas output changed during replay')
    if renew:
        # Match the owner's actual serialization. All original science is exact.
        text = (producer.canonical(admitted) if hasattr(producer, 'canonical') else
                json.dumps(admitted, indent=2, sort_keys=True, ensure_ascii=False) + '\n')
        output.write_text(text)
    elif canonical(current) != canonical(admitted):
        raise ValueError('canonical atlas custody is stale; run --renew')
    return {'owner': owner, 'receipt': spec['receipt'],
            'retained_path': spec['retained_path'], 'retained_sha256': spec['retained_sha256'],
            'current_sha256': digest(output),
            'unchanged_producer': spec['producer_source'],
            'actual_producer_replayed': True, 'scientific_payload_exactly_preserved': True,
            'source_scope_exactly_preserved': not bool(spec.get('reviewed_source_scope_append')), 'physical_statuses_changed': False,
            'reviewed_source_scope_append': spec.get('reviewed_source_scope_append', {}),
            'typed_nested_receipt_pointers': spec['reviewed_nested_receipt_pointers'],
            'unchanged_rendered_outputs': rendered}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--owner', choices=ORDER, action='append')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--renew', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    selected = set(args.owner or ORDER)
    reports = [execute(name, renew=args.renew) for name in ORDER if name in selected]
    print(json.dumps(reports, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
