"""Replay the complete finite ordered-cover protection certificate."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..exchange_sector import statistics_packet
from . import protection, service


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    result = statistics_packet.source_bindings()
    for path in sorted(HERE.iterdir()):
        if path.suffix == '.py' or path.name == 'README.md':
            result[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(result.items()))


_IMPORTED_SOURCE_BINDINGS = source_bindings()


def report():
    before = source_bindings()
    if before != _IMPORTED_SOURCE_BINDINGS:
        raise RuntimeError('protection sources changed since import; replay in a fresh interpreter')
    result = {'schema': 'oph.native_ordered_cover_protection.v1',
              'domain': 'binary native H66 tensor one coherent orientation bit; complete C/E raw records',
              'repair_law': 'native endpoint update, followed by record-selected transported-order XOR',
              'memory_law': 'orientation bit retained coherently; only X-commuting additional services admitted',
              'charge': 'I66 tensor X',
              'algebra': {law: protection.algebra_certificate(law) for law in ('C', 'E')},
              'memory_and_localization': protection.memory_and_localization_certificate(),
              'execution': service.certificate(),
              'theorem_boundary': {
                  'exact_discrete_algebra_and_charge_proof': True,
                  'new_Lean_proof_claimed': False,
                  'existing_topology_source': 'model/exchange_sector/STATISTICS.md',
                  'existing_pi1_Z2_uses_declared_triangle_and_commuting_square_fillings': True,
                  'branch_invariance_requires_no_IID_or_all_seam_assumption': True,
                  'full_algebra_and_fixed_space_classification_uses_all_seams': True,
                  'uniform_channel_fixed_space_uses_faithful_uniform_reference': True,
                  'coherence_between_charge_sectors_preserved': False,
                  'localized_source_defect_or_physical_particle_identified': False,
                  'collision_sector_or_refinement_covered': False},
              'source_sha256': before}
    if source_bindings() != before:
        raise RuntimeError('source files changed during protection replay')
    return result


def compare(expected, actual):
    encode = lambda x: json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False)
    if encode(expected) != encode(actual):
        raise ValueError('ordered-cover protection receipt differs from complete replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = report()
    if args.write:
        OUTPUT.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    elif args.check:
        compare(strict_json.load(OUTPUT), value)
    print('Verified native ordered-cover charge, complete instruments and explicit memory boundary.')


if __name__ == '__main__':
    main()
