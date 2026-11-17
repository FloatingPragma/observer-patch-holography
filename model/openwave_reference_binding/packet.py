"""Reproduce source-bound finite reference-binding certificates."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from . import binding, strings

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def source_bindings():
    seeds = [str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py', '.md', '.lean')]
    seeds += ['model/coherent_defect_wave/chiral_joint/interacting_charge_localization/PROOF.md',
              'model/coherent_defect_wave/chiral_joint/interacting_charged_response/PROOF.md',
              'model/coherent_defect_wave/chiral_joint/peter_weyl_links/PROOF.md',
              'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json',
              'model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/results.json']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def report():
    pins = source_bindings()
    result = {'schema': 'oph.finite_reference_binding.v1', 'input_sha256': pins,
              'admissible_native_strings': strings.certificate(), **binding.certificate(),
              'verdict': 'PASS_CONSTRAINED_REFERENCE_BINDING_WITH_SCREENED_COMPLEMENT'}
    if pins != source_bindings(): raise ValueError('binding inputs changed during replay')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true'); group.add_argument('--check', action='store_true')
    args = parser.parse_args(); path = HERE/'results.json'; result = report()
    if args.write: path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    elif json.loads(path.read_text()) != result: raise ValueError('reference-binding receipt differs from replay')
    print(result['verdict'])
