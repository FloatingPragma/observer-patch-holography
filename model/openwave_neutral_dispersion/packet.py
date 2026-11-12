"""Replay the exact neutral flatness theorem on its complete source."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from . import symmetry,flatness

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def source_bindings():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update('model/'+name+'/'+file for name in ('openwave_neutral_stationary','openwave_neutral_candidate')
                 for file in ('PROOF.md','results.json'))
    seeds.add('model/requirements.txt')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def report():
    pins=source_bindings()
    result={'schema':'oph.native_neutral_exact_flatness.v1','input_sha256':pins,
        'complete_carrier_symmetry':symmetry.certificate(),
        'neutral_band_flatness_and_records':flatness.certificate(),
        'verdict':'PASS_EXACT_PERIODIC_SIDE_TWO_ODD_SECTOR_MOMENTUM_FLATNESS'}
    if pins!=source_bindings():raise ValueError('neutral symmetry inputs changed during replay')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();result=report();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=result:raise ValueError('neutral flatness receipt differs from replay')
    print(result['verdict'])
