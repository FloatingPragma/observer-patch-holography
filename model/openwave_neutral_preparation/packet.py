"""Replay local-product preparation and its complete finite source instrument."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from . import local,filtering,records,replay

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def source_bindings():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update('model/'+name+'/'+file for name in ('openwave_neutral_stationary','openwave_neutral_dispersion')
                 for file in ('PROOF.md','results.json'))
    seeds.add('model/requirements.txt')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def report():
    pins=source_bindings()
    result={'schema':'oph.native_neutral_local_product_preparation.v1','input_sha256':pins,
        'explicit_local_product':local.certificate(),
        'finite_filter_and_polar_service':filtering.certificate(),
        'complete_preparation_records':records.certificate(),
        'native_atomic_polynomial_replay':replay.certificate(),
        'verdict':'PASS_EXPLICIT_LOCAL_PRODUCT_HIGH_SUCCESS_APPROXIMATE_STATIONARY_NEUTRAL_PREPARATION'}
    if pins!=source_bindings():raise ValueError('neutral preparation inputs changed during replay')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();result=report();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=result:raise ValueError('neutral preparation receipt differs from replay')
    print(result['verdict'])
