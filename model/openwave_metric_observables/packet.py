"""Reproduce finite metric, separated-clock and scalar-wave observations."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from . import geometry,clocks,waves

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def source_bindings():
    seeds=[str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md','.lean')]
    seeds+=['model/openwave_native_force/PROOF.md','model/openwave_native_force/results.json',
            'model/finite_area_response/README.md','model/matched_einstein_limit/PROOF.md',
            'model/maxwell_local_selection/PROOF.md']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def report():
    pins=source_bindings()
    value={'schema':'oph.finite_reciprocal_metric_observables.v1','input_sha256':pins,
           'metric':geometry.metric_certificate(),'native_clocks':clocks.certificate(),
           'spatial_scalar_wave':waves.certificate(),'Einstein_premise_test':geometry.einstein_test(),
           'verdict':'PASS_FINITE_METRIC_RESPONSE_WITH_CLOCK_AND_EINSTEIN_OBSTRUCTIONS'}
    if pins!=source_bindings():raise ValueError('metric source changed during replay')
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();value=report();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=value:raise ValueError('metric observer receipt differs from replay')
    print(value['verdict'])
