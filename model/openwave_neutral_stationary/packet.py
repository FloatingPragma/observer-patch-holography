"""Replay the neutral source, spectral preparation and complete readouts."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from . import source,isolation,ground,translations,local_profile,filter,records,dynamics,replay

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DEPENDENCIES=('local_reference','atomic_reference_control','atomic_even_remainder',
              'atomic_spin_geometry','ground_composite_response','peter_weyl_links')


def source_bindings():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update('model/coherent_defect_wave/chiral_joint/'+name+'/'+file
                 for name in DEPENDENCIES for file in ('PROOF.md','results.json'))
    seeds.update(('model/openwave_neutral_candidate/PROOF.md',
                  'model/openwave_neutral_candidate/results.json','model/requirements.txt'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def report():
    pins=source_bindings()
    result={'schema':'oph.native_neutral_stationary.v1','input_sha256':pins,
        'source':source.certificate(),'isolation':isolation.certificate(),
        'ground':ground.certificate(),'translations':translations.certificate(),
        'local_profile':local_profile.certificate(),'filter':filter.certificate(),
        'complete_stationary_records':records.certificate(),'dynamics':dynamics.certificate(),
        'native_atomic_polynomial_replay':replay.certificate(),
        'verdict':'PASS_DECLARED_PERIODIC_KAPPA_TWO_NEUTRAL_STATIONARY_BAND_AND_COMPLETE_RECORDS'}
    if pins!=source_bindings():raise ValueError('neutral source changed during replay')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();result=report();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=result:raise ValueError('neutral stationary receipt differs from replay')
    print(result['verdict'])
