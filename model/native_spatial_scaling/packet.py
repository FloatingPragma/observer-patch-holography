"""Strict source-bound replay of finite native spatial motion and clocks."""
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
from model.native_program import receipt as native_receipt
from model.native_sector_readout import packet as detector_receipt
from . import geometry,propagation,records

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple('model/native_spatial_scaling/'+name for name in (
    '__init__.py','geometry.py','propagation.py','records.py','packet.py',
    'test_spatial_scaling.py','PROOF.md','README.md'))+(
    'model/federated_observer/SCALING.md','model/federated_observer/scaling.py',
    'model/physical_space_time/OCCUPATION_FIELD.md',
    'model/observer_time_fields/REFINEMENT.md','model/native_local_response/PROOF.md',
    'model/native_local_response/bounds.py','model/native_local_response/cone.py',
    'model/native_sector_readout/readout.py','model/native_sector_readout/PROOF.md',
    'tools/strict_json.py','model/requirements.txt')


def source_bindings():
    paths=set(SOURCE_FILES)|set(geometry.affine.SOURCE_PATHS)
    paths.update(native_receipt.source_bindings());paths.update(detector_receipt.source_bindings())
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)
            if not name.startswith('docs/')}


def produce():
    before=source_bindings()
    value={'schema':'oph.native_spatial_scaling.v1','geometry':geometry.geometry_certificate(),
           'native_displacement':propagation.certificate(),'complete_record_access':records.certificate(),
           'scope':{'stationary_one_unit_preparation_declared':True,
                    'independent_uniform_global_selector_required':True,
                    'physical_spatial_coordinate_or_clock_identified':False,
                    'additional_Poisson_clock_inserted':False,
                    'arbitrary_quantum_input_displacement_bound_claimed':False,
                    'existing_density_continuum_reproved_as_new':False,
                    'all_level_floating_constructor_return_claimed':False,
                    'host_limits_changed':False,'particle_mass_or_relativistic_QFT_claimed':False},
           'source_bindings':before}
    if source_bindings()!=before:raise ArithmeticError('source files changed during spatial replay')
    return value


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('native spatial-scaling receipt differs from source replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:OUTPUT.write_text(canonical(produce()));print('Native spatial-scaling receipt written.')
    else:verify();print('Native spatial-scaling receipt verified.')
