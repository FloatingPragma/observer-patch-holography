"""Reproduce source-bound neutral threshold and boundary comparisons."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from . import source,threshold

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def source_bindings():
    seeds=[str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md','.lean')]
    seeds+=['model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/results.json',
            'model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/PROOF.md',
            'model/coherent_defect_wave/chiral_joint/atomic_even_remainder/PROOF.md',
            'model/coherent_defect_wave/chiral_joint/ground_composite_response/PROOF.md',
            'model/coherent_defect_wave/chiral_joint/peter_weyl_links/PROOF.md',
            'model/dark_spectral_source/README.md']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def report():
    pins=source_bindings()
    result={'schema':'oph.native_neutral_threshold.v1','input_sha256':pins,
            'currents':source.current_certificate(),'boundary':source.boundary_certificate(),
            'threshold':threshold.certificate(),'dark_source_control':source.dark_number_control(),
            'verdict':'PASS_NATIVE_NEUTRAL_GAP_AND_SEPARATE_PERIODIC_ZERO_MOMENTUM_CANDIDATE'}
    if pins!=source_bindings():raise ValueError('neutral source changed during replay')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();value=report();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=value:raise ValueError('neutral threshold differs from replay')
    print(value['verdict'])
