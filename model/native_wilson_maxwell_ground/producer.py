"""Frozen-source ground obstruction, declared local completion and preparation."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
sys.path.insert(0,str(RER))
from model.native_wilson_maxwell_ground import ground
from model.particle_common_contract.provenance import import_closure
RESEARCH=HERE/'research_reference.json'
RESEARCH_SHA='34df20f3caec2161b84c2fcc36c450ba8dbb9eda78f8db2ecaa98db1425ca4c5'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def certificate():
    value=ground.calculate();parent=ground.admit()
    seeds=set(parent['input_sha256'])|{ground.PARENT,ground.GEOMETRY,
        'model/native_wilson_maxwell_refinement/bounds.py',
        'model/native_block_dirac/algebra.py','model/particle_common_contract/provenance.py',
        'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/PROOF.md'}
    seeds.add(str(RESEARCH.relative_to(RER)))
    seeds.update(str(path.relative_to(RER)) for path in HERE.iterdir() if path.suffix in ('.py','.md'))
    value['input_sha256']={name:digest(RER/name) for name in sorted(import_closure(seeds))}
    value['owned_source_sha256']={path.name:digest(path) for path in sorted(HERE.iterdir()) if path.suffix in ('.py','.md')}
    if digest(RESEARCH)!=RESEARCH_SHA:raise ValueError('immutable ground comparison differs')
    ignore={'input_sha256','owned_source_sha256'}
    science=lambda data:{k:v for k,v in data.items() if k not in ignore}
    if science(value)!=science(json.loads(RESEARCH.read_text())):
        raise ValueError('public ground science differs from frozen research')
    value['independent_research_comparison']={'reference_sha256':RESEARCH_SHA,
        'all_scientific_fields_equal':True,'ignored_only_provenance_fields':sorted(ignore)}
    return value


def main():
    parser=argparse.ArgumentParser();m=parser.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true')
    args=parser.parse_args();value=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:
        if json.loads(path.read_text())!=value:raise ValueError('ground certificate differs')
    print('PASS: full-source obstruction, positive macrocell completion and complete ground records')


if __name__=='__main__':main()
