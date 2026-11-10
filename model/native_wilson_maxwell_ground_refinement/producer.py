"""Portable uniform interacting-ground and complete-record certificate."""
import argparse,hashlib,json
from pathlib import Path
from . import refinement
from model.particle_common_contract.provenance import import_closure
HERE=Path(__file__).resolve().parent
RER=refinement.RER
RESEARCH=HERE/'research_reference.json'
RESEARCH_SHA='b159c06176137e5f0389bd7befe76a3b037bb1150365f5621a6a325e47d28953'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def certificate():
    value=refinement.scientific();parent,_=refinement.admit()
    cutoff=json.loads((RER/refinement.CUTOFF_PARENT).read_text())
    seeds=set(parent['input_sha256'])|set(cutoff['input_sha256'])|{refinement.PARENT,refinement.CUTOFF_PARENT,
        'model/native_wilson_maxwell_ground/ground.py','model/native_wilson_maxwell_refinement/bounds.py',
        'model/particle_common_contract/provenance.py',str(RESEARCH.relative_to(RER))}
    seeds.update(str(p.relative_to(RER)) for p in HERE.iterdir() if p.suffix in ('.py','.md'))
    value['input_sha256']={p:digest(RER/p) for p in sorted(import_closure(seeds))}
    value['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    if digest(RESEARCH)!=RESEARCH_SHA:raise ValueError('immutable research comparison differs')
    ignore={'input_sha256','owned_source_sha256'}
    science=lambda data:{k:v for k,v in data.items() if k not in ignore}
    if science(value)!=science(json.loads(RESEARCH.read_text())):raise ValueError('public scientific fields differ')
    value['independent_research_comparison']={'reference_sha256':RESEARCH_SHA,
        'all_scientific_fields_equal':True,'ignored_only_provenance_fields':sorted(ignore)}
    return value

def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    args=parser.parse_args();value=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    elif value!=json.loads(path.read_text()):raise ValueError('ground refinement receipt differs')
    print('PASS: uniform fixed-cube interacting ground and complete photon-capacity records')
