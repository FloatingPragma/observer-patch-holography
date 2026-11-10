"""Fail-closed custody for the finite algebraic preparation schedule."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import certificate as preparation_certificate

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).parent
DEPENDENCIES = tuple('model/coherent_defect_wave/chiral_joint/'+name+'/'+file
    for name in ('collective_spin_response','atomic_spin_geometry',
                 'atomic_reference_control','local_reference','neutral_local_spectrum')
    for file in ('PROOF.md','results.json'))


def input_hashes():
    seeds = {str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(DEPENDENCIES)
    seeds.add('model/coherent_defect_wave/chiral_joint/neutral_local_spectrum/equitable_partition.npz')
    seeds.add('model/requirements.txt')
    for name in seeds:
        if not (ROOT/name).is_file():
            raise FileNotFoundError('missing declared interacting preparation input: '+name)
    paths = import_closure(seeds)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema':'oph.finite_algebraic_interacting_ground_preparation.v1',
        'preparation':preparation_certificate(),'input_sha256':input_hashes()}
