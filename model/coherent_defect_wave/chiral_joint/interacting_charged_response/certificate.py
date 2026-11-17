"""Source custody for odd, compensated records on the actual interacting ground."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import operator_certificate
from .local import certificate as local_certificate
from .transfer import certificate as transfer_certificate

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
DEPENDENCIES=tuple('model/coherent_defect_wave/chiral_joint/'+name+'/'+file
    for name in ('charged_ground_exchange','atomic_reference_control','local_reference',
                 'collective_spin_response','interacting_ground_preparation')
    for file in ('PROOF.md','results.json'))


def input_hashes():
    seeds={str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(DEPENDENCIES);seeds.add('model/requirements.txt')
    for name in seeds:
        if not (ROOT/name).is_file():raise FileNotFoundError('missing declared interacting charged input: '+name)
    closure=import_closure(seeds)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(closure)}


def certificate():
    return {'schema':'oph.same_H_interacting_odd_charged_response.v1',
        'complete_probe_and_reference':operator_certificate(),
        'exact_local_channels':{family:local_certificate(family) for family in ('ec','L')},
        'same_H_interacting_channels':{family:transfer_certificate(family) for family in ('ec','L')},
        'input_sha256':input_hashes()}
