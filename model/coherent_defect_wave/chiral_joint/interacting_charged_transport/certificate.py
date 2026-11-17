"""Custody of same-H charge motion, preparation and all route records."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .algebra import certificate as algebra_certificate
from .transport import certificate as transport_certificate

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
DEPENDENCIES=tuple('model/coherent_defect_wave/chiral_joint/'+name+'/'+file
    for name in ('interacting_charged_response','interacting_ground_preparation',
                 'collective_spin_response','atomic_reference_control','atomic_spin_geometry','local_reference')
    for file in ('PROOF.md','results.json'))


def input_hashes():
    seeds={str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(DEPENDENCIES);seeds.add('model/requirements.txt')
    for name in seeds:
        if not (ROOT/name).is_file():raise FileNotFoundError('missing declared charged transport input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def certificate():
    return {'schema':'oph.same_H_interacting_charged_transport.v1',
        'actual_route_and_charge_algebra':algebra_certificate(),
        'complete_physical_channels':{family:transport_certificate(family) for family in ('ec','L')},
        'input_sha256':input_hashes()}
