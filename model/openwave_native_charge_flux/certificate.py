"""Source-bound complete current/flux attachment receipt."""
from pathlib import Path
import hashlib
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from .algebra import endpoint_certificate,gauss_cut_census,route_certificate
from .records import certificate as record_certificate

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DEPENDENCIES=tuple('model/coherent_defect_wave/chiral_joint/'+name+'/'+file
    for name in ('interacting_charged_transport','interacting_charged_response',
                 'interacting_ground_preparation','interacting_charge_localization',
                 'peter_weyl_links','collective_spin_response','atomic_reference_control',
                 'atomic_spin_geometry','local_reference')
    for file in ('PROOF.md','results.json'))


def input_hashes():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(DEPENDENCIES);seeds.update(SOURCE_FILES);seeds.add('model/requirements.txt')
    for name in seeds:
        if not (ROOT/name).is_file():raise FileNotFoundError('missing native charge/flux input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def certificate():
    return {'schema':'oph.native_joint_charge_flux.v1',
        'endpoint_and_full_current_algebra':endpoint_certificate(),
        'actual_Gauss_compiler_census':gauss_cut_census(),
        'actual_route_attachment':route_certificate(),
        'complete_joint_records':record_certificate(),
        'input_sha256':input_hashes(),
        'verdict':'PASS_SAME_SOURCE_UNIT_CHARGE_GAUSS_FLUX_ATTACHMENT'}


def validate_pins(value):
    for name,expected in value['input_sha256'].items():
        path=Path(name)
        if path.is_absolute() or '..' in path.parts or hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
            raise ValueError('native charge/flux input changed: '+name)
