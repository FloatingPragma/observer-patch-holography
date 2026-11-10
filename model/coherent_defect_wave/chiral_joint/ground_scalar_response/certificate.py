"""Fail-closed scalar response custody on the actual physical local ground."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from .model import source_identities,matrix_certificate
from .spectrum import certificate as spectrum_certificate
from .apparatus import certificate as apparatus_certificate

ROOT=Path(__file__).resolve().parents[4]
PACKAGE='model/coherent_defect_wave/chiral_joint/ground_scalar_response/'


def input_hashes():
    seeds={PACKAGE+name for name in ('__init__.py','__main__.py','model.py','spectrum.py',
        'apparatus.py','certificate.py','test_response.py','PROOF.md','README.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/local_reference/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/local_reference/results.json',
        'model/coherent_defect_wave/chiral_joint/ground_composite_response/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/ground_composite_response/results.json'))
    for name in sorted(seeds):
        if not (ROOT/name).is_file():raise FileNotFoundError('missing declared ground scalar input: '+name)
    paths=import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT/name).is_file():raise FileNotFoundError('missing ground scalar import input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema':'oph.actual_local_ground_scalar_response.v1',
        'source':{'Hamiltonian':'all native local Yukawa channels plus scalar excitation projector, all coefficients one',
            'full_CAR_modes':30,'scalar_states':5,'clock':'exp(-it H_local), hbar=1',
            'reference':'P_ground/2 on the full zero-boundary-charge local physical ground',
            'operators':['N_H: actual scalar excitation projector','Y: sum of all three native Yukawa channels'],
            'individual_Yukawa_channel_response_claimed':False},
        'full_source_identities':source_identities(),'exact_energy_and_current_matrices':matrix_certificate(),
        'complete_apparatus':apparatus_certificate(),'scalar_response':spectrum_certificate(),
        'scope':{'finite_clock_scalar_transition':True,'relativistic_rest_mass_identification':False,
            'many_cell_stationary_reference_or_dispersion':False,'arbitrary_input_ground_filter_protocol':False},
        'input_sha256':input_hashes()}
