"""Fail-closed custody for the actual local ground composite readout."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from .spectrum import certificate as spectral_certificate
from .apparatus import certificate as apparatus_certificate

ROOT=Path(__file__).resolve().parents[4]
PACKAGE='model/coherent_defect_wave/chiral_joint/ground_composite_response/'


def input_hashes():
    seeds={PACKAGE+name for name in ('__init__.py','__main__.py','model.py','spectrum.py',
        'apparatus.py','certificate.py','test_response.py','PROOF.md','README.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/local_reference/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/local_reference/results.json',
        'model/coherent_defect_wave/chiral_joint/neutral_composite/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/higgs_composites/PROOF.md'))
    for name in sorted(seeds):
        if not (ROOT/name).is_file():raise FileNotFoundError('missing declared ground response input: '+name)
    paths=import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT/name).is_file():raise FileNotFoundError('missing ground response import input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema':'oph.global_physical_local_ground_composite_response.v1',
        'source':{'Hamiltonian':'all native local Yukawa channels plus scalar excitation projector',
            'coefficients':{'up':'1','down':'1','charged_lepton':'1','scalar_penalty':'1'},
            'clock':'exp(-it H_local), hbar=1','full_CAR_modes':30,'scalar_states':5,
            'reference':'normalized global zero-boundary-charge physical ground doublet projector',
            'reference_selection_is_declared_preparation_law':True,
            'operator':'J_alpha=phi_0 c_alpha,11-phi_1 c_alpha,10 and its actual adjoint',
            'spatial_or_composite_hopping_added':False},
        'complete_apparatus':apparatus_certificate(),
        'actual_ground_composite_spectrum':spectral_certificate(),
        'scope':{'actual_interacting_local_ground_response':True,
            'blank_input_composite_preparation_substituted':False,
            'unconstrained_or_boundary_compensated_reference':False,
            'interacting_many_cell_vacuum_or_dispersion_proved':False,
            'observed_relativistic_particle_identification_proved':False},
        'input_sha256':input_hashes()}
