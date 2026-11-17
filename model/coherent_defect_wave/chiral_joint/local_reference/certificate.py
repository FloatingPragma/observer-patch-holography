"""Fail-closed custody for the complete physical local reference theorem."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from . import algebra,characters,physical,spectrum,reference

ROOT=Path(__file__).resolve().parents[4]
PACKAGE='model/coherent_defect_wave/chiral_joint/local_reference/'


def input_hashes():
    seeds={PACKAGE+name for name in ('__init__.py','__main__.py','algebra.py','characters.py',
        'reduction.py','physical.py','spectrum.py','reference.py','certificate.py',
        'test_reference.py','PROOF.md','README.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/higgs_composites/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/family_cell/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/neutral_local_spectrum/PROOF.md'))
    for name in sorted(seeds):
        if not (ROOT/name).is_file():raise FileNotFoundError('missing declared local reference input: '+name)
    paths=import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT/name).is_file():raise FileNotFoundError('missing local reference import input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema':'oph.full_physical_local_reference.v1',
        'source':{'Hamiltonian':'all native local Yukawa channels plus scalar excitation projector',
            'coefficients':{'up':'1','down':'1','charged_lepton':'1','scalar_penalty':'1'},
            'full_CAR_modes':30,'full_scalar_states':5,'full_unconstrained_dimension':5*2**30,
            'domain':'one cell, all twelve Gauss generators zero, no compensating boundary charge',
            'clock':'exp(-it H_local), hbar=1','spatial_or_composite_hopping_added':False},
        'exact_particle_hole_algebras':algebra.certificate(),
        'complete_quark_character_bound':characters.certificate(),
        'complete_maximal_module_physical_reduction':physical.certificate(),
        'global_physical_ground_and_gap':spectrum.certificate(),
        'selected_stationary_reference':reference.certificate(),
        'input_sha256':input_hashes()}
