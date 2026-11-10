"""Fail-closed source custody and exact physical local energy records."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from .model import spectrum
from .symmetries import certificate as symmetry_certificate

ROOT = Path(__file__).resolve().parents[4]
PACKAGE = 'model/coherent_defect_wave/chiral_joint/neutral_local_spectrum/'


def input_hashes():
    seeds = {PACKAGE + name for name in ('__init__.py', '__main__.py', 'model.py',
        'symmetries.py', 'certificate.py', 'test_spectrum.py', 'PROOF.md', 'README.md',
        'equitable_partition.npz')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/family_cell/equitable_partition.npz',
        'model/coherent_defect_wave/chiral_joint/family_cell/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/higgs_composites/PROOF.md'))
    for name in sorted(seeds):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing declared neutral local input: ' + name)
    paths = import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing neutral local import input: ' + name)
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema': 'oph.full_native_neutral_local_spectrum.v1',
            'source': {'Hamiltonian': 'all three native local Yukawa channels plus scalar excitation projector',
                       'coefficients': {'up': '1', 'down': '1', 'charged_lepton': '1', 'scalar_penalty': '1'},
                       'full_CAR_modes': 30, 'scalar_states': 5,
                       'input': 'D_up^dagger Omega, with D=J/sqrt(2)',
                       'clock': 'exp(-it H_local), hbar=1',
                       'spatial_link_or_composite_hopping_added': False},
            'physical_symmetries': symmetry_certificate(),
            'complete_physical_spectrum': spectrum(),
            'input_sha256': input_hashes()}
