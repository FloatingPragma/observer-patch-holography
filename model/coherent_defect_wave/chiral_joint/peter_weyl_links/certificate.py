"""Fail-closed custody and complete finite Peter-Weyl source records."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from .model import descriptor
from .checks import (clebsch_checks, matrix_checks, path_checks,
                     rotation_word_checks, actual_hamiltonian_checks)
from .experiment import preparation_checks, transfer_checks

ROOT = Path(__file__).resolve().parents[4]
PACKAGE = 'model/coherent_defect_wave/chiral_joint/peter_weyl_links/'


def input_hashes():
    seeds = {PACKAGE + name for name in ('__init__.py', '__main__.py', 'clebsch.py',
        'model.py', 'checks.py', 'experiment.py', 'certificate.py', 'test_links.py',
        'README.md', 'PROOF.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/periodic_family/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/oriented_links/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/higgs_composites/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/neutral_mobility/PROOF.md'))
    for name in sorted(seeds):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing declared Peter-Weyl input: ' + name)
    paths = import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing Peter-Weyl import input: ' + name)
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema': 'oph.chiral_peter_weyl_links.receipt.v1',
        'source': descriptor(),
        'complete_native_Clebsch_derivation': clebsch_checks(),
        'full_link_matrices': matrix_checks(),
        'full_domain_path_controls': path_checks(),
        'complete_cubic_operator_words': rotation_word_checks(),
        'actual_occupied_Hamiltonian_controls': actual_hamiltonian_checks(),
        'physical_source_reference_preparation': preparation_checks(),
        'complete_neutral_transfer_experiment': transfer_checks(),
        'input_sha256': input_hashes()}
