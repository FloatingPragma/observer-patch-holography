"""Stable algebraic receipt for the complete oriented-link construction."""
import hashlib
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from .model import descriptor
from .checks import (matrix_checks, peter_weyl_check, path_checks,
                     rotation_word_checks, actual_hamiltonian_checks)
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[4]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[4]


def input_hashes():
    here = Path(__file__).parent
    seeds = {str(path.relative_to(ROOT)) for path in here.iterdir()
             if path.suffix in ('.py', '.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
                  'model/coherent_defect_wave/chiral_source/PROOF.md',
                  'model/coherent_defect_wave/chiral_higgs/PROOF.md',
                  'model/coherent_defect_wave/chiral_joint/periodic_family/PROOF.md'))
    for name in sorted(seeds):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing declared oriented-link input: ' + name)
    paths = import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing oriented-link import input: ' + name)
    return bound_pins({name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)})


def certificate():
    return {
        'schema': 'oph.chiral_oriented_quantum_links.receipt.v1',
        'source': descriptor(),
        'full_link_matrices': matrix_checks(),
        'Peter_Weyl_distinction': peter_weyl_check(),
        'path_norms': path_checks(),
        'proper_cubic_rotations': rotation_word_checks(),
        'actual_full_carrier_executions': actual_hamiltonian_checks(),
        'input_sha256': bound_pins(input_hashes()),
    }
