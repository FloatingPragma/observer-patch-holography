"""Bind the actual collective operator, exact residues and full-source transfer."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import operator_certificate
from .residue import certificate as residue_certificate
from .transfer import certificate as transfer_certificate

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).parent
DEPENDENCIES = tuple('model/coherent_defect_wave/chiral_joint/'+name+'/'+file
    for name in ('atomic_spin_geometry', 'atomic_reference_control',
                 'atomic_even_remainder', 'atomic_momentum_response',
                 'local_reference', 'charged_ground_exchange')
    for file in ('PROOF.md', 'results.json'))


def input_hashes():
    seeds = {str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    seeds.update(DEPENDENCIES)
    seeds.add('model/requirements.txt')
    for name in seeds:
        if not (ROOT/name).is_file():
            raise FileNotFoundError('missing declared collective spin input: '+name)
    closure = import_closure(seeds)
    for name in closure:
        if not (ROOT/name).is_file():
            raise FileNotFoundError('missing collective spin import: '+name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(closure)}


def certificate():
    return {'schema': 'oph.actual_interacting_collective_spin_response.v1',
        'actual_native_operator': operator_certificate(),
        'exact_atomic_residue': residue_certificate(),
        'same_source_interacting_records': transfer_certificate(),
        'input_sha256': input_hashes()}
