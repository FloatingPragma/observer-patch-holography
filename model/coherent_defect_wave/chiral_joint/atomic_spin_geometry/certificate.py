"""Recompute the complete spin spectrum, actual momentum and source transfer."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .cube import full_operators, spectrum
from .momentum import certificate as momentum_certificate
from .transfer import certificate as transfer_certificate

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).parent
DEPENDENCIES = tuple('model/coherent_defect_wave/chiral_joint/'+name+'/'+file
                     for name in ('atomic_reference_control', 'atomic_even_remainder',
                                  'atomic_momentum_response', 'local_reference',
                                  'charged_ground_exchange')
                     for file in ('PROOF.md', 'results.json'))


def input_hashes():
    seeds = {str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    seeds.update(DEPENDENCIES)
    seeds.add('model/requirements.txt')
    for name in seeds:
        if not (ROOT/name).is_file():
            raise FileNotFoundError('missing declared spin geometry input: '+name)
    closure = import_closure(seeds)
    for name in closure:
        if not (ROOT/name).is_file():
            raise FileNotFoundError('missing spin geometry import: '+name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(closure)}


def certificate():
    full_operators()
    return {'schema': 'oph.coherent_defect_wave.chiral_joint.atomic_spin_geometry.v1',
            'full_spin_spectrum': spectrum(), 'actual_source_momentum': momentum_certificate(),
            'full_interacting_transfer': transfer_certificate(), 'input_sha256': input_hashes()}
