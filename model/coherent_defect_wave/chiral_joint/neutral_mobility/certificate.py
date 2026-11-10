"""Recomputed exact receipt and source ancestry for neutral odd arrival."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from .model import descriptor
from .checks import local_preparation_checks, full_domain_checks, rotation_checks, experiment

ROOT = Path(__file__).resolve().parents[4]


def input_hashes():
    here = Path(__file__).parent
    seeds = {str(path.relative_to(ROOT)) for path in here.iterdir()
             if path.suffix in ('.py', '.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/oriented_links/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/higgs_composites/PROOF.md'))
    for name in sorted(seeds):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing declared neutral-mobility input: ' + name)
    paths = import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing neutral-mobility import input: ' + name)
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema': 'oph.chiral_neutral_mobility.receipt.v1',
            'source': descriptor(),
            'physical_preparation_and_local_detector': local_preparation_checks(),
            'full_domain_operator_identities': full_domain_checks(),
            'spatial_symmetries': rotation_checks(),
            'complete_physical_arrival_experiment': experiment(),
            'input_sha256': input_hashes()}
