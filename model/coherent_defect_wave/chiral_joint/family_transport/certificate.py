"""Complete common-source family records with fail-closed custody."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from .model import descriptor
from .checks import dictionary_checks, preparation_and_detector_checks, full_powers

ROOT = Path(__file__).resolve().parents[4]
PACKAGE = 'model/coherent_defect_wave/chiral_joint/family_transport/'


def input_hashes():
    seeds = {PACKAGE + name for name in ('__init__.py', '__main__.py', 'model.py',
        'checks.py', 'certificate.py', 'test_transport.py', 'README.md', 'PROOF.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/gauge_defect_transport/PROOF.md',
        'model/coherent_defect_wave/chiral_source/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/peter_weyl_links/PROOF.md'))
    for name in sorted(seeds):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing declared family-transport input: ' + name)
    paths = import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing family-transport import input: ' + name)
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema': 'oph.common_source_native_family_transport.receipt.v1',
        'source_and_static_reference': descriptor(),
        'actual_native_family_dictionary': dictionary_checks(),
        'physical_preparations_and_complete_detectors': preparation_and_detector_checks(),
        'same_clock_full_H_probability_records': full_powers(),
        'input_sha256': input_hashes()}
