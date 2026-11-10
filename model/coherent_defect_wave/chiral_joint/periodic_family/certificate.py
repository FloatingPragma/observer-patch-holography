"""Recomputed exact periodic source receipt and explicit input custody."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from .checks import certificate_data

ROOT = Path(__file__).resolve().parents[4]


def input_hashes():
    here = Path(__file__).parent
    seeds = {str(p.relative_to(ROOT)) for p in here.iterdir() if p.suffix in ('.py', '.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_source/PROOF.md',
        'model/coherent_defect_wave/chiral_regulator/PROOF.md',
        'model/coherent_defect_wave/chiral_higgs/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/SYMMETRIES.md',
        'Lean/QFT/PeriodicMomentum.lean', 'Lean/lean-toolchain',
        'Lean/lakefile.lean', 'Lean/lake-manifest.json'))
    paths = import_closure(seeds)
    missing = [name for name in paths if not (ROOT / name).is_file()]
    if missing:
        raise ValueError('missing declared source inputs: ' + repr(sorted(missing)))
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {**certificate_data(), 'input_sha256': input_hashes()}
