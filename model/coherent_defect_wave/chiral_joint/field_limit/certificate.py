"""Source custody and exact finite low-energy interaction receipt."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from .model import certificate_data

ROOT=Path(__file__).resolve().parents[4]


def input_hashes():
    seeds={str(p.relative_to(ROOT)) for p in Path(__file__).parent.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/two_cell/SPECTRUM_PROOF.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/REFINEMENT_PROOF.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/PROPAGATION_PROOF.md',
        'model/coherent_defect_wave/chiral_source/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/periodic_family/model.py',
        'model/coherent_defect_wave/chiral_joint/periodic_family/PROOF.md',
        'Lean/QFT/EffectiveQuantumComparison.lean','Lean/lean-toolchain','Lean/lakefile.lean','Lean/lake-manifest.json'))
    missing=sorted(p for p in seeds if not (ROOT/p).is_file())
    if missing:raise ValueError('missing declared field-limit input: '+', '.join(missing))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def certificate():
    result=certificate_data()
    result['schema']='oph.coherent_defect_wave.chiral_joint.field_limit.v1'
    result['input_sha256']=input_hashes()
    return result
