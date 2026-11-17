"""Bind exact source energy branches to their proof and actual gates."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import source_certificate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins

ROOT=Path(__file__).resolve().parents[3]


def certificate():
    local=Path(__file__).parent
    paths={p for p in local.glob('*.py') if not p.name.startswith('test_')}
    paths.add(local/'PROOF.md')
    paths.update(ROOT/name for name in (
        'model/coherent_defect_wave/energy/model.py',
        'model/coherent_defect_wave/energy/PROOF.md',
        'model/coherent_seam_repair/PROOF.md',
        'model/coherent_seam_repair/results.json',
        'model/requirements.txt',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    names=import_closure({str(p.relative_to(ROOT)) for p in paths})
    return {'schema':'oph.coherent_defect_wave.dispersion.v1',**source_certificate(),
            'input_sha256':bound_pins({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                            for name in sorted(names) if (ROOT/name).is_file()})}
