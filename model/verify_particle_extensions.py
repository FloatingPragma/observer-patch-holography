"""Additional live particle producers, independent of retained replay ordering."""
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
PACKAGES=(
    'openwave_beta_angular', 'openwave_lepton_widths', 'openwave_neutrino_matter',
    'openwave_spin_gradient_nuclei', 'openwave_rest_spectroscopy', 'openwave_colour_string',
    'openwave_common_electromagnetism', 'openwave_boosted_records',
    'openwave_single_valley_leptons', 'openwave_mass_catalogue',
    'openwave_longitudinal_review', 'openwave_boost_review',
)

def source_bindings():
    """Reconstruct each producer's live import/receipt inventory, fail closed."""
    for name in PACKAGES:
        module=importlib.import_module(f'model.{name}.__main__')
        expected=module.pins()
        receipt=json.loads((ROOT/'model'/name/'results.json').read_text())
        if not expected or receipt.get('source_sha256')!=expected:
            raise ValueError(f'particle extension source inventory changed: {name}')
        if not str(receipt.get('verdict','')).startswith('PASS'):
            raise ValueError(f'particle extension producer failed: {name}')
    return len(PACKAGES)

def selected(index=None,count=None):
    if index is None and count is None:return PACKAGES
    if index is None or count is None or count<1 or not 0<=index<count:
        raise ValueError('both valid shard index and count required')
    return tuple(name for i,name in enumerate(PACKAGES)if i%count==index)

def replay(index=None,count=None):
    names=selected(index,count)
    environment=dict(os.environ)
    for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):environment.setdefault(k,'1')
    for name in names:
        argv=[sys.executable,'-m',f'model.{name}','--check']
        if name=='openwave_boosted_records':
            upstream=environment.get('OPH_OPENWAVE_ROOT')
            python=environment.get('OPH_TAICHI_PYTHON')
            if not upstream or not python:
                raise ValueError('boosted replay requires pinned OPH_OPENWAVE_ROOT and OPH_TAICHI_PYTHON; no numerical fallback')
            argv+=['--openwave-root',upstream,'--python',python]
        subprocess.run(argv,cwd=ROOT,env=environment,check=True)
    return len(names)
