"""Production execution and rejection gates must survive Python optimization."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


def optimized(script):
    process = subprocess.run(
        [sys.executable, '-B', '-O', '-c', script],
        cwd=Path(__file__).resolve().parents[3],
        env={**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'},
        capture_output=True, text=True, timeout=60,
    )
    assert process.returncode == 0, process.stdout + process.stderr


def test_optimized_layers_cover_cell_and_federation_with_matchings():
    optimized('''
from model.executable_model import carrier
from model.candidate_laws.contrast_memory_spherical_v1 import sectors
for car, count in ((carrier.cell(), 5), (carrier.federation(0), 6)):
    layers = sectors.layers(car)
    if len(layers) != count:
        raise RuntimeError('wrong layer count')
    if sorted(i for layer in layers for i in layer) != list(range(len(car.seams))):
        raise RuntimeError('layers do not cover every seam exactly once')
    for layer in layers:
        ports = [p for i in layer for p in car.seams[i]]
        if len(ports) != len(set(ports)):
            raise RuntimeError('a layer is not a matching')
''')


def test_optimized_carrier_rejects_tampered_degree_table():
    optimized('''
from model.executable_model import carrier
carrier._federation_arrays = lambda level: (12, carrier.CELL_SEAMS, (6,) * 12)
try:
    carrier.federation(0)
except ValueError as error:
    if str(error) != 'stored degrees disagree with the seams':
        raise
else:
    raise RuntimeError('tampered degree table accepted')
''')


def test_optimized_exact_law_rejects_changed_write():
    optimized('''
from model.candidate_laws.contrast_memory_spherical_v1 import law
law.native_event = lambda a, b, c: (a, b)
try:
    law.exact_checks()
except ValueError as error:
    if not isinstance(error.args[0], dict) or not error.args[0]:
        raise
else:
    raise RuntimeError('changed write accepted by exact law gate')
''')
