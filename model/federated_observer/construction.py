"""Join the retained federation, current predictions and refinement checks."""
from __future__ import annotations

import hashlib
from pathlib import Path

from . import correlations, memory, refinement, scaling, source


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((
    *source.SOURCE_PATHS, *correlations.SOURCE_PATHS,
    *refinement.SOURCE_PATHS, *scaling.SOURCE_PATHS, *memory.SOURCE_PATHS,
    'model/federated_observer/__init__.py',
    'model/federated_observer/__main__.py',
    'model/federated_observer/construction.py',
    'model/observer/scalar_response/experiment.py',
    'tools/strict_json.py', 'model/requirements.txt',
))))


def run():
    return {
        'schema': 'oph.federated-observer.v1',
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                          for path in SOURCE_PATHS},
        'construction': source.run(),
        'current_predictions': correlations.run(),
        'refinement': refinement.run(),
        'coarse_prediction_with_memory': memory.run(),
        'scaling': scaling.run(),
        'scope': {
            'primitive': 'the same twelve-port observer and binary endpoint repair',
            'larger_objects': 'coupled federations, not enlarged elementary carriers',
            'screen_levels_executed': [0, 1, 2],
            'gluing': 'retained port_pair convention, not uniquely selected by OPH',
            'physical_global_clock_required': False,
            'controlled_relativistic_spacetime_limit_proved': False,
            'new_empirical_measurements_used': False,
            'new_lean_formalization_claimed': False,
        },
    }
