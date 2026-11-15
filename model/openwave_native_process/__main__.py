"""Write or cold-replay the complete native process receipt."""
from pathlib import Path
import json
from model.openwave_pure_su2_source.custody import pin_package, execute
from .source import run, federation

HERE = Path(__file__).resolve().parent


def pins():
    return pin_package(HERE, (*federation.SOURCE_PATHS,
        'evidence/native_geometric_source_20260925/spec.json',
        'Lean/EventAlgebra/SlotLocalitySurface.lean',
        'Lean/ObserverPatchHolography/Locality/DependencyCone.lean'))


if __name__ == '__main__':
    execute(HERE, pins, lambda: json.loads(json.dumps(run())))
