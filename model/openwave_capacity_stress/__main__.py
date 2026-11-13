from pathlib import Path
from model.openwave_pure_su2_source.custody import execute, pin_package
from .source import run

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def pins():
    extra = [str(p.relative_to(ROOT)) for p in (HERE/'inputs').iterdir() if p.is_file()]
    extra += ['model/openwave_capacity_stress/comparison_protocol.json',
              'model/openwave_capacity_stress/Transfer.lean',
              'Lean/ObserverPatchHolography/CapacityStress.lean',
              'Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean',
              'Lean/ObserverPatchHolography/CapacityNonidentifiability.lean']
    return pin_package(HERE,extra)


if __name__ == '__main__':
    execute(HERE,pins,run)
