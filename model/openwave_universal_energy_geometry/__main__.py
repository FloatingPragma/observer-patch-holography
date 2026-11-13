from pathlib import Path
from model.openwave_pure_su2_source.custody import execute, pin_package
from .source import run
HERE = Path(__file__).resolve().parent
def pins():
    return pin_package(HERE, ('model/openwave_universal_energy_geometry/profile_fixture.json',
                             'model/openwave_native_force/PROOF.md',
                             'model/openwave_native_force/results.json',
                             'model/openwave_energy_lapse/PROOF.md',
                             'model/finite_observer_readout_closure/PROOF.md'))
if __name__ == '__main__':
    execute(HERE, pins, run)
