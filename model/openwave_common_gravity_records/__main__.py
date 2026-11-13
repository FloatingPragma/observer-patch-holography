from pathlib import Path
from model.openwave_pure_su2_source.custody import execute, pin_package
from .records import run
HERE = Path(__file__).resolve().parent
def pins():
    return pin_package(HERE, ('model/openwave_universal_energy_geometry/results.json',
                             'model/openwave_universal_energy_geometry/profile_fixture.json',
                             'model/openwave_universal_energy_geometry/PROOF.md',
                             'model/openwave_joint_native_clocks/PROOF.md',
                             'model/openwave_joint_native_clocks/clocks.py',
                             'model/finite_observer_readout_closure/PROOF.md'))
if __name__ == '__main__':
    execute(HERE, pins, run)
