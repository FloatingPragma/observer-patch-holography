from pathlib import Path
from model.openwave_pure_su2_source.custody import execute, pin_package
from .spectrum import run
HERE = Path(__file__).resolve().parent
def pins():
    return pin_package(HERE, ('model/openwave_su2_spectral/numerical_witnesses.json',
                             'model/openwave_pure_su2_source/results.json',
                             'model/openwave_pure_su2_source/PROOF.md',
                             'model/particle_yang_mills_gap_contract/PROOF.md'))
if __name__ == '__main__':
    execute(HERE, pins, run)
