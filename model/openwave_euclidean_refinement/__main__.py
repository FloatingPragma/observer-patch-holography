from pathlib import Path
from model.openwave_pure_su2_source.custody import execute, pin_package
from .response import run
HERE = Path(__file__).resolve().parent
def pins():
    return pin_package(HERE, ('model/openwave_coupled_haar_spectrum/results.json',
                             'model/openwave_coupled_haar_spectrum/numerical_witnesses.json',
                             'model/openwave_coupled_haar_spectrum/PROOF.md',
                             'model/particle_yang_mills_gap_contract/PROOF.md',
                             'model/openwave_su2_transfer/PROOF.md',
                             'model/openwave_su2_propagation/PROOF.md',
                             'Lean/ObserverPatchHolography/DilationGapObstruction.lean'))
if __name__ == '__main__':
    execute(HERE, pins, run)
