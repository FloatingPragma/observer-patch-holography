from pathlib import Path
from .source import run
from .custody import execute, pin_package
HERE = Path(__file__).resolve().parent
def pins():
    return pin_package(HERE, ('model/particle_yang_mills_gap_contract/PROOF.md',
                             'model/particle_yang_mills_gap_contract/transfer.py',
                             'Lean/ObserverPatchHolography/RepairGapChain.lean'))
if __name__ == '__main__':
    execute(HERE, pins, run)
