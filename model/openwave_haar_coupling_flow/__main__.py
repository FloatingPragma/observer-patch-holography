from pathlib import Path
from model.openwave_pure_su2_source.custody import execute,pin_package
from .flow import run
HERE=Path(__file__).resolve().parent
def pins():
    return pin_package(HERE,('model/openwave_haar_coupling_flow/numerical_witnesses.json', 'model/openwave_coupled_haar_spectrum/PROOF.md', 'model/openwave_coupled_haar_source/results.json'))
if __name__=="__main__":execute(HERE,pins,run)
