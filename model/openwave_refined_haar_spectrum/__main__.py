from pathlib import Path
from model.openwave_pure_su2_source.custody import execute,pin_package
from .certificate import run
HERE=Path(__file__).resolve().parent
def pins():
    return pin_package(HERE,('model/openwave_refined_haar_spectrum/numerical_witnesses.json','model/openwave_refined_haar_source/PROOF.md','model/openwave_cubic_haar_source/PROOF.md','model/openwave_su2_transfer/PROOF.md'))
if __name__=='__main__':execute(HERE,pins,run)
