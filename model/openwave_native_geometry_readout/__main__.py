from pathlib import Path
from model.openwave_pure_su2_source.custody import execute,pin_package
from .source import run
HERE=Path(__file__).resolve().parent
def pins():
    return pin_package(HERE,{'code/source_selection_model/response.json',
                            'code/m1_source_realization/receipt.json'})
if __name__=='__main__':execute(HERE,pins,run)
