import json
from pathlib import Path
from model.openwave_pure_su2_source.custody import execute,pin_package,ROOT
from .source import run,PARENT
HERE=Path(__file__).resolve().parent
def pins():
    return pin_package(HERE,{PARENT,*json.loads((ROOT/PARENT).read_text())['source_sha256']})
if __name__=='__main__':execute(HERE,pins,run)
