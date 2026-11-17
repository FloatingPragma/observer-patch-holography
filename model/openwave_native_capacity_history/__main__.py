import json
from pathlib import Path
from model.openwave_pure_su2_source.custody import execute,pin_package,ROOT
from .source import run,PARENT
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent
def pins():
    return pin_package(HERE,{PARENT,*bound_pins(json.loads((ROOT/PARENT).read_text())['source_sha256'])})
if __name__=='__main__':execute(HERE,pins,run)
