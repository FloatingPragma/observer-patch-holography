"""Source-bound writer and complete cold replay for native electric links."""
from pathlib import Path
import json
from model.openwave_pure_su2_source.custody import pin_package,execute
from .source import run as source_run
HERE=Path(__file__).resolve().parent
def pins():
    return pin_package(HERE,('model/openwave_active_native_charge/PROOF.md',
      'model/openwave_native_process/PROOF.md','model/openwave_local_gauss_quantum/PROOF.md'))
def run():return json.loads(json.dumps(source_run()))
if __name__=='__main__':execute(HERE,pins,run)
