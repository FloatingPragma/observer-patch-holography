"""Original source-service writer and cold replay with complete source pins."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run
HERE=Path(__file__).resolve().parent

def pins():
    paths={'model/common_source/repair_continuation.py',
           'code/source_selection_model/response.json','code/m1_source_realization/receipt.json'}
    for parent in ('model/openwave_active_native_charge/results.json',
                   'model/openwave_coherent_native_architecture/results.json'):
        paths.add(parent)
        for name,h in json.loads((ROOT/parent).read_text())['source_sha256'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,('stale charge or architecture',name)
            paths.add(name)
    return pin_package(HERE,paths)

if __name__=='__main__':
    execute(HERE,pins,lambda:json.loads(json.dumps(run())))
