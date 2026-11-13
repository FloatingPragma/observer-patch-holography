from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,execute,pin_package
from .source import run as source_run
HERE=Path(__file__).resolve().parent
PARENTS=('model/native_gauss_walk/results.json','model/native_electric_gauss/results.json',
         'model/native_observer_program/results.json','model/openwave_m10_review_gate/results.json')
def pins():
    paths=set(PARENTS)
    for name in PARENTS:
        for path,digest in json.loads((ROOT/name).read_text())['source_sha256'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
            paths.add(path)
    return pin_package(HERE,paths)
def run():return json.loads(json.dumps(source_run()))
if __name__=='__main__':execute(HERE,pins,run)
