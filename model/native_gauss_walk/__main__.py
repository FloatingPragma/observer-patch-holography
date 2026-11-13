"""Strict source-bound replay of the joint charged native field state."""
from pathlib import Path
import hashlib,json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run as source_run
HERE=Path(__file__).resolve().parent
def pins():
    paths=set()
    for parent in('model/native_electric_gauss/results.json','model/native_observer_program/results.json'):
        paths.add(parent)
        for p,h in json.loads((ROOT/parent).read_text())['source_sha256'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale parent',p)
            paths.add(p)
    return pin_package(HERE,paths)
def run():return json.loads(json.dumps(source_run()))
if __name__=='__main__':execute(HERE,pins,run)
