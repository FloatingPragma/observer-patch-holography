"""Original full-source writer and independent cold comparison."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run
HERE=Path(__file__).resolve().parent
PARENTS=['model/openwave_lambda_spectral_measure/results.json','model/openwave_baryon_poles/results.json','model/openwave_cluster_weyl_sequence/results.json']
def pins():
    paths=set(PARENTS)|{str((HERE/'witnesses.json').relative_to(ROOT))}
    for name in PARENTS:
        for path,expected in json.loads((ROOT/name).read_text())['source_sha256'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==expected,('stale parent',path)
            paths.add(path)
    return pin_package(HERE,paths)
if __name__=='__main__':execute(HERE,pins,run)
