from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT, execute, pin_package
from .source import run

HERE = Path(__file__).resolve().parent
PARENTS = ('model/openwave_native_response/results.json',
           'model/openwave_native_event_clock/results.json')


def pins():
    paths = set(PARENTS)
    for name in PARENTS:
        for path, digest in json.loads((ROOT/name).read_text())['source_sha256'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
            paths.add(path)
    return pin_package(HERE, paths)


if __name__ == '__main__':
    execute(HERE, pins, run)
