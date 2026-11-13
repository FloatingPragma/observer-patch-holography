"""Original producer and strict source-bound cold replay."""
from pathlib import Path
import hashlib
import json
from model.openwave_pure_su2_source.custody import ROOT, pin_package, execute
from .source import certificate

HERE = Path(__file__).resolve().parent
PARENTS = ('model/openwave_clock_phase_resource/results.json',
           'model/openwave_passive_clock_records/results.json')


def pins():
    extra = list(PARENTS)
    for parent in PARENTS:
        for name, digest in json.loads((ROOT/parent).read_text())['source_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, ('stale parent', name)
            extra.append(name)
    return pin_package(HERE, extra)


if __name__ == '__main__':
    execute(HERE, pins, certificate)
