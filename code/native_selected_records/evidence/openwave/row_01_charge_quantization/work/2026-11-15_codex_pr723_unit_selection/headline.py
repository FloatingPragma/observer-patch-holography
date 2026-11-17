"""Offline native regeneration, byte comparison, and independent verification."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def main():
    registration = json.loads((HERE / 'PREREG.json').read_text())
    for name, digest in registration['files'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Changed registered input: ' + name)
    with tempfile.TemporaryDirectory(prefix='row01-unit-selection-') as tmp:
        output = Path(tmp)
        subprocess.run([sys.executable, str(HERE / 'experiment.py'), '--output', tmp],
                       check=True, stdout=subprocess.DEVNULL)
        for name in ('RECEIPT.json', 'EVIDENCE.json.gz'):
            if (output / name).read_bytes() != (HERE / name).read_bytes():
                raise ValueError('Replay changed: ' + name)
        subprocess.run([sys.executable, str(HERE / 'verify.py'), tmp], check=True)
    print('Byte-identical native replay; finite unit-selection claim only; row remains untested.')


if __name__ == '__main__':
    main()
