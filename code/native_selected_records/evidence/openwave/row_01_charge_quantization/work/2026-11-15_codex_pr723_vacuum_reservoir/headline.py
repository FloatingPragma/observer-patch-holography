"""Offline byte replay into temporary output; no sealed receipt is modified."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


if __name__ == '__main__':
    registration = json.loads((HERE/'PREREG.json').read_text())
    for name, digest in registration['files'].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest, name
    with tempfile.TemporaryDirectory(prefix='row01-reservoir-') as output:
        subprocess.run([sys.executable, str(HERE/'experiment.py'), '--output', output],
                       check=True, stdout=subprocess.DEVNULL)
        for name in ('RECEIPT.json', 'EVIDENCE.json.gz'):
            assert (HERE/name).read_bytes() == (Path(output)/name).read_bytes(), name
    receipt = json.loads((HERE/'RECEIPT.json').read_text())
    print(json.dumps(receipt, indent=2))
    assert receipt['passed']
