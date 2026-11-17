"""Complete byte replay and independent record verification in under60seconds."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


if __name__ == '__main__':
    for registration in ('VERIFY_PREREG.json','MUTATION_PREREG.json'):
        for name,digest in json.loads((HERE/registration).read_text())['files'].items():
            assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest,name
    subprocess.run([sys.executable,str(HERE/'headline.py')],check=True)
    subprocess.run([sys.executable,str(HERE/'verify.py')],check=True)
    controls = subprocess.run([sys.executable,str(HERE/'controls.py')],check=True,
                              capture_output=True,text=True)
    assert json.loads(controls.stdout) == json.loads((HERE/'CONTROLS.json').read_text())
