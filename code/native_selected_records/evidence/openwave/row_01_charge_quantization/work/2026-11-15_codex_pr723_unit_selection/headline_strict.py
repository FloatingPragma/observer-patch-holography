"""Current small replay entry point; preserves the first registered replay."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


if __name__ == '__main__':
    registration = json.loads((HERE / 'VERIFY_PREREG.json').read_text())
    for name, digest in registration['files'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Changed verification input: ' + name)
    subprocess.run([sys.executable, str(HERE / 'headline.py')], check=True)
    subprocess.run([sys.executable, str(HERE / 'verify_strict.py')], check=True)
