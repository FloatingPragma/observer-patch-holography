"""Only frozen producer bytes are imported by the local executable."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE/'inputs'
PINS = json.loads((HERE/'INPUTS.json').read_text())['files_sha256']


def verify():
    for path, expected in PINS.items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
            raise ArithmeticError('frozen producer bytes changed: '+path)


verify()
sys.path[:0] = [str(ROOT/'openwave-tuning/common_dynamics/interacting_gaussian'),
               str(ROOT/'openwave-tuning/next/CHAIN2/local12_limit_audit')]
import branches as source
import certified_bounds as bounds
import fock_transfer as total_cap
import refinement_adapter as native
