"""Bind the exact correspondence calculation and its imported operators."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import result

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent


def certificate():
    seeds = {str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.add('model/requirements.txt')
    value = result()
    value['input_sha256'] = {
        p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
        for p in sorted(import_closure(seeds))}
    return value
