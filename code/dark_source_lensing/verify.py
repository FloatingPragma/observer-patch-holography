"""Verify custody and independently replay the complete mathematical certificate."""
import argparse
import hashlib
import json
from pathlib import Path
from . import check
from .format import digest, equal, keys, load, need
import pin_scope

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CLAIM = 'OPH-DARK-SOURCE-DYNAMICS-LENSING'
OWN = ('__init__.py', 'format.py', 'model.py', 'check.py', 'global_ray.py', 'verify.py', 'build.py',
       'test_lensing.py', 'README.md', 'CONTRACT.md')
SOURCES = ['code/dark_source_lensing/'+p for p in OWN]+[
    'docs/research/DARK_SOURCE_DYNAMICS_LENSING.md', '.github/workflows/dark-source-lensing.yml',
    'requirements.txt', '.gitattributes',
    'Lean/ObserverPatchHolography/EinsteinBranch/DarkSector.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/DeepProfileClosure.lean']
SOURCES = pin_scope.bound(SOURCES)


def pins():
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(SOURCES)}


def claim_digest():
    rows = json.loads((ROOT/'claims/claim_registry.yaml').read_text(encoding='utf-8'))['claims']
    own = [row for row in rows if row['claim_id'] == CLAIM]
    need(len(own) == 1, 'one registered source/lensing claim')
    return digest(own[0])


def verify(path=HERE/'receipt.json'):
    row = load(path)
    keys(row, 'schema sources claim evidence')
    equal(row['schema'], 'oph-dark-source-lensing-v1', 'receipt schema')
    equal(pin_scope.bound_pins(row['sources']), pins(), 'complete source custody')
    equal(row['claim'], claim_digest(), 'registered claim custody')
    check.verify(row['evidence'])
    return row


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified exact source/stress bounds, null observables and global nonuniqueness')
