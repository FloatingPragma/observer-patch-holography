"""Bounded source custody followed by independent semantic replay."""
import argparse
import hashlib
import json
from pathlib import Path
from .format import need, keys, equal, load, canonical
from . import source_check, observable_check

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OWN = ('__init__.py', 'format.py', 'source.py', 'source_check.py', 'observables.py',
       'observable_check.py', 'verify.py', 'build.py', 'test_paired_gravity.py',
       'README.md', 'CONTRACT.md', 'DATA.md', 'measurements.json')
SOURCES = ['code/paired_gravity/'+p for p in OWN]+[
    'requirements.txt', '.gitattributes',
    'extra/PAIRED_SOURCE_GRAVITY.md', '.github/workflows/paired-gravity.yml',
    'code/source_scalar_execution/source_scalar_execution_receipt.json',
    'code/source_scalar_execution/scalar_execution_algebra.py',
    'code/source_scalar_execution/verify_source_scalar_execution.py',
    'Lean/Screen/SeamCurrentCarrierQuotient.lean',
    'Lean/Screen/PrimitivePortFrameQuotient.lean', 'Lean/Screen/PortFrameGram.lean']
CLAIM = 'OPH-SOURCE-PAIRED-GRAVITY-NONIDENTIFICATION'
# Reviewed primary numerical transcription, independent of the producer/receipt.
DATA_SHA256 = '6a06eaef2cd6e28ca964557a7de29e5ece6b696c3af47839b2b31fe7fcfa9aa1'


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(SOURCES)}


def claim_digest():
    rows = json.loads((ROOT/'claims/claim_registry.yaml').read_text(encoding='utf-8'))['claims']
    found = [row for row in rows if row['claim_id'] == CLAIM]
    need(len(found) == 1, 'one registered paired-gravity claim')
    return hashlib.sha256(canonical(found[0])).hexdigest()


def data():
    path = HERE/'measurements.json'
    need(hashlib.sha256(path.read_bytes()).hexdigest() == DATA_SHA256, 'reviewed primary measurement transcription')
    return load(path)


def evidence(packet):
    keys(packet, 'source observables')
    # Cheap observable census precedes the full source replay.
    observable_check.verify(data(), packet['observables'], packet['source'])
    source_check.verify(packet['source'])


def verify(path=HERE/'receipt.json'):
    row = load(path)
    keys(row, 'schema sources claim evidence')
    equal(row['schema'], 'oph-paired-gravity-v1', 'receipt schema')
    equal(row['sources'], pins(), 'complete source custody')
    equal(row['claim'], claim_digest(), 'registered claim custody')
    evidence(row['evidence'])
    return row


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(p.parse_args().path)
    print('Verified paired-gravity source nonidentifiability and retrospective comparison')
