"""Bounded primary input, complete custody and independent semantic replay."""
import argparse
import hashlib
import json
from pathlib import Path

from .check import need, keys, same, verify_evidence


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CLAIM = 'OPH-SOURCE-ADMISSIBLE-CLOCK-SPATIAL-NONSELECTION'
OWN = ('__init__.py', 'model.py', 'check.py', 'verify.py', 'build.py',
       'test_admissibility.py', 'README.md', 'CONTRACT.md')
SOURCES = ['code/source_gravity_admissibility/'+p for p in OWN]+[
    'code/source_selection_model/response.json',
    'code/source_selection_model/verify_response.py',
    'code/source_selection_model/DERIVATION.md',
    'code/source_selection_model/RECORD_GLUING.md',
    'code/m1_source_realization/model.py', 'code/m1_source_realization/check.py',
    'code/a5_closure/port_current_inner_certificate.py',
    'claims/axiom_registry.yaml']


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('ascii')


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(SOURCES)}


def claim_pin():
    data = json.loads((ROOT/'claims/claim_registry.yaml').read_text(encoding='utf-8'))
    rows = [r for r in data['claims'] if r['claim_id'] == CLAIM]
    need(len(rows) == 1, 'one scientific claim')
    # Bind the consumed definitions, without invalidating this receipt when
    # an unrelated assumption is added elsewhere in the shared dictionary.
    lines = (ROOT/'claims/assumption_dictionary.md').read_text(encoding='utf-8').splitlines()
    definitions = {}
    for name in rows[0]['assumptions']:
        found = [line for line in lines if line.startswith('| `'+name+'` |')]
        need(len(found) == 1, 'one definition for each consumed assumption')
        definitions[name] = found[0]
    return hashlib.sha256(canonical(dict(claim=rows[0], assumptions=definitions))).hexdigest()


def load(path):
    need(Path(path).stat().st_size < 750_000, 'bounded receipt')
    raw = Path(path).read_bytes()
    need(len(raw) < 750_000, 'bounded receipt')
    def unique(items):
        out = {}
        for k, v in items:
            need(k not in out, 'duplicate JSON key')
            out[k] = v
        return out
    def forbidden(value):
        raise ValueError('nonfinite JSON token')
    def integer(value):
        need(len(value) < 100, 'bounded JSON integer')
        return int(value)
    try:
        return json.loads(raw.decode('ascii'), object_pairs_hook=unique,
                          parse_constant=forbidden, parse_int=integer)
    except (RecursionError, UnicodeError) as exc:
        raise ValueError('invalid receipt encoding or nesting') from exc


def verify(path=HERE/'receipt.json'):
    packet = load(path)
    keys(packet, 'schema sources claim evidence')
    same(packet['schema'], 'oph-source-gravity-admissibility-v1')
    same(packet['sources'], pins(), 'source custody')
    same(packet['claim'], claim_pin(), 'scientific claim custody')
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified admitted source controls, identical clocks and distinct finite spatial reads')
