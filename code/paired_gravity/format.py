"""Bounded, canonical receipt input; checks survive python -O."""
import hashlib
import json
import re
from fractions import Fraction as F
from pathlib import Path


def need(ok, message):
    if not ok:
        raise ValueError(message)


def keys(obj, names):
    need(type(obj) is dict and set(obj) == set(names.split()), 'exact object keys')


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('ascii')


def digest(obj):
    return hashlib.sha256(canonical(obj)).hexdigest()


def equal(actual, expected, label):
    need(canonical(actual) == canonical(expected), label)


def rational(value):
    need(type(value) is str and len(value) <= 160 and
         re.fullmatch(r'-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?', value) is not None,
         'bounded rational string')
    try:
        x = F(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError('invalid rational') from exc
    need(str(x) == value and max(x.numerator.bit_length(), x.denominator.bit_length()) <= 512,
         'canonical bounded rational')
    return x


def load(path):
    need(Path(path).stat().st_size <= 1_000_000, 'receipt size limit')
    with Path(path).open('rb') as stream:
        raw = stream.read(1_000_001)
    need(len(raw) <= 1_000_000, 'receipt size limit')
    def pairs(items):
        out = {}
        for k, v in items:
            need(k not in out, 'duplicate JSON key')
            out[k] = v
        return out
    def forbidden(value):
        raise ValueError('floating/nonfinite JSON token forbidden')
    def integer(value):
        need(len(value) <= 100, 'bounded JSON integer')
        return int(value)
    return json.loads(raw.decode('ascii'), object_pairs_hook=pairs,
                      parse_int=integer, parse_float=forbidden, parse_constant=forbidden)
