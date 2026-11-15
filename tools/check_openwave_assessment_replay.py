"""Compare complete local receipts across numerical runtimes with bounded roundoff."""
import argparse
import json
import math
import re
from pathlib import Path


# These fields declare the experiment, even when embedded below a measured
# result. Runtime roundoff never authorizes changing a preparation or clock.
DECLARED_FLOAT_FIELDS = frozenset({
    'mass', 'm', 'energy_mass', 'energy_coupling', 'coupling', 'time', 'times',
    'evolution_time', 'time_per_step', 'step_time', 'times_over_gamma',
    'commit_probability', 'commit_probabilities', 'gamma', 'gammas', 'rate',
    'rates', 'seed', 'cutoff', 'radius', 'tolerance', 'atol', 'rtol',
})


def measured_float_path(path):
    tokens = [token for token in re.split(r'[.\[\]]+', path) if token and not token.isdigit()]
    if 'supplied' in tokens or 'preparation' in tokens:
        return False
    if any(token in DECLARED_FLOAT_FIELDS or token.endswith(('_time', '_mass', '_coupling')) for token in tokens):
        return False
    return ('.measurement.' in path or path.startswith('receipt.observer_protocol.')
            or path.startswith('receipt.record_extensions.'))


def compare(expected, actual, path='receipt'):
    if path=='receipt.versions':
        if set(expected)!=set(actual) or any(not isinstance(v,str) for v in actual.values()):
            raise ValueError('invalid runtime provenance')
        return
    if type(expected) is not type(actual):
        raise ValueError('type mismatch: '+path)
    if isinstance(expected,dict):
        if set(expected)!=set(actual):raise ValueError('key mismatch: '+path)
        for key in expected:compare(expected[key],actual[key],path+'.'+key)
    elif isinstance(expected,list):
        if len(expected)!=len(actual):raise ValueError('length mismatch: '+path)
        for i,(a,b) in enumerate(zip(expected,actual)):compare(a,b,path+f'[{i}]')
    elif isinstance(expected,float):
        if not measured_float_path(path) and expected != actual:
            raise ValueError('declared float mismatch: '+path)
        if not math.isfinite(expected) or not math.isfinite(actual) or not math.isclose(expected,actual,rel_tol=1e-10,abs_tol=1e-10):
            raise ValueError('numeric mismatch: '+path)
    elif expected!=actual:raise ValueError('value mismatch: '+path)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('expected',type=Path);parser.add_argument('actual',type=Path)
    args=parser.parse_args();compare(json.loads(args.expected.read_text()),json.loads(args.actual.read_text()))
    print('PASS_COMPLETE_SOURCE_BOUND_ASSESSMENT_REPLAY')
