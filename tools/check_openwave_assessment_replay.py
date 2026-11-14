"""Compare complete local receipts across numerical runtimes with bounded roundoff."""
import argparse
import json
import math
from pathlib import Path


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
        if not ('.measurement.' in path or path.startswith('receipt.observer_protocol.')) and expected != actual:
            raise ValueError('declared float mismatch: '+path)
        if not math.isfinite(expected) or not math.isfinite(actual) or not math.isclose(expected,actual,rel_tol=1e-10,abs_tol=1e-10):
            raise ValueError('numeric mismatch: '+path)
    elif expected!=actual:raise ValueError('value mismatch: '+path)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('expected',type=Path);parser.add_argument('actual',type=Path)
    args=parser.parse_args();compare(json.loads(args.expected.read_text()),json.loads(args.actual.read_text()))
    print('PASS_COMPLETE_SOURCE_BOUND_ASSESSMENT_REPLAY')
