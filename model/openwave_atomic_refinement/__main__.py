"""Regenerate all source Green operators and finite atom eigenpairs."""
import argparse,hashlib,json
from pathlib import Path
from .refinement import run
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def produce():
    result=run()
    paths=[*sorted(HERE.glob('*.py')),*sorted(HERE.glob('*.md')),
        ROOT/'model/openwave_field_observables/orbitals.py']
    result['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    return result


def equal(a,b,key=''):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v,k) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y,key) for x,y in zip(a,b))
    if isinstance(b,float):return isinstance(a,(int,float)) and abs(a-b)<=max(1e-9,abs(b)*1e-6)
    return type(a) is type(b) and a==b


def main():
    ap=argparse.ArgumentParser();g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    args=ap.parse_args();value=produce();target=HERE/'results.json'
    if args.write:target.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    elif not equal(json.loads(target.read_text()),value):raise SystemExit('atomic source result differs')
    print(value['verdict'])


if __name__=='__main__':main()
