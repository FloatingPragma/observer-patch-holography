import argparse,hashlib,json
from pathlib import Path
from .analysis import run
from model.particle_common_contract.provenance import import_closure
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def produce():
    r=run();paths=set(HERE.glob('*.py'))|set(HERE.glob('*.md'))|set(HERE.glob('*.txt'))|set((HERE/'data').glob('*.npz'))
    paths.update(ROOT/p for p in import_closure([str((HERE/'analysis.py').relative_to(ROOT))]))
    r['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    return r


def equal(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    if isinstance(b,float):return isinstance(a,(int,float)) and abs(a-b)<=1e-10*max(1,abs(b))
    return type(a) is type(b) and a==b


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();r=produce();target=HERE/'results.json'
    if args.write:target.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    elif not equal(json.loads(target.read_text()),r):raise SystemExit('scale-coupling source result differs')
    print(r['verdict'])


if __name__=='__main__':main()
