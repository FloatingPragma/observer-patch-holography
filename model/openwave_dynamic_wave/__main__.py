import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .certificate import run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
 seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')};seeds.add('model/requirements.txt')
 return bound_pins({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))})
def equivalent(a,b):
 if isinstance(b,dict):return isinstance(a,dict)and set(a)==set(b)and all(equivalent(a[k],v)for k,v in b.items())
 if isinstance(b,list):return isinstance(a,list)and len(a)==len(b)and all(equivalent(x,y)for x,y in zip(a,b))
 if isinstance(b,float):return isinstance(a,(float,int))and abs(a-b)<=max(2e-12,abs(b)*2e-10)
 return type(a)is type(b)and a==b
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();before=pins();r=json.loads(json.dumps(run()));assert before==pins();r['source_sha256']=bound_pins(before);path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:assert equivalent(json.loads(path.read_text()),r),'dynamic source replay differs'
 print(r['verdict'])
if __name__=='__main__':main()
