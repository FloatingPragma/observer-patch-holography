import argparse,hashlib,json,math
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .spectrum import certificate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
 seeds=[str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md','.tex')]
 from model.common_source.chiral_candidate.action import SOURCE_FILES
 seeds+=list(SOURCE_FILES)
 return bound_pins({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))})
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();before=pins();r=certificate();assert before==pins();r['source_sha256']=bound_pins(before);path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:
  old=json.loads(path.read_text())
  # Validate both retained and recomputed residuals before omitting them from
  # the equality comparison of deterministic rounded observer records.
  def residual_bounds(x):
   if isinstance(x,dict):
    for k,v in x.items():
     if k.endswith('residual'):assert isinstance(v,(int,float))and math.isfinite(v)and 0<=v<1e-10,'invalid composite residual'
     else:residual_bounds(v)
   elif isinstance(x,list):
    for v in x:residual_bounds(v)
  residual_bounds(old);residual_bounds(r)
  def clean(x):
   if isinstance(x,dict):return {k:clean(v) for k,v in x.items() if not k.endswith('residual')}
   if isinstance(x,list):return [clean(v) for v in x]
   return x
  assert clean(old)==clean(r),'composite replay differs'
 print(r['verdict'])
if __name__=='__main__':main()
