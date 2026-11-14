"""Copy a verified source closure. Re-freezing is an explicit research action."""
from pathlib import Path
import hashlib, importlib.util, json, shutil, sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
PINS={
 'openwave-tuning/next/CHAIN2/local12_limit_audit/MANIFEST.json':'997916c3b216f3006b9c88b8ef39419a213df5075e8159f47298d87ad295e2ad',
 'openwave-tuning/common_dynamics/clock/MANIFEST.json':'6f80b0fe899cdeffd59d3742c318fe627817c19685d8f164abed69721f0b0992',
 'openwave-tuning/common_dynamics/continuum/MANIFEST.json':'3067f7bdea43ca5137d6b115f7e1bd92b4d753ba99684ed198ba2da170240adf'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def freeze():
 paths=set(); queue=[ROOT/p for p in PINS]
 queue += [ROOT/'openwave-tuning/auditor/native_walk_bridge/COMMON_SOURCE_CUSTODY.json',ROOT/'openwave-tuning/auditor/common_walk_particle/MANIFEST.json']
 for p,v in PINS.items():assert sha(ROOT/p)==v,p
 while queue:
  p=queue.pop()
  if p in paths:continue
  paths.add(p)
  if p.suffix!='.json' or not (p.name.endswith('MANIFEST.json') or p.name in ('COMMON_SOURCE_CUSTODY.json','COMMON_WALK_RECEIPT.json','SOURCE_ERROR.json')):continue
  d=json.loads(p.read_text())
  for k in ('files_sha256','source_sha256'):
   for rel,v in d.get(k,{}).items():
    candidates=[ROOT/rel,ROOT/'reverse-engineering-reality'/rel,p.parent/rel]
    candidate=next((x for x in candidates if x.is_file()),candidates[-1])
    assert sha(candidate)==v,str(candidate)
    queue.append(candidate)
 sys.path[:0]=[str(ROOT/'openwave-tuning/next/CHAIN2/local12_limit_audit'),str(ROOT/'openwave-tuning/common_dynamics/continuum'),str(ROOT/'openwave-tuning/common_dynamics/particles')]
 import refinement_adapter as a, continuum, mass_coin
 from model.native_interface import pointer
 from model.native_clock import configuration
 a.source_constants(48); pointer.exact_instrument_audit()
 for mod in list(sys.modules.values()):
  f=getattr(mod,'__file__',None)
  if f:
   p=Path(f).resolve()
   if p.is_file() and p.is_relative_to(ROOT) and not p.is_relative_to(HERE):paths.add(p)
 for rel in ('openwave-tuning/common_dynamics/particles/mass_coin.py','openwave-tuning/common_dynamics/particles/MASS_COIN.json','openwave-tuning/auditor/common_walk_review/source_error.py','openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json'):
  paths.add(ROOT/rel)
 for mod in (pointer,configuration):paths.update(ROOT/'reverse-engineering-reality'/p for p in mod.SOURCE_PATHS)
 # Package initializers preserve the imported module semantics.
 for p in list(paths):
  for parent in p.parents:
   if parent==ROOT:break
   init=parent/'__init__.py'
   if init.is_file():paths.add(init)
 pins={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}
 for rel,v in pins.items():
  dest=HERE/'inputs'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,dest);assert sha(dest)==v
 assert all(sha(ROOT/p)==v for p,v in pins.items())
 (HERE/'INPUTS.json').write_text(json.dumps({'frozen_manifests':PINS,'files_sha256':pins},indent=2)+'\n')
 print('frozen',len(pins),'files')
if __name__=='__main__':freeze()
