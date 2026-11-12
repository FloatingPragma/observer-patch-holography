import argparse,json,hashlib
from pathlib import Path
from .review import review,TARGETS,ROOT
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_larmor_dynamics.certificate import equivalent
HERE=Path(__file__).resolve().parent

def pins():
 seeds=set(SOURCE_FILES)
 for p in HERE.iterdir():
  if p.suffix in('.py','.md'):seeds.add(str(p.relative_to(ROOT)))
 for name in TARGETS:
  path=ROOT/'model'/name/'results.json';d=json.loads(path.read_text());seeds.add(str(path.relative_to(ROOT)));seeds.update(d['source_sha256'])
 return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}
def main():
 parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=parser.parse_args();before=pins();r=review();assert pins()==before;r['source_sha256']=before;p=HERE/'results.json'
 if a.write:p.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:assert equivalent(r,json.loads(p.read_text()))
 print(r['verdict'])
if __name__=='__main__':main()
