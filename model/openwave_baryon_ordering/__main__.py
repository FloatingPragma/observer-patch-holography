import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from .calculation import certificate
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
 seeds=[str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md','.tex')]
 seeds+=list(SOURCE_FILES)+['model/lepton_cycle/experimental_inputs.json','model/openwave_vectorlike_flavour/results.json','model/openwave_native_composites/results.json']
 return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();before=pins();r=certificate();assert pins()==before;r['source_sha256']=before;path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:assert json.loads(path.read_text())==r
 print(r['verdict'])
if __name__=='__main__':main()
