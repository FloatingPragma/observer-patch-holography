import argparse,json,hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .analysis import run
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
 seeds=[str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')]
 seeds+=['model/openwave_scale_coupling/PROOF.md','model/openwave_scale_coupling/results.json']
 files=import_closure(seeds)|{str(p.relative_to(ROOT))for folder in[HERE/'data',HERE.parent/'openwave_scale_coupling/data']for p in folder.glob('*.npz')}
 return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(files)}
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();before=pins();v=run();assert pins()==before;v['source_sha256']=before;path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
 elif json.loads(path.read_text())!=v:raise ValueError('larger-volume replay differs')
 print(v['verdict'])
if __name__=='__main__':main()
