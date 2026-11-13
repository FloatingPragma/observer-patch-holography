import argparse, hashlib, json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .catalogue import certificate, PARENTS, read_parent
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
    seeds=[str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in ('.py','.md','.tex')]
    result={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}
    for name in PARENTS:
        _,parent=read_parent(name)
        result.update(parent)
        path=f'model/{name}/results.json';result[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    return dict(sorted(result.items()))
def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    before=pins();v=certificate();assert pins()==before;v['source_sha256']=before;path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
    elif not equivalent(json.loads(path.read_text()),v):raise ValueError('mass catalogue replay differs')
    print(v['verdict'])
if __name__=='__main__':main()
