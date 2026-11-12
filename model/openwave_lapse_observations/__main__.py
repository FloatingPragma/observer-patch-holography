import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .observations import run
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def bindings():
    files={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')}
    files|={'model/openwave_energy_lapse/PROOF.md','model/openwave_energy_lapse/results.json','model/openwave_metric_observables/PROOF.md','model/matched_einstein_limit/PROOF.md'}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(files))}
def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();before=bindings();v=run();assert bindings()==before;v['source_sha256']=before;path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
    elif json.loads(path.read_text())!=v:raise ValueError('lapse observer replay differs')
    print(v['verdict'])
if __name__=='__main__':main()
