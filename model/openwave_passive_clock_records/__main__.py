"""Produce and independently replay the declared native rest-clock records."""
import argparse,hashlib,json
from pathlib import Path
from .clock import certificate
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def pins():
    parent='model/openwave_common_electromagnetism/results.json'
    paths=set(json.loads((ROOT/parent).read_text())['source_sha256'])|{parent}
    paths.add('model/openwave_lepton_widths/inputs.json')
    paths|={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in ('.py','.md','.tex')}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(paths))}

def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    before=pins();r=certificate();assert before==pins();r['source_sha256']=before
    path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    else:assert equivalent(json.loads(path.read_text()),r)
    print(r['verdict'])
if __name__=='__main__':main()
