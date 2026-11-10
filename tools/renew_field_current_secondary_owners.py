"""Execute two original downstream owners after the exact documentary replay."""
from __future__ import annotations
import argparse,hashlib,importlib,json
from pathlib import Path
from tools import strict_json
from tools import renew_canonical_live_receipts as live
ROOT=Path(__file__).resolve().parents[1]
AUTHORITY='tools/field_current_secondary_owners.json'
AUTHORITY_SHA256='95ae02d71a43f02b57c5ab961033491312807084ad5d048ae236b02dd3e857e4'
ORDER=('spin_joint_field_limit','protected_defect_catalog')

def execute(name,*,renew=False):
    if live.digest(ROOT/AUTHORITY)!=AUTHORITY_SHA256:raise ValueError('secondary owner authority changed')
    cfg=strict_json.load(ROOT/AUTHORITY)
    if tuple(cfg['owners'])!=ORDER:raise ValueError('secondary owner scope changed')
    live.verify_executable_source_pins(cfg['unchanged_executable_source_sha256'],ROOT)
    spec=cfg['owners'][name]
    if live.digest(ROOT/spec['retained_path'])!=spec['retained_sha256']:raise ValueError('retained downstream source changed')
    old=strict_json.load(ROOT/spec['retained_path']);path=ROOT/spec['receipt'];before=path.read_bytes();current=strict_json.loads(before.decode())
    if live.canonical(live.carry_metadata(spec['receipt'],old,current))!=live.canonical(current):raise ValueError('downstream science changed before replay')
    producer=importlib.import_module(spec['module']);fresh=producer.build()
    admitted=live.carry_metadata(spec['receipt'],old,fresh)
    if live.canonical(admitted)!=live.canonical(fresh):raise ValueError('downstream scientific replay changed')
    for _,relative,h,p in live.bindings(ROOT,spec['receipt'],fresh):
        if live.digest(p)!=h:raise ValueError('downstream source changed: '+relative)
    if name=='protected_defect_catalog' and producer.render(fresh)!=(producer.HERE/'CATALOG.md').read_text():raise ValueError('catalogue scientific prose changed')
    live.verify_executable_source_pins(cfg['unchanged_executable_source_sha256'],ROOT)
    if path.read_bytes()!=before:raise ValueError('downstream receipt changed during replay')
    if renew:path.write_text(producer.canonical(admitted))
    elif live.canonical(current)!=live.canonical(admitted):raise ValueError('downstream custody stale')
    return {'receipt':spec['receipt'],'current_sha256':live.digest(path),'retained_path':spec['retained_path'],'retained_sha256':spec['retained_sha256'],'actual_unchanged_producer_executed':spec['module']+':build','canonical_scientific_payload_exactly_preserved':True,'scientific_statuses_changed':False}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--renew',action='store_true');a=p.parse_args()
    print(json.dumps({n:execute(n,renew=a.renew)for n in ORDER},indent=2))
if __name__=='__main__':main()
