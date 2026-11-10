"""Execute the fixed internal parent closure of the particle completion owner."""
from __future__ import annotations
import hashlib,importlib,json
from pathlib import Path
from tools import strict_json
from tools import renew_canonical_live_receipts as live
ROOT=Path(__file__).resolve().parents[1]
AUTHORITY='tools/field_current_interface_parents.json'
AUTHORITY_SHA256='a3cbc7f435f34e423dec0420b88fdd83db3a66c734bafcaaa99bca0dcee3e369'
RECORD='model/completion/field_current_interface_parent_replays.json'

def execute(*,renew=False):
    if live.digest(ROOT/AUTHORITY)!=AUTHORITY_SHA256:raise ValueError('internal parent authority changed')
    cfg=strict_json.load(ROOT/AUTHORITY)
    live.verify_executable_source_pins(cfg['unchanged_executable_source_sha256'],ROOT)
    receipts={row['receipt']:name for name,row in cfg['owners'].items()}
    done={};pending=set(cfg['owners'])
    while pending:
        ready=[]
        for name in sorted(pending):
            spec=cfg['owners'][name];old=strict_json.load(ROOT/spec['retained_path'])
            dependencies={receipts[p]for _,p,_,_ in live.bindings(ROOT,spec['receipt'],old)if p in receipts and receipts[p]!=name}
            if dependencies<=set(done):ready.append(name)
        if not ready:raise ValueError('internal source dependency cycle')
        for name in ready:
            spec=cfg['owners'][name];path=ROOT/spec['receipt'];before=path.read_bytes()
            if live.digest(ROOT/spec['retained_path'])!=spec['retained_sha256']:raise ValueError('retained internal parent changed')
            old=strict_json.load(ROOT/spec['retained_path']);current=strict_json.loads(before.decode())
            if live.canonical(live.carry_metadata(spec['receipt'],old,current))!=live.canonical(current):raise ValueError('internal parent science changed before execution')
            producer=importlib.import_module(spec['module']);fresh=producer.produce()
            admitted=live.carry_metadata(spec['receipt'],old,fresh)
            if live.canonical(admitted)!=live.canonical(fresh):raise ValueError('internal original producer science changed: '+name)
            for _,p,h,q in live.bindings(ROOT,spec['receipt'],fresh):
                if live.digest(q)!=h:raise ValueError('internal owner input changed: '+p)
            if path.read_bytes()!=before:raise ValueError('internal parent changed during execution')
            if renew:path.write_text(producer.canonical(admitted))
            elif live.canonical(current)!=live.canonical(admitted):raise ValueError('internal parent metadata stale')
            done[name]={'receipt':spec['receipt'],'current_sha256':live.digest(path),'retained_sha256':spec['retained_sha256'],'actual_unchanged_producer_executed':spec['module']+':produce','scientific_payload_exactly_preserved':True,'source_inventory_exactly_preserved':True}
            pending.remove(name);print('PASS internal parent '+name,flush=True)
    live.verify_executable_source_pins(cfg['unchanged_executable_source_sha256'],ROOT)
    value={'schema':'oph.field_current_internal_parent_executions.v1','authority_sha256':AUTHORITY_SHA256,'executed':done,'physical_statuses_changed':False}
    if renew:(ROOT/RECORD).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    return value
