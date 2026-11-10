"""Replay the two existing class consumers after the class custody seal.

Their receipts bind the current class. The class does not bind them back. This
ordering avoids a circular content hash while preserving both original owner
checks, all scientific payloads and the complete current source maps.
"""
from __future__ import annotations
import argparse,hashlib,json,subprocess,sys,tempfile
from pathlib import Path
from tools import strict_json
from tools import renew_canonical_live_receipts as live
ROOT=Path(__file__).resolve().parents[1]
ARCH='model/completion/historical/field_current_records_inputs/'
RECORD='model/completion/field_current_class_consumer_replays.json'
RETAINED={'model/particle_identification/completion_results.json': 'fc32cb97292b90edb3e333a39551740217efc71dc0e124409f168f654aa8ae79', 'model/regional_clock/results.json': '6fe58a8991dd9fb54b0723e76b7f428d84c14250e563b030b2570fc66f4f734f'}

def execute(*,renew=False):
    from model.completion import field_current_records_custody as custody
    custody.validate()
    klass=live.digest(ROOT/custody.CLASS)
    for name,h in RETAINED.items():
        if live.digest(ROOT/(ARCH+name))!=h:raise ValueError('retained class consumer changed')
    report={}
    from tools.renew_field_current_interface_parents import execute as replay_internal
    report['internal_parents']=replay_internal(renew=renew)
    report['particle_completion']=live.execute('particle_completion',renew=renew)
    name='model/regional_clock/results.json';path=ROOT/name;before=path.read_bytes()
    old=strict_json.load(ROOT/(ARCH+name));current=strict_json.loads(before.decode())
    if live.canonical(live.carry_metadata(name,old,current))!=live.canonical(current):raise ValueError('regional science changed before replay')
    from model.regional_clock import verify as owner
    fresh=owner.report();candidate=live.carry_metadata(name,old,fresh)
    with tempfile.TemporaryDirectory(prefix='oph-current-regional-owner-')as folder:
        p=Path(folder)/'candidate.json';p.write_text(json.dumps(candidate,indent=2,sort_keys=True)+'\n')
        done=subprocess.run([sys.executable,'-B','-W','error','-m','model.regional_clock.verify','--check',str(p)],cwd=ROOT,capture_output=True,text=True)
        if done.returncode:raise ValueError('unchanged regional owner rejected:\n'+done.stdout+done.stderr)
    if path.read_bytes()!=before or live.digest(ROOT/custody.CLASS)!=klass:raise ValueError('class consumer inputs changed during replay')
    if renew:path.write_text(json.dumps(candidate,indent=2,sort_keys=True)+'\n')
    elif live.canonical(current)!=live.canonical(candidate):raise ValueError('regional source custody stale')
    report['regional_clock']={'receipt':name,'current_sha256':live.digest(path),'retained_sha256':RETAINED[name],'actual_unchanged_owner_checker_executed':True,'scientific_payload_exactly_preserved':True,'numerical_diagnostics_admitted_only_by_original_owner':True}
    custody.validate()
    value={'schema':'oph.field_current_class_consumer_replays.v1','current_class_sha256':klass,'executed':report,'old_class_review_reinterpreted':False,'class_does_not_bind_consumers_back':True}
    if renew:(ROOT/RECORD).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:
        old_record=strict_json.load(ROOT/RECORD)
        if old_record!=value:raise ValueError('class consumer execution record differs')
    return value

def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True);g.add_argument('--renew',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();print(json.dumps(execute(renew=a.renew),indent=2))
if __name__=='__main__':main()
