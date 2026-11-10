"""Run four unchanged downstream owners with authenticated custody-only inputs.

The original canonical numbers and complete NPZ are retained. Each actual owner
checker recomputes its scientific output and applies its existing exact or
numerical admission; that admission is not a physical error bound.
"""
from __future__ import annotations
import argparse,copy,hashlib,importlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tools import strict_json
from tools import renew_canonical_live_receipts as live
AUTHORITY='tools/source_record_readout_owners.json'
AUTHORITY_SHA256='5ce926726f445fba143bd4b18936596fa72f9e12481185031e58184bb73b87dc'
OWNERS=('particle_identification','protected_particle_readout','resolved_particle_spectrum','conditional_standard_model')


def require(condition,message):
    if not condition:raise ValueError(message)


def authority(root=ROOT):
    root=Path(root);require(live.digest(root/AUTHORITY)==AUTHORITY_SHA256,'readout owner authority changed')
    value=strict_json.load(root/AUTHORITY)
    require(set(value['owners'])==set(OWNERS),'readout owner inventory changed')
    live.verify_executable_source_pins(value['unchanged_executable_source_sha256'],root)
    return value


def project(spec,old,pins,*,root=ROOT):
    field=spec['source_field'];original=old[field]
    require(type(pins) is dict and set(pins)==set(original),'readout source scope changed')
    changed={p for p in original if pins[p]!=original[p]}
    require(changed<=set(spec['only_renewable_source_paths']),'unreviewed readout source changed')
    for name,expected in pins.items():require(live.digest(Path(root)/name)==expected,'readout source pin differs: '+name)
    result=copy.deepcopy(old);result[field]=pins
    return result


def execute(owner,*,renew=False):
    approved=authority();require(owner in OWNERS,'unknown readout owner');spec=approved['owners'][owner]
    require(live.digest(ROOT/spec['retained_path'])==spec['retained_sha256'],'retained readout changed')
    old=strict_json.load(ROOT/spec['retained_path']);output=ROOT/spec['receipt'];before=output.read_bytes()
    current=strict_json.loads(before.decode());field=spec['source_field']
    # Authenticate canonical science even when its live pins are superseded.
    comparison=copy.deepcopy(current);comparison[field]=old[field]
    require(live.canonical(comparison)==live.canonical(old),'readout science changed before replay')
    producer=importlib.import_module(spec['module']);raw=getattr(producer,spec['source_entry'])()
    pins=raw[field] if spec['source_entry']=='run' else raw
    candidate=project(spec,old,pins)
    array_hash=None
    with tempfile.TemporaryDirectory(prefix='oph-current-readout-owner-') as folder:
        path=Path(folder)/output.name;path.write_text(json.dumps(candidate,indent=2,sort_keys=True)+'\n')
        if 'complete_record_archive' in old:
            name=old['complete_record_archive']['filename'];require(Path(name).name==name,'array path is not adjacent')
            source=output.parent/name;array_hash=live.digest(source)
            require(array_hash==old['complete_record_archive']['sha256'],'retained complete record array changed')
            shutil.copyfile(source,path.parent/name)
        run=subprocess.run([sys.executable,'-B','-m',spec['module'],'--check',str(path)],cwd=ROOT,capture_output=True,text=True)
        require(run.returncode==0,'unchanged owner checker failed:\n'+run.stdout+run.stderr)
    authority();require(output.read_bytes()==before,'readout changed during replay')
    if renew:output.write_text(json.dumps(candidate,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:require(live.canonical(current)==live.canonical(candidate),'readout custody stale; run --renew')
    return {'owner':owner,'receipt':spec['receipt'],'retained_sha256':spec['retained_sha256'],
            'current_sha256':live.digest(output),'actual_unchanged_owner_checker_executed':True,
            'canonical_scientific_payload_exactly_preserved':True,'unchanged_complete_array_sha256':array_hash,
            'fresh_binary64_science_byte_equality_claimed':False,'physical_error_bound_from_replay_admission':False,
            'scientific_statuses_changed':False}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--owner',choices=OWNERS,action='append')
    g=p.add_mutually_exclusive_group(required=True);g.add_argument('--renew',action='store_true');g.add_argument('--check',action='store_true')
    a=p.parse_args();selected=set(a.owner or OWNERS)
    for name in OWNERS:
        if name in selected:print(json.dumps(execute(name,renew=a.renew),sort_keys=True),flush=True)


if __name__=='__main__':main()
