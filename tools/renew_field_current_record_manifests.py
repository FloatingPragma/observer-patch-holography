"""Renew seven declared artifact slots after actual live owner replays.

Every prior manifest and canonical receipt is retained. The exact seven-path
scope, owner replay record, executable ancestry and scientific payload equality
are checked; other manifest values and artifact names cannot change.
"""
from __future__ import annotations
import argparse,copy,hashlib,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tools import strict_json
from tools import renew_canonical_live_receipts as live

AUTHORITY='tools/field_current_record_manifest_authority.json'
AUTHORITY_SHA256='7d1bbcd60e8b8f5b311d10ca63c7e7b45b20d0903a94ecdbc4a8b1d9f0842f03'
ARCHIVE='model/completion/historical/field_current_records_inputs/'
EXPECTED_TARGETS={
 'model/common_source/chiral_candidate/results.json',
 'model/common_source/full_repair_response/results.json',
 'model/common_source/hierarchy_diagnosis/results.json',
 'model/common_source/native_defect_decision_live_results.json',
 'model/common_source/native_defect_decision_results.json',
 'model/common_source/native_repair_mass_results.json',
 'model/empirical_diagnostics/cross_scale_results.json'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition,message):
    if not condition:raise ValueError(message)


def project(name,old,approved,*,root=ROOT):
    result=copy.deepcopy(old)
    for row in approved['only_live_artifact_entries']:
        if row['manifest']!=name:continue
        key=row['artifact_key'];target=row['target']
        require(result['files'][key]==row['retained_pointer_sha256'],'retained manifest target changed')
        require(str(Path(name).parent/key)==target,'manifest target path changed')
        require(sha(Path(root)/target)==row['executed_current_sha256'],'executed canonical artifact changed')
        result['files'][key]=row['executed_current_sha256']
    return result


def execute(*,renew=False,root=ROOT):
    root=Path(root)
    require(sha(root/AUTHORITY)==AUTHORITY_SHA256,'manifest transition authority changed')
    approved=strict_json.load(root/AUTHORITY)
    require({r['target'] for r in approved['only_live_artifact_entries']}==EXPECTED_TARGETS
            and len(approved['only_live_artifact_entries'])==7,'manifest transition scope changed')
    replay=approved['executed_owner_replays']
    require(sha(root/replay['path'])==replay['sha256'],'actual owner execution record changed')
    live.load_catalog()  # Authenticate all unchanged scientific executable ancestry.
    for row in approved['only_live_artifact_entries']:
        target=row['target']
        require(sha(root/(ARCHIVE+target))==approved['retained_canonical_receipt_sha256'][target],
                'retained canonical receipt changed: '+target)
        old=strict_json.load(root/(ARCHIVE+target));current=strict_json.load(root/target)
        require(live.canonical(live.carry_metadata(target,old,current,root=root))==live.canonical(current),
                'canonical scientific artifact changed: '+target)
        for _,relative,expected,path in live.bindings(root,target,current):
            require(sha(path)==expected,'canonical source custody changed: '+relative)
    outputs={}
    for name,row in approved['retained_manifests'].items():
        require(sha(root/row['path'])==row['sha256'],'retained manifest bytes changed')
        old=strict_json.load(root/row['path']);expected=project(name,old,approved,root=root)
        current=strict_json.load(root/name)
        require(current==expected or (renew and current==old),'unlisted manifest field changed')
        outputs[name]=expected
    if renew:
        for name,value in outputs.items():(root/name).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    return {'schema':'oph.field_current_record_manifest_validation.v1','exact_artifact_slots':7,
            'scientific_statuses_changed':False,'actual_replay_record_sha256':replay['sha256'],
            'current_manifests':{p:sha(root/p) for p in outputs}}


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--renew',action='store_true');g.add_argument('--check',action='store_true')
    a=p.parse_args();print(json.dumps(execute(renew=a.renew),indent=2,sort_keys=True))


if __name__=='__main__':main()
