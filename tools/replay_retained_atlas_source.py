#!/usr/bin/env python3
"""Execute unchanged primary owners in their exact historical input view.

All files are copied into an isolated temporary repository, with six paths
resolved to their authenticated historical receipt bytes. No module or guard
is monkeypatched, and no canonical primary receipt is rewritten.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.particle_atlas_explainer import retained_source_equivalence as admission  # noqa: E402

SCRIPT = '''
import hashlib,importlib,json,platform
from pathlib import Path
approved=json.loads(Path('replay_authority.json').read_text())
receipts={}
for name,spec in approved['retained_receipts'].items():
    module,function=spec['producer'].split(':')
    print('Replaying '+module,flush=True)
    producer=importlib.import_module(module)
    fresh=getattr(producer,function)()
    old=json.loads(Path(name).read_text())
    if fresh!=old: raise ValueError('historical scientific replay differs: '+name)
    receipts[name]={'producer':spec['producer'],'retained_sha256':spec['sha256'],
                    'complete_semantic_json_exactly_equal':True}
from model.particle_common_contract import contract
from hashlib import sha256
members=[]
bound_contract=json.loads(Path('model/particle_common_contract/contract.json').read_text())
for level in (0,1):
    member=contract.member(level);contract.verify_member(member)
    key='fixed_primary_member' if level==0 else 'next_declared_member'
    if member!=bound_contract[key]:raise ValueError('retained member identity changed')
    physical=json.loads(json.dumps(member))
    del physical['source']['implementation_sha256']
    for section in ('state','reference','clock'):del physical[section]['source_id']
    del physical['reference']['state_id'];del physical['ids'];del physical['member_id']
    members.append({'level':level,'member_id':member['member_id'],**member['ids'],
                    'physical_specification_sha256':contract.digest(physical)})
for name,expected in approved['retained_input_sha256'].items():
    if sha256(Path(name).read_bytes()).hexdigest()!=expected:
        raise ValueError('historical input changed during replay: '+name)
Path('replay_result.json').write_text(json.dumps({'receipts':receipts,'members':members,
     'python':platform.python_version(),'all_retained_receipts_exactly_replayed':True},indent=2)+'\\n')
'''


def run():
    approved = admission.load_authority()
    inputs = admission.historical_input_paths(approved=approved)
    with tempfile.TemporaryDirectory(prefix='oph-retained-primary-') as directory:
        view = Path(directory)
        for name, expected in approved['retained_input_sha256'].items():
            original = inputs[name]
            target = view / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original, target)
            if admission.sha(target) != expected:
                raise ValueError('historical view does not match retained bytes: ' + name)
        (view / 'replay_authority.json').write_text(json.dumps(approved))
        env = dict(os.environ, PYTHONPATH=str(view), PYTHONDONTWRITEBYTECODE='1')
        subprocess.run([sys.executable, '-c', SCRIPT], cwd=view, env=env, check=True)
        result = json.loads((view / 'replay_result.json').read_text())
    admission.historical_input_paths(approved=approved)
    result.update({'schema': 'oph.retained_primary_replay.v1',
                   'authority_sha256': admission.AUTHORITY_SHA256,
                   'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'historical_view_guard_unchanged': True,
                   'canonical_primary_receipts_written': False})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', required=True)
    parser.parse_args()
    result = run()
    (ROOT / admission.REPLAY).write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print('PASS: unchanged historical owners and original physical IDs; separate current navigation')


if __name__ == '__main__':
    main()
