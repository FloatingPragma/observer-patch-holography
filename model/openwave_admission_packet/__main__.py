"""Execute three actual OpenWave overlays without altering tracked upstream files."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from model.openwave_platform_kernel.__main__ import UPSTREAM
from model.openwave_platform_kernel.export import source_packet
from model.openwave_pure_su2_source.custody import ROOT,pin_package,equivalent

HERE=Path(__file__).resolve().parent
KERNEL=ROOT/'model/openwave_platform_kernel'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def pins():
    return pin_package(HERE, ['model/openwave_platform_kernel/engine.py',
        'model/openwave_platform_kernel/export.py',
        'model/openwave_platform_kernel/source_packet.json',
        'model/openwave_admission_packet/prediction_audit.json'])


def run(upstream, python):
    env={**os.environ,'GIT_OPTIONAL_LOCKS':'0','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=upstream,env=env,text=True).strip()
    assert head==UPSTREAM
    assert not subprocess.check_output(['git','diff','--name-only','HEAD','--'],cwd=upstream,env=env).strip()
    original=KERNEL/'source_packet.json'
    assert json.loads(original.read_text())==source_packet()
    before=pins()
    overlay=upstream/'openwave/xperiments/oph_admission'
    overlay.mkdir(exist_ok=True)
    files={'__init__.py':HERE/'__init__.py','platform_tests.py':HERE/'platform_tests.py',
        'engine.py':KERNEL/'engine.py','source_packet.json':original}
    assert {p.name for p in overlay.iterdir()} <= set(files)|{'__pycache__'}
    for name,source in files.items():
        target=overlay/name
        if target.exists():assert target.read_bytes()==source.read_bytes(),('different existing overlay',name)
        else:shutil.copyfile(source,target)
    outputs=[]
    with tempfile.TemporaryDirectory(prefix='oph-admission-replays-') as temporary:
        for test in ('stability','rest_spectrum','mass_shell'):
            output=Path(temporary)/(test+'.json')
            command=[python,'-m','openwave.xperiments.oph_admission.platform_tests',
                '--source','openwave/xperiments/oph_admission/source_packet.json',
                '--test',test,'--output',str(output)]
            result=subprocess.run(command,cwd=upstream,env={**env,'PYTHONPATH':str(upstream)},check=True)
            value=json.loads(output.read_text())
            assert value['upstream_version']=='26.12'
            assert value['source_packet_sha256']==sha(original)
            outputs.append({'command':'python -m openwave.xperiments.oph_admission.platform_tests --source openwave/xperiments/oph_admission/source_packet.json --test '+test+' --output <receipt.json>',
                'actual_exit_status':result.returncode,'result':value})
    assert pins()==before
    assert not subprocess.check_output(['git','diff','--name-only','HEAD','--'],cwd=upstream,env=env).strip()
    return {'schema':'oph.openwave_admission_replay_packet.v1',
        'verdict':'PASS_THREE_EXECUTED_FINITE_PLATFORM_REPLAYS',
        'source_sha256':before,'upstream_commit':head,
        'unchanged_tracked_upstream_source':True,
        'three_actual_platform_replays':outputs,
        'scope':{'same_complete_180_coordinate_native_free_source':True,
            'literal_rows_4_11_30_fully_validated':False,
            'genuine_unexposed_prediction_demonstrated':False,
            'official_admission_or_award':False,
            'new_independent_nature_matches':0},
        'admission_gate':'The phase2/9 analyzer and electron calibration have retrospective target exposure. Stable finite free evolution and sine-dispersion checks do not establish a novel observational prediction or a physical interacting continuum.'}


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--openwave-root',type=Path,required=True)
    p.add_argument('--python',default=sys.executable)
    a=p.parse_args();value=run(a.openwave_root.resolve(),a.python)
    path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:assert equivalent(json.loads(path.read_text()),value)
    print(value['verdict'])


if __name__=='__main__':main()
