#!/usr/bin/env python3
"""Reproduce frozen review checks, retaining stdout/stderr even on failure."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
NAMES=('20261115_codex_intrinsic_phase_process','20261115_codex_retained_family_oracle')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    assert not sys.flags.optimize
    prereg=json.loads((HERE/'PREREG.json').read_text())
    for path,value in prereg['review_files'].items():
        assert sha(HERE/path)==value, f'review freeze mismatch: {path}'
    args.output.mkdir(parents=True,exist_ok=False)
    env={k:v for k,v in os.environ.items() if not k.startswith('PYTHON')}
    runs=[]
    start=time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='row11_completion_') as directory:
        staged=Path(directory)
        for name in NAMES:
            packet=HERE/'paused_originals'/name
            original=ROOT/'openwave/row_11_lepton_masses/work'/name
            # Review reads preserved originals, not concurrently changing docs.
            sources=json.loads((packet/'SOURCES.json').read_text())
            for relative,value in sources.items():
                frozen=HERE/'frozen_sources'/relative
                source=frozen if frozen.exists() else ROOT/relative
                assert sha(source)==value, f'missing original bytes: {relative}'
                target=staged/relative
                target.parent.mkdir(parents=True,exist_ok=True)
                if target.exists():
                    assert sha(target)==value
                else:
                    shutil.copyfile(source,target)
            target=staged/original.relative_to(ROOT)
            shutil.copytree(packet,target,dirs_exist_ok=True)
        supplementary=json.loads((HERE/'SUPPLEMENTARY_SOURCES.json').read_text())
        for relative,value in supplementary.items():
            source=ROOT/relative
            assert sha(source)==value
            target=staged/relative
            target.parent.mkdir(parents=True,exist_ok=True)
            if target.exists():
                assert sha(target)==value
            else:
                shutil.copyfile(source,target)
        jobs=[('phase_independent',[sys.executable,'-I','-B',str(HERE/'phase_check.py')])]
        for name in NAMES:
            packet=staged/'openwave/row_11_lepton_masses/work'/name
            jobs.append((name+'_replay',[sys.executable,'-I','-B',str(HERE/'replay_child.py'),
                         '--root',str(staged),'--packet',str(packet),
                         '--extra',str(HERE/'SUPPLEMENTARY_SOURCES.json')]))
        oracle=staged/'openwave/row_11_lepton_masses/work'/NAMES[1]
        jobs.append(('oracle_independent',[sys.executable,'-I','-B',str(HERE/'oracle_check.py'),
                     '--packet',str(oracle)]))
        for label,command in jobs:
            began=datetime.now(timezone.utc).isoformat()
            tick=time.perf_counter()
            process=subprocess.run(command,cwd=staged,env=env,text=True,capture_output=True)
            (args.output/(label+'.stdout')).write_text(process.stdout)
            (args.output/(label+'.stderr')).write_text(process.stderr)
            run={'job':label,'started_at_utc':began,'command':command,
                 'exit_code':process.returncode,'seconds':time.perf_counter()-tick,
                 'stdout_sha256':sha(args.output/(label+'.stdout')),
                 'stderr_sha256':sha(args.output/(label+'.stderr'))}
            if process.returncode==0:
                result=json.loads(process.stdout)
                run['all_pass']=result['all_pass']
            runs.append(run)
            (args.output/'RUNS.json').write_text(json.dumps(runs,indent=2,sort_keys=True)+'\n')
            print(json.dumps(run),flush=True)
            if process.returncode:
                raise SystemExit(f'{label} failed; original logs preserved in {args.output}')
    summary={'schema':'row11.completion_run.v1','all_pass':all(r['all_pass'] for r in runs),
             'prereg_sha256':sha(HERE/'PREREG.json'),'runs':runs,
             'seconds':time.perf_counter()-start}
    (args.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    main()
