from pathlib import Path
import json,os,subprocess,sys
from model.openwave_m10_review_gate.source import authenticate,parse,CONTRACT
from .headless import run as native_run

ROOT=Path(__file__).resolve().parents[2]
def run(openwave_root):
    upstream=authenticate(openwave_root)
    working=upstream/'openwave/xperiments/m10_oph'
    command=[sys.executable,'-m','model.openwave_native_platform_gate.headless']
    completed=subprocess.run(command,cwd=working,
        env={**os.environ,'PYTHONPATH':str(ROOT),'GIT_OPTIONAL_LOCKS':'0'},capture_output=True,text=True,check=True)
    actual=json.loads(completed.stdout);assert actual==native_run()
    status=parse((upstream/'MODELS.md').read_text());authenticate(upstream)
    return {'schema':'oph.native_platform_gate.v1',
        'verdict':'PASS_AUTHOR_HEADLESS_NATIVE_REPLAY_WITH_PINNED_PLATFORM_GATE',
        'upstream_commit':CONTRACT['upstream_commit'],'upstream_sha256':CONTRACT['upstream_sha256'],
        'official_status':status,'native_execution':actual,
        'execution':{'command':'python -m model.openwave_native_platform_gate.headless',
            'cwd':'<pinned-openwave>/openwave/xperiments/m10_oph',
            'source_import':'PYTHONPATH=<RER>; source code executes, no result receipt is imported',
            'exit_code':completed.returncode,'owner':'author-side local platform-layout replay',
            'maintainer_run_or_upstream_acceptance':False},
        'criterion_gate':CONTRACT['gate'],
        'scope':{'native_complete_channel_actually_executed':True,
            'criterion_specific_physical_comparison_completed':False,
            'new_official_awards':0,'upstream_files_changed':False,'upstream_submission_sent':False},
        'next_acceptance':'Join a specific measured quantity to this process with an explicit preparation, complete observer instrument and physical clock/unit map; test the literal criterion and seek maintainer review separately.'}
