#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys, json
ROOT=Path(__file__).resolve().parent
LOG=ROOT/'outputs'/'latest_validation.log'
cmds=[
 [sys.executable,str(ROOT/'proofs'/'symbolic_completion_proofs.py')],
 [sys.executable,str(ROOT/'producers'/'sm_eft_action_producer.py')],
 [sys.executable,str(ROOT/'checkers'/'check_completion_bundle.py')],
 [sys.executable,str(ROOT/'checkers'/'check_sm_eft_action.py')],
 [sys.executable,str(ROOT/'producers'/'eft_matching_producer.py')],
 [sys.executable,str(ROOT/'checkers'/'check_eft_matching.py')],
 [sys.executable,str(ROOT/'producers'/'rule_engine_a.py')],
 [sys.executable,str(ROOT/'producers'/'rule_engine_b.py')],
 [sys.executable,str(ROOT/'producers'/'enumerate_diagram_universe.py')],
 [sys.executable,str(ROOT/'checkers'/'check_rule_engines.py')],
 [sys.executable,str(ROOT/'producers'/'counterterm_producer.py')],
 [sys.executable,str(ROOT/'checkers'/'check_counterterms.py')],
 [sys.executable,str(ROOT/'producers'/'fj_direct_engine.py')],
 [sys.executable,str(ROOT/'producers'/'fj_scalar_blocks.py')],
 [sys.executable,'-m','pytest','-q',str(ROOT/'tests')],
]
lines=[]
for c in cmds:
    p=subprocess.run(c,cwd=ROOT,text=True,capture_output=True)
    lines.append('$ '+' '.join(c)+'\n'+p.stdout+p.stderr)
    if p.returncode!=0:
        LOG.write_text('\n'.join(lines),encoding='utf-8')
        raise SystemExit(p.returncode)
status='DRAFT_SUFFICIENCY_STACK_DEFINED__SIMULATION_RECEIPTS_OPEN__NO_OPH_NATIVE_POLE_PROMOTION'
lines.append(status+'\n')
LOG.write_text('\n'.join(lines),encoding='utf-8')
print('\n'.join(lines))
