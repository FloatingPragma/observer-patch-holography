"""Local reviewer harness; executes only the two already-read checker scripts."""
from pathlib import Path
import copy
import json
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
base_r=json.loads((HERE/'RECEIPT.json').read_text())
base_d=json.loads((HERE/'RECORDS.json').read_text())
base_t=json.loads((HERE/'TRAJECTORIES.json').read_text())
base_tr=json.loads((HERE/'TRAJECTORY_RECEIPT.json').read_text())


def case(name,reader,mutate):
    with tempfile.TemporaryDirectory(prefix='mutation_',dir=HERE) as temp:
        p=Path(temp)
        (p/'inputs').mkdir()
        shutil.copyfile(HERE/'inputs/graph.json',p/'inputs/graph.json')
        shutil.copyfile(HERE/reader,p/reader)
        r,d,t,tr=map(copy.deepcopy,(base_r,base_d,base_t,base_tr))
        mutate(r,d,t,tr)
        for file,obj in [('RECEIPT.json',r),('RECORDS.json',d),('TRAJECTORIES.json',t),('TRAJECTORY_RECEIPT.json',tr)]:
            (p/file).write_text(json.dumps(obj))
        result=subprocess.run([sys.executable,'-I',str(p/reader)],capture_output=True,text=True)
        return {'name':name,'exit_code':result.returncode,'rejected':result.returncode!=0,
                'last_output':(result.stderr or result.stdout).splitlines()[-1:]}


def delete(r,d,t,tr):
    d['columns'].pop()
    r['direct_native_calls']-=1
    r['record_columns']-=1


def forged_bond(r,d,t,tr):
    for e in t['trials'][0]['events']:
        if e[4][1]is not None:
            e[4][1]=(e[4][1]+1)%12
            return


cases=[
 ('native_current_plus_one','independent_check.py',lambda r,d,t,tr:d['columns'][0].__setitem__(6,d['columns'][0][6]+1)),
 ('invalid_native_coin','independent_check.py',lambda r,d,t,tr:d['columns'][0].__setitem__(5,2)),
 ('deleted_column_adjusted_count','independent_check.py',delete),
 ('false_compact_fidelity','independent_check.py',lambda r,d,t,tr:r['baseline'].__setitem__('compact_fidelity_after_one_unread_event','1')),
 ('false_stationary_reference','independent_check.py',lambda r,d,t,tr:r['selection'].__setitem__('stationary_reference',['1/3','2/3'])),
 ('false_work_probability','independent_check.py',lambda r,d,t,tr:r['work']['probabilities'].__setitem__('0','2/3')),
 ('false_product_kernel','independent_check.py',lambda r,d,t,tr:r['separation']['independent_epoch_kernel'][0].__setitem__(0,'1')),
 ('forged_discovered_bond','trajectory_check.py',forged_bond),
 ('invalid_trajectory_coin','trajectory_check.py',lambda r,d,t,tr:t['trials'][0]['events'][0].__setitem__(1,2)),
 ('swapped_pair_selector','trajectory_check.py',lambda r,d,t,tr:t['pairs'][1]['selector_records'][1].__setitem__(0,1-t['pairs'][1]['selector_records'][1][0])),
 ('first_pair_before_state','trajectory_check.py',lambda r,d,t,tr:t['pairs'][0]['joint'][0][0].__setitem__(0,1-t['pairs'][0]['joint'][0][0][0])),
 ('first_pair_selector','trajectory_check.py',lambda r,d,t,tr:t['pairs'][1]['selector_records'][0].__setitem__(0,1-t['pairs'][1]['selector_records'][0][0])),
 ('first_pair_event_boundary','trajectory_check.py',lambda r,d,t,tr:t['pairs'][0]['joint'][0].__setitem__(2,t['pairs'][0]['joint'][0][2]+2)),
]
results=[case(*x) for x in cases]
clean=[case('unchanged_'+reader,reader,lambda r,d,t,tr:None) for reader in ('independent_check.py','trajectory_check.py')]
out={'mutations':results,'clean_controls':clean,'all_mutations_rejected':all(x['rejected'] for x in results),
     'all_clean_controls_pass':all(not x['rejected'] for x in clean)}
with (HERE/'MUTATION_RECEIPT.json').open('x') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps(out,sort_keys=True))
raise SystemExit(0 if out['all_mutations_rejected'] and out['all_clean_controls_pass'] else 1)
