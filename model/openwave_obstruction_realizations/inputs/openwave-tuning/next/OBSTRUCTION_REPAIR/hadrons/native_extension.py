"""Execute every input column of the added current-square native block.

This does not transfer the old preparation or complete old row-13 program.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import numpy as np
import audit

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
BRIDGE=ROOT/'openwave-tuning/common_dynamics/color_chain_native'


def main():
    pins={str((BRIDGE/f).relative_to(ROOT)):hashlib.sha256((BRIDGE/f).read_bytes()).hexdigest()
          for f in ['native_path.py','run_native.py','source.py','MANIFEST.json']}
    sys.path.insert(0,str(BRIDGE))
    sys.modules['source']=audit.source
    spec=importlib.util.spec_from_file_location('colour_native_execution_helper',BRIDGE/'run_native.py')
    runner=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    results=[]
    for weakup in [1,2]:
        d=audit.source.data(2,F(1),weakup)
        values=audit.extension_values(d,'current_square')
        time=F(1,8)
        word=runner.p.diagonal_word(values,(F(),F(),time),tuple(range(4)),((4,1),))
        target=np.eye(32,dtype=complex)
        target[16:,16:]=np.diag(np.exp(-1j*float(time)*np.array(list(map(float,values)))))
        result,_,_=runner.run('new current-square controlled block weakup='+str(weakup),word,5,target)
        result['source_diagonal_exact']=list(map(str,values))
        result['source_time_exact']=str(time)
        result['false_control_and_padding_identity']=True
        results.append(result)
    for path,digest in pins.items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    out=dict(schema='oph.obstruction_repair.current_square_native_block.v1',live_helper_hashes=pins,
        blocks=results,total_executed_pulses=sum(r['pulses'] for r in results),
        all_columns=True,all_first_export_complements_retained=True,
        full_amended_preparation_or_decay_executed=False,
        scope='Native implementation of the new rational local-current diagonal only; the old row13 filter is not transferred.')
    (HERE/'NATIVE_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Executed',out['total_executed_pulses'],'native pulses; full amended experiment remains open.')


if __name__=='__main__':
    main()
