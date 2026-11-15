"""Actual zero-defect source history and native zero-pointer intervention."""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import sys
import topology
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PROVIDER=ROOT/'reverse-engineering-reality/code/openwave_row1_integrity'
AUDIT=ROOT/'openwave-tuning/next/CHAIN3/row1_integrity_audit'
sys.path[:0]=[str(PROVIDER),str(AUDIT)]
import integrity as it
from observer_only import admit_record_only

def main():
    paths=[Path(__file__),HERE/'topology.py',PROVIDER/'integrity.py',PROVIDER/'native_kernel.py',
           PROVIDER/'native_cycle_meter.py',PROVIDER/'strict_tape_validator.py',PROVIDER/'closing_sweep.py',AUDIT/'observer_only.py']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    G=(1,)+(0,)*11;region=(0,);order=list(range(7,30))+list(range(7))
    # This is one finite prescribed branch with explicit probability, not
    # an assertion that random source evolution stays at its preparation.
    def hold(e):return next(c for c in (0,1) if it.native(G,it.EDGES[e],c)[1][1]==0)
    prefix=[(0,hold(0),'C'),(7,hold(7),'E')]
    received,audit=it.source_execute(G,prefix,[hold(e) for e in order],order)
    assert audit['source_L1']==list(G)
    assert all(r['current']==0 for r in received['source_records']+received['sweep_records'])
    geo=topology.construct(it.EDGES)
    assert all(topology.discrete_curve_degree(geo,[0]*30,(v,))==0 for v in range(12))
    pointer=it.meter.execute(11,(0,),it.meter.favorable_coins(11,(0,),2),2)
    claim=it.assemble_claim(G,received,region,pointer,2);del claim['meter']['final_loads']
    verdict=admit_record_only(claim,G,region,2,producer=it)
    assert verdict['admitted'] and verdict['readout']==0
    assert topology.raw_pointer_degree(pointer,it.meter.CYCLE)==0 and len(pointer['records'])>0
    # A zero answer is not allowed to bypass physical discovery and reset.
    fake=deepcopy(claim);fake['meter']['records']=[]
    assert not admit_record_only(fake,G,region,2,producer=it)['admitted']
    out={'schema':'oph.row1.zero_defect.control.v1','G':G,'received':received,'audit':audit,
         'claim':claim,'admission':verdict,'all_twelve_spatial_degrees':list([0]*12),
         'zero_readout_requires_actual_native_reset':True,'source_fine_history_probability':'1/4294967296',
         'conditional_selected_branch':True,'source_sha256':pins}
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in pins.items())
    (HERE/'VACUUM.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'native_pointer_events':len(pointer['records']),'readout':verdict['readout'],'empty_zero_echo_rejected':True}))

if __name__=='__main__':main()
