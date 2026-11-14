"""Native finite budgets, sharp-record and flat-sector controls."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import math
import numpy as np
import native_spectroscopy as ns
from run import bloch_source,peak,pins

HERE=Path(__file__).resolve().parent


def main():
    before=pins();budgets=[]
    for size in (2,4):
        for mass in (F(0),F(1,2),F(1)):
            budgets.append(ns.certificate(size,4,mass,precision=96))
    source=ns.walk.Walk(2,mass=F(0))
    g=source.geometry
    columns=np.zeros((g.graph.ports,*g.shape),complex)
    for v,address in enumerate(g.addresses):columns[v,address,:12]=1/math.sqrt(12)
    laws=[]
    for mass in (F(0),F(1,2),F(1)):
        w=ns.walk.Walk(2,h=F(1,2),mass=mass)
        out=w.step(columns)
        law=(abs(out[:,g.addresses])**2).sum(axis=2).T
        laws.append(law)
    independence=max(peak(x-laws[0]) for x in laws)
    assert independence<1e-14
    T=ns.a.common.native_transfer()
    assert peak(T-laws[0])<1e-14 and peak(T-T.T)<1e-14
    # Construct an actual internal triangle circulation from the carrier's
    # seam/event slots. Every vertex has one incoming and one outgoing edge.
    lookup={tuple(sorted(edge)):direction for direction,matching in enumerate(ns.a.common.INTERNAL_MATCHINGS) for edge in matching}
    triangle=next(t for t in combinations(range(12),3) if all(e in lookup for e in combinations(t,2)))
    cycle=np.zeros(144,complex)
    for x,y in zip(triangle,triangle[1:]+triangle[:1]):
        event=2*lookup[tuple(sorted((x,y)))]+1
        cycle[12*x+event]+=1;cycle[12*y+event]-=1
    cycle/=np.linalg.norm(cycle)
    flat=[]
    for k in ((.2,0,0),(.1,.2,.3)):
        W,S,R=bloch_source(k,.25)
        err=float(np.linalg.norm(W@cycle-np.exp(.25j)*cycle))
        assert err<1e-14
        flat.append({'k':list(k),'mass_coin_angle':.25,'phase':.25,'residual':err})
    # An omitted I generator is observable against the inactive control.
    theta=.37;rho=np.ones((2,2))/2
    from run import readout,contrast
    c=complex(np.exp(1j*theta))
    controlled=np.array([[1,c.conjugate()],[c,1]])/2
    difference=abs(contrast(readout(controlled,'Y','C'))-contrast(readout(rho,'Y','C')))
    assert difference>.3
    # Complete local first-export records are mandatory even at tiny weight.
    bad=ns.b.complete_run(ns.b.H(0),1,np.array([1,0]),96)
    bad['failure_complement_counts']=[]
    try:ns.b.require_complete(bad)
    except ValueError:missing_rejected=True
    else:raise AssertionError('missing failure accounting accepted')
    assert pins()==before
    result={'finite_native_budgets':budgets,
      'sharp_fresh_event_every_step':{'all_initial_ports':g.graph.ports,'masses':['0','1/2','1'],
        'mass_independence_residual':independence,'actual_native_transfer_residual':peak(T-laws[0]),
        'real_symmetric_transfer':True,'eigenvalue_interval':[float(np.linalg.eigvalsh(T)[0]),float(np.linalg.eigvalsh(T)[-1])],
        'coherent_KG_transfer_claimed':False},
      'retained_flat_cycle':{'triangle':list(triangle),'cases':flat,'not_a_propagating_KG_branch':True},
      'controlled_scalar_phase_omission_observable_difference':difference,
      'missing_export_accounting_rejected':missing_rejected,'all_provider_pins_unchanged':True}
    (HERE/'ADDITIONAL_CONTROLS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'finite_budgets':len(budgets),'sharp_record_ports':g.graph.ports,'flat_cycle_cases':len(flat),'scalar_phase_control':difference}),flush=True)


if __name__=='__main__':main()
