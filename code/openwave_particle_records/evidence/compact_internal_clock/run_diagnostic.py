#!/usr/bin/env python3
"""Complete native arc laws for the frozen compact internal-reference study.

S calls the literal source. R and coherent preparation are declared memory
services, not a claim that the prior 12-port compiled word covers this graph.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
STATIC=HERE.parent/'compact_static_adjudication'
PROVIDER=HERE.parent/'localized_rest_join'
sys.path.insert(0,str(PROVIDER/'standalone'/'reverse-engineering-reality'))
from model.observer.native import native

V=240; P=256; M=630; K=1260; E=2048; DIM=P*E

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,indent=2)+'\n'
def pins():
    for name in ('SOURCE_PINS.json','PREREG.json','PREREG_DIAGNOSTIC.json'):
        obj=json.loads((HERE/name).read_text())
        for path,digest in obj['files'].items():
            if sha(ROOT/path)!=digest:raise RuntimeError('changed pin: '+path)

def experiment(record=False):
    started=time.monotonic();pins()
    graph=json.loads((PROVIDER/'evidence'/'graph.json').read_text())
    seams=[tuple(x) for x in graph['seams']]
    assert len(seams)==M and graph['ports']==V
    edge_id={frozenset(e):i for i,e in enumerate(seams)}
    f=np.zeros(P);f[[4,5,6,7]]=[.5,-.5,-.5,.5]
    fi=(2*f).astype(np.int64)
    dest=np.arange(DIM,dtype=np.int32).reshape(P,E)
    currents=np.zeros((V,K),dtype=np.int8)
    counts=np.zeros((V,V),dtype=np.int64)
    contrast_change=0;calls=0
    for p in range(V):
        for ei,(a,b) in enumerate(seams):
            for c in (0,1):
                out,(coin,j)=native((int(p==a),int(p==b)),(0,1),c);calls+=1
                assert coin==c and out==(int(p==a)-j,int(p==b)+j)
                q=p if p not in (a,b) else (a if out[0]==1 else b)
                ec=2*ei+c
                dest[p,ec]=q*E+2*ei+(1-c if j else c)
                currents[p,ec]=j;counts[q,p]+=1
                contrast_change+=int(fi[q]-fi[p])**2
    dest=dest.ravel();idx=np.arange(DIM)
    assert np.array_equal(dest[dest],idx)
    assert np.array_equal(np.sort(dest),idx)
    assert np.array_equal(counts@fi[:V],1254*fi[:V])
    delta=Fraction(contrast_change,2*K*int(fi@fi));lam=1-delta
    assert delta==Fraction(1,210)

    def S(x):return x[dest]
    def R(x):
        y=-x.reshape(P,E).copy()
        y[:V,:K]+=2*np.sum(x.reshape(P,E)[:V,:K],axis=1)[:,None]/K
        y[V:]=x.reshape(P,E)[V:]
        return y.ravel()
    def W(x):return S(R(x))
    u=np.zeros((P,E));u[:,:K]=f[:,None]/math.sqrt(K);u=u.ravel()
    su=S(u);s=math.sqrt(1-float(lam)**2);w=(su-float(lam)*u)/s
    theta=math.acos(float(lam))
    triangle=next(t for t in itertools.combinations(range(12),3)
                  if any(fi[x] for x in t) and
                  all(frozenset(e) in edge_id for e in itertools.combinations(t,2)))
    vi=np.zeros(DIM,dtype=np.int8);cycle_arcs=[]
    for a,b in zip(triangle,triangle[1:]+triangle[:1]):
        ei=edge_id[frozenset((a,b))];left,right=seams[ei]
        c=0 if a==left else 1
        z=a*E+2*ei+c;zr=int(dest[z]);assert zr//E==b
        vi[z]=1;vi[zr]=-1
        cycle_arcs.append({'from':a,'to':b,'event':2*ei+c,'current':int(currents[a,2*ei+c]),'reverse_arc':zr})
    assert np.array_equal(vi[dest],-vi)
    assert np.count_nonzero(vi)==6
    assert np.array_equal(vi.reshape(P,E).sum(axis=1),np.zeros(P,dtype=int))
    v=vi/math.sqrt(6);chi=(u+v)/math.sqrt(2)
    B=np.array([u,w,v]);orth=float(np.max(abs(B@B.T-np.eye(3))))
    assert orth<2e-11
    state_error=max(float(np.max(abs(W(u)-(float(lam)*u+s*w)))),
                    float(np.max(abs(W(w)-(-s*u+float(lam)*w)))),
                    float(np.max(abs(W(v)-v))))
    assert state_error<2e-11
    support=set(np.flatnonzero(fi))
    envelope=support|set(triangle)
    for a,b in seams:
        if a in support:envelope.add(b)
        if b in support:envelope.add(a)
    outside=np.array([p not in envelope for p in range(P)])
    kept=np.flatnonzero((u!=0)|(su!=0)|(vi!=0))
    omitted_count=DIM-len(kept)
    witness=next(int(z) for z in np.flatnonzero(vi) if fi[z//E])
    a,b,d=u[kept],w[kept],v[kept]
    coef=np.array([d*d/2+(a*a+b*b)/4,d*a,d*b,(a*a-b*b)/4,a*b/2])
    ko=int(np.searchsorted(kept,witness))
    assert abs(coef[1,ko])>1e-8
    states={'compact':u.copy(),'internal_reference':chi.copy()}
    probabilities={key:[] for key in states};energy={key:[] for key in states}
    energy['birth']=[];mass_sums=[];max_leak=0.;formula_error=0.
    hold=chi.copy();phase_plus=chi.copy();changed_dist={'held_source':0.,'phase_plus':0.,'double_frequency':0.,'dephased':0.}
    held_laws=[];plus_laws=[];dephased_laws=[];double_laws=[]
    def energy_click(x):
        amplitudes=f@x.reshape(P,E)
        return float(amplitudes@amplitudes)
    for n in range(129):
        tr=np.array([1,math.cos(n*theta),math.sin(n*theta),math.cos(2*n*theta),math.sin(2*n*theta)])
        predicted=tr@coef
        ideal=(math.cos(n*theta)*u+math.sin(n*theta)*w+v)/math.sqrt(2)
        state_error=max(state_error,float(np.max(abs(states['internal_reference']-ideal))))
        pchi=states['internal_reference'][kept]**2
        formula_error=max(formula_error,float(np.max(abs(pchi-predicted))))
        for key,x in states.items():
            probabilities[key].append(x[kept]**2)
            energy[key].append([energy_click(x),1-energy_click(x)])
            mass_sums.append(float(x@x))
            max_leak=max(max_leak,float(np.sum(x.reshape(P,E)[outside]**2)))
        energy['birth'].append([energy_click(v),1-energy_click(v)])
        dep=(states['compact'][kept]**2+v[kept]**2)/2
        double=((math.cos(2*n*theta)*u+math.sin(2*n*theta)*w+v)/math.sqrt(2))[kept]**2
        for key,law in [('held_source',hold[kept]**2),('phase_plus',phase_plus[kept]**2),('double_frequency',double),('dephased',dep)]:
            changed_dist[key]=max(changed_dist[key],float(np.sum(abs(law-pchi))/2))
        held_laws.append(hold[kept]**2);plus_laws.append(phase_plus[kept]**2)
        dephased_laws.append(dep);double_laws.append(double)
        states={key:W(x) for key,x in states.items()}
        hold=R(hold);phase_plus=S(phase_plus)
    assert state_error<2e-11 and formula_error<2e-11
    assert max_leak==0 and max(abs(x-1) for x in mass_sums)<2e-11, {"max_outside_probability":max_leak,"norm_residual":max(abs(x-1) for x in mass_sums),"state_residual":state_error,"record_residual":formula_error}
    assert all(d>1e-7 for d in changed_dist.values())
    # No hidden support can enter the controls: S and R preserve this union.
    assert np.count_nonzero(hold[np.setdiff1d(idx,kept)])==0
    onset=chi[kept]**2
    onset_average=float(onset.sum()/(V*K))
    vacuum_records=[]
    for ei in range(M):
        for c in (0,1):
            out,rec=native((0,0),(0,1),c);calls+=1
            assert out==(0,0) and rec==(c,0)
            vacuum_records.append([ei,c,0])
    elog=-math.log(float(lam));E6=K*elog
    ratios={'theta_over_E6_word_count':theta/E6,'two_theta_over_E6_word_count':2*theta/E6,
            'theta_over_log_energy_attempt_count':theta/elog,'two_theta_over_log_energy_attempt_count':2*theta/elog}
    scientific={
        'candidate':'A-signed-compact-arc-clock-v1',
        'source':{'native_calls':calls,'charged_native_columns':V*K,'vacuum_native_columns':K,
                  'full_arc_dimension':DIM,'body_arcs':V*K,'event_padding_arcs':V*(E-K),
                  'invalid_position_arcs':(P-V)*E,'source_current_sum_squares':contrast_change,
                  'D_over_Var':str(delta),'lambda':str(lam),'involution_columns':DIM},
        'cycle':{'triangle':list(triangle),'arcs':cycle_arcs,'integer_divergence_zero':True,
                 'native_reversal_sign':-1,'walk_eigenvalue':1},
        'clock':{'theta':theta,'return_probability_frequency':2*theta,'internal_reference_frequency':theta,
                 'witness_arc':witness,'witness_position':witness//E,'witness_event':witness%E,
                 'witness_coefficients':coef[:,ko].tolist(),'all_position_event_labels_retained':DIM,
                 'explicit_nonzero_support_labels':len(kept),'default_zero_labels':omitted_count,
                 'principal_walk_energy_branches':[0,-theta,theta]},
        'localization':{'all_integer_times_proved':True,'position_envelope':sorted(envelope),
                        'within_original_twelve_port_cell':all(p<12 for p in envelope),'max_outside_probability':max_leak},
        'energy':{'sealed_E6':E6,'sealed_per_attempt_log_energy':elog,'mass_clock_ratios':ratios,
                  'mass_clock_gates':{k:abs(val-1)<=1e-8 for k,val in ratios.items()},
                  'compact_ray_click_min':min(x[0] for x in energy['compact']),
                  'compact_ray_click_max':max(x[0] for x in energy['compact']),
                  'internal_reference_initial_ray_click':energy['internal_reference'][0][0],
                  'birth_ray_click':energy['birth'][0][0],
                  'unresolved_is_not_zero_energy':True},
        'onset':{'rank_one_filter_success_on_uniform_valid_arc':onset_average,
                 'rank_one_filter_failure_on_uniform_valid_arc':1-onset_average,
                 'prepared_chi_success':float(chi@chi)**2,'native_vacuum_charged_onset_probability':0,
                 'autonomous_onset_established':False,'filter_is_active_preparation':True},
        'controls_max_complete_total_variation':changed_dist,
        'diagnostics':{'orthonormal_residual':orth,'walk_state_residual':state_error,
                       'full_record_harmonic_residual':formula_error,
                       'norm_residual':max(abs(x-1) for x in mass_sums)},
        'scope':{'row3_full':False,'row4_full':False,'declared_coherent_memory_law':True,
                 'new_240_port_reflection_compiler_executed':False,'new_19_bit_pointer_executed':False,
                 'self_start_physical_energy_and_free_law_attachment_open':True},
    }
    if record:
        evidence=HERE/'evidence';evidence.mkdir(exist_ok=True)
        np.savez_compressed(evidence/'native_arcs.npz',destination=dest,current=currents,transfer_counts=counts,
                            vacuum_records=np.asarray(vacuum_records,dtype=np.int32))
        np.savez_compressed(evidence/'complete_records.npz',labels=kept,default_zero_count=omitted_count,
                            compact=np.asarray(probabilities['compact']),internal_reference=np.asarray(probabilities['internal_reference']),
                            dephased=np.asarray(dephased_laws),held_source=np.asarray(held_laws),phase_plus=np.asarray(plus_laws),
                            double_frequency=np.asarray(double_laws),energy_compact=np.asarray(energy['compact']),
                            energy_internal_reference=np.asarray(energy['internal_reference']),energy_birth=np.asarray(energy['birth']),
                            onset_success=onset,onset_failure=1-onset,theta_harmonic_coefficients=coef)
        artifacts={str(p.relative_to(HERE)):sha(p) for p in sorted(evidence.iterdir())}
        receipt={'scientific':scientific,'artifacts':artifacts,'elapsed_seconds':time.monotonic()-started,
                 'prereg_sha256':sha(HERE/'PREREG.json')}
        (HERE/'RECEIPT.json').write_text(canonical(receipt))
    return scientific,time.monotonic()-started

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args()
    result,elapsed=experiment(args.record)
    if not args.record:
        old=json.loads((HERE/'RECEIPT.json').read_text())
        assert old['scientific']==result,'scientific replay changed'
        for path,digest in old['artifacts'].items():assert sha(HERE/path)==digest
    print(canonical({'scientific':result,'elapsed_seconds':elapsed}))
    if elapsed>=60:raise SystemExit('headline exceeded60 seconds')
