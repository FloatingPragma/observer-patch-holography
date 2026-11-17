#!/usr/bin/env python3
"""Actual native loops and complete common current clock; replay is read-only."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy import sparse as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
provider=HERE.parent/'native_current_clock/probe.py'
spec=importlib.util.spec_from_file_location('field_cycle_current_provider',provider)
cc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=cc;spec.loader.exec_module(cc)
sg=cc.sg

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return json.loads(json.dumps(x,sort_keys=True,allow_nan=False))
def maxabs(x):return cc.maxabs(x)

def loop_columns(states,q,arrays):
    n=len(states);idx={s:i for i,s in enumerate(states)}
    targets=np.zeros((2,n),int);signs=np.ones((2,n),int)
    increments=np.zeros((2,n,3),int);returns=np.zeros((2,n),bool)
    rows=[];words=(sg.EDGES,tuple(reversed(sg.EDGES)))
    for direction,word in enumerate(words):
        for col,initial in enumerate(states):
            s=initial;phase=1
            for step,seam in enumerate(word):
                before=s;hop,sign,_=sg.joined(s,(0,1),seam,limit=1)
                phase*=sign;s=hop.state
                a,b=seam;coin=0 if a in before.pair and b not in before.pair else 1
                loads=tuple(int(v in before.pair) for v in range(12))
                _,native_record=cc.native.native(loads,seam,coin)
                assert native_record[1]==hop.proposed_number_current
                rows.append((direction,col,step,sg.EDGES.index(seam),idx[before],idx[s],
                             hop.proposed_number_current,hop.realized_signed_current,
                             int(hop.overflow),sign,coin,native_record[1]))
                assert sg.g.gauss(before)==sg.g.gauss(s)
            existing,history=sg.g.native_word(initial,seams=word,limit=1)
            assert existing==s and len(history)==3
            targets[direction,col]=idx[s];signs[direction,col]=phase
            increments[direction,col]=[s.flux[sg.g.EDGE_INDEX[e]]-initial.flux[sg.g.EDGE_INDEX[e]] for e in sg.EDGES]
            returns[direction,col]=initial.tagged_port==s.tagged_port
            if returns[direction,col]:
                full=[b-a for a,b in zip(initial.flux,s.flux)]
                assert not any(sg.g.divergence(full))
    matrices=[sp.coo_matrix((signs[d],(targets[d],np.arange(n))),shape=(n,n)).tocsr() for d in range(2)]
    U,V=matrices;I=sp.eye(n,format='csr')
    assert maxabs(U.getH()-V)==0 and maxabs(U.getH()@U-I)==0
    P=(6*I+U+V)/8
    eig=np.linalg.eigvalsh(P.toarray())
    tag0=np.array([i for i,s in enumerate(states) if s.tagged_port==0])
    K=U[tag0][:,tag0];P0=(6*sp.eye(len(tag0))+K+K.getH())/8
    eig0=np.linalg.eigvalsh(P0.toarray())
    assert eig.min()>=.5-1e-9 and eig.max()<=1+1e-9
    assert eig0.min()>=.5-1e-9 and eig0.max()<=1+1e-9
    # Complete lazy source instrument: all eight branch and terminal coordinates.
    branch_targets=np.tile(np.arange(n),(8,1));branch_signs=np.ones((8,n),int)
    branch_targets[6:]=targets;branch_signs[6:]=signs
    branch_return=np.zeros((8,len(tag0)),bool)
    for b in range(8):branch_return[b]=[states[t].tagged_port==0 for t in branch_targets[b,tag0]]
    assert np.all(np.sum(np.full((8,n),1/8),axis=0)==1)
    # Interior cycle prediction and whole-boundary translation control.
    interior=[];illegal_shift=0
    for col in tag0:
        s=states[col];local=np.array([s.flux[sg.g.EDGE_INDEX[e]] for e in sg.EDGES])
        shifted=local+np.array([-q,-q,q])
        if np.all(abs(shifted)<=1):
            interior.append(int(col));assert returns[0,col]
            assert np.array_equal(increments[0,col],[-q,-q,q])
        else:illegal_shift+=1
    arrays.update({f'q{q}_loop_targets':targets,f'q{q}_loop_signs':signs,
                   f'q{q}_loop_increments':increments,f'q{q}_loop_returns':returns,
                   f'q{q}_loop_native_records':np.array(rows,dtype=np.int64),
                   f'q{q}_lazy_branch_targets':branch_targets,f'q{q}_lazy_branch_signs':branch_signs,
                   f'q{q}_lazy_tag0_return_flags':branch_return,
                   f'q{q}_lazy_kernel_eigenvalues':eig,f'q{q}_lazy_return_eigenvalues':eig0})
    return U,targets,signs,increments,returns,{
        'native_path_calls':len(rows),'full_loop_columns_per_direction':n,
        'returned_columns_per_direction':returns.sum(axis=1).tolist(),
        'nonreturned_columns_per_direction':(~returns).sum(axis=1).tolist(),
        'tag0_source_columns':len(tag0),'tag0_return_counts_by_lazy_address':branch_return.sum(axis=1).tolist(),
        'lazy_kernel_spectral_interval':[float(eig.min()),float(eig.max())],
        'lazy_return_spectral_interval':[float(eig0.min()),float(eig0.max())],
        'interior_tag0_forward_columns':len(interior),'boundary_translation_mutant_illegal_columns':illegal_shift,
        'returned_loop_divergence_exact_zero':True,'random_unitary_channel_equals_amplitude_kernel':False}

def common_clock(target,current,sign,spins,U,mode='native'):
    n=target.shape[1];N=n+1
    valid=sp.diags([1]*n+[0],format='csr');invalid=sp.diags([0]*n+[1],format='csr')
    active_spins=spins if mode!='spin_blind' else np.tile(np.eye(2),(3,1,1))
    first=cc.select_matrix(target,current,sign,active_spins,'zero_current' if mode=='zero_current' else 'native')[:12*N,:12*N]
    Uext=sp.block_diag((U,sp.csr_matrix((1,1))),format='csr')
    up=sp.kron(sp.eye(2),Uext,format='csr');bad=sp.kron(sp.eye(2),invalid,format='csr')
    loop=sp.bmat([[bad,up.getH()],[up,bad]],format='csr')
    if mode=='omit_loop':loop=sp.eye(4*N,format='csr')
    S=sp.block_diag((first,loop,*[sp.eye(4*N)]*4),format='csr')
    seed=np.zeros(16);seed[[0,2,4,8,10,12,14]]=1/np.sqrt(8);seed[[6,7]]=1/4
    reflect=2*np.outer(seed,seed)-np.eye(16)
    R=sp.kron(sp.csr_matrix(reflect),sp.kron(sp.eye(2),valid))+sp.kron(sp.eye(32),invalid)
    inclusion=sp.coo_matrix((np.ones(n),(np.arange(n),np.arange(n))),shape=(N,n)).tocsr()
    T=sp.kron(sp.csr_matrix(seed[:,None]),sp.kron(sp.eye(2),inclusion),format='csr')
    J=cc.current_operators(target,current if mode!='zero_current' else np.zeros_like(current),sign)
    sumcurrent=sum(sp.kron(sp.csr_matrix(active_spins[e]),J[e]) for e in range(3))
    loopreal=sp.eye(n) if mode=='omit_loop' else (U+U.getH())/2
    A=(4*sp.eye(2*n)+sumcurrent+sp.kron(sp.eye(2),loopreal))/8
    return S,R,S@R,T,A,J

def record_diagonals(states,q,current,targets,increments,returns,mode):
    n=len(states);N=n+1
    currents=[];returned=[];nonreturned=[];electric=[]
    for f,seam in enumerate(sg.EDGES):
        vals=[s.flux[sg.g.EDGE_INDEX[seam]] for s in states]+[0]
        electric.append(np.tile(vals,32))
        cq=np.zeros((8,2,2,N));lr=cq.copy();ln=cq.copy()
        if mode!='zero_current':cq[f,:,:,:n]=q*current[f]
        if mode!='omit_loop':
            for d in range(2):
                lr[3,d,:,:n]=increments[d,:,f]*returns[d]
                ln[3,d,:,:n]=increments[d,:,f]*(~returns[d])
        currents.append(cq.ravel());returned.append(lr.ravel());nonreturned.append(ln.ravel())
    return np.array(electric),np.array(currents),np.array(returned),np.array(nonreturned)

def spin_resolved_states(A,S,T,B,normal,dim):
    small=B.conj().T@(A@B)
    assert np.linalg.norm(A@B-B@small)<1e-9
    eig,V=np.linalg.eigh(small)
    Snormal=np.einsum('i,iab->ab',normal,sg.SIGMA)/2
    N=dim//32;spinop=sp.kron(sp.eye(16),sp.kron(sp.csr_matrix(Snormal),sp.eye(N)))
    groups=[];start=0
    while start<len(eig):
        end=start+1
        while end<len(eig) and abs(eig[end]-eig[start])<1e-10:end+=1
        groups.append((start,end));start=end
    results=[]
    for lo,hi in groups:
        lam=float(np.mean(eig[lo:hi]));signal=B@V[:,lo:hi];u=T@signal
        beta=np.sqrt(max(0.,1-lam*lam))
        phase_signs=(1,-1) if beta>1e-10 else (1,)
        for phase_sign in phase_signs:
            psi=(u-1j*phase_sign*(S@u-lam*u)/beta)/np.sqrt(2) if beta>1e-10 else u
            polarized,rotation=np.linalg.eigh(psi.conj().T@(spinop@psi))
            psi=psi@rotation;signal_rotated=signal@rotation
            for k in range(hi-lo):
                results.append((lam,complex(lam,phase_sign*beta),psi[:,k],signal_rotated[:,k],float(polarized[k]),[lo,hi],phase_sign))
    return eig,results,spinop

def stationary_diagnostics(states,q,mode,S,R,W,T,A,J,diags,normal,weights,arrays):
    E,C,LR,LN=diags;B=cc.physical_signal_basis(states)
    eig,prepared,spinop=spin_resolved_states(A,S,T,B,normal,W.shape[0])
    entries=[];vectors=[];maxres=0.;maxbalance=0.
    for index,(lam,phase,psi,signal,spin,band,phase_sign) in enumerate(prepared):
        after=R@psi;prob=abs(after)**2
        charge=C@prob;ret=LR@prob;nonret=LN@prob
        drift=E@(abs(W@psi)**2-abs(psi)**2)
        balance=-charge+ret+nonret
        residual=float(np.linalg.norm(W@psi-phase*psi))
        maxres=max(maxres,residual);maxbalance=max(maxbalance,float(max(abs(balance))))
        assert residual<1e-9 and max(abs(balance))<1e-9 and max(abs(drift))<1e-9
        coherent=[float(np.vdot(signal,sp.kron(sp.eye(2),j)@signal).real) for j in J]
        moment=weights@charge/2
        entries.append({'index':index,'band':band,'lambda':lam,'phase_sign':phase_sign,
                        'phase_radians':float(np.angle(phase)),'spin_normal':spin,
                        'realized_charge_current':charge.tolist(),'returned_loop_electric_increment':ret.tolist(),
                        'nonreturned_loop_electric_increment':nonret.tolist(),'balance':balance.tolist(),
                        'stationary_electric_increment':drift.tolist(),'coherent_J_quadratures':coherent,
                        'axial_transport_moment':moment.tolist(),'normal_transport_moment':float(normal@moment),
                        'eigenstate_residual':residual})
        vectors.append(psi)
    arrays[f'q{q}_{mode}_stationary_vectors']=np.array(vectors)
    arrays[f'q{q}_{mode}_signal_eigenvalues']=eig
    maxcurrent=max(max(abs(np.array(x['realized_charge_current']))) for x in entries)
    maxret=max(max(abs(np.array(x['returned_loop_electric_increment']))) for x in entries)
    maxmoment=max(abs(x['normal_transport_moment']) for x in entries)
    return {'signal_dimension':len(eig),'stationary_state_count':len(entries),
            'amplitude_spectral_interval':[float(eig.min()),float(eig.max())],
            'max_stationary_residual':maxres,'max_electric_balance_residual':maxbalance,
            'max_realized_charge_current':maxcurrent,'max_returned_loop_electric_response':maxret,
            'max_normal_transport_moment':maxmoment,'states':entries}

def run_case(q,frame,spins,normal,weights,arrays):
    states,target,current,sign,overflow,proposed,raw=cc.source_columns(q)
    n=len(states);N=n+1
    U,lt,ls,inc,ret,loop_info=loop_columns(states,q,arrays)
    arrays.update({f'q{q}_source_targets':target,f'q{q}_source_currents':current,
                   f'q{q}_source_signs':sign,f'q{q}_source_overflow':overflow,
                   f'q{q}_source_native_records':raw})
    modes={};reference={};allerrors={}
    for mode in ('native','omit_loop','zero_current','spin_blind'):
        S,R,W,T,A,J=common_clock(target,current,sign,spins,U,mode)
        dim=W.shape[0];I=sp.eye(dim,format='csr')
        errors={'select_hermiticity':maxabs(S-S.getH()),'select_involution':maxabs(S@S-I),
                'reflection_involution':maxabs(R@R-I),'walk_unitarity':maxabs(W.getH()@W-I),
                'compression':maxabs(T.getH()@S@T-A)}
        diags=record_diagonals(states,q,current,lt,inc,ret,mode)
        E,C,LR,LN=diags
        for f in range(3):
            ee=sp.diags(E[f],format='csr');delta=sp.diags(-C[f]+LR[f]+LN[f],format='csr')
            errors[f'electric_balance_{f}']=maxabs(W.getH()@ee@W-ee-R.getH()@delta@R)
        for vertex in range(12):
            gv=[sg.g.gauss(s)[vertex] for s in states]+[0]
            G=sp.diags(np.tile(gv,32),format='csr')
            errors[f'gauss_{vertex}']=maxabs(W@G-G@W)
        assert max(errors.values())<1e-9
        invalid=np.zeros(dim,complex);invalid[n]=1;assert np.array_equal(W@invalid,invalid)
        signal=np.zeros(2*n,complex)
        for memory in range(6):signal[states.index(sg.g.pair_string((0,1),q,memory))]=1/np.sqrt(6)
        initial=T@signal;v=initial.copy();laws=[];recur=[];t0=signal.copy();t1=A@signal
        for tick in range(9):
            law=np.array([cc.reference_read(initial,v,z) for z in (1,1j)])
            assert max(abs(law.sum(axis=(1,2))-1))<1e-9
            laws.append(law)
            expected=t0 if tick==0 else t1
            recur.append(float(np.linalg.norm(T.getH()@v-expected)))
            if tick>=1:t0,t1=t1,2*(A@t1)-t0
            v=W@v
        assert max(recur)<1e-9
        reference[mode]=np.array(laws);arrays[f'q{q}_{mode}_reference_laws']=reference[mode]
        modes[mode]=stationary_diagnostics(states,q,mode,S,R,W,T,A,J,diags,normal,weights,arrays)
        modes[mode]['chebyshev_error']=max(recur)
        allerrors[mode]=errors
    changes={k:float(np.max(abs(v-reference['native']))) for k,v in reference.items() if k!='native'}
    assert changes['omit_loop']>1e-6 and changes['zero_current']>1e-6
    assert modes['omit_loop']['max_realized_charge_current']<1e-9
    return {'charge_number':q,'source_columns':n,'full_walk_columns':32*N,
            'source_native_calls':len(raw),'loop':loop_info,'matrix_errors':allerrors,
            'changed_law_reference_differences':changes,'stationary_diagnostics':modes,
            'invalid_failure_probability':1,'spectator_spin':'two-state full tensor identity',
            'no_far_field_dipole_g_or_larmor_inferred':True}

def plant_reader(weights,normal):
    cycle=np.array([1,1,-1]);out=[]
    for current in (cycle,-cycle,np.zeros(3)):
        moment=weights@current/2
        # Separate polygon area from the record-current sum.
        frame,_,_=cc.geometric_gates();a,b,c=(frame[x] for x in sg.TRIANGLE)
        expected=np.cross(b-a,c-a)/2*current[0]
        assert np.linalg.norm(moment-expected)<1e-12
        out.append({'current':current.tolist(),'moment':moment.tolist(),'normal':float(normal@moment)})
    assert abs(out[0]['normal'])>1e-3 and out[1]['normal']==-out[0]['normal']
    return {'cases':out,'zero_reader_mutant_error':abs(out[0]['normal']),
            'scope':'circulation reader only; this is not a planted spatial dipole detector'}

def pins():
    inherited=json.loads((HERE.parent/'native_current_clock/SOURCE_PINS.json').read_text())
    paths=set(inherited)
    for name in ('probe.py','PROOF.md','LOOP_CONTINUATION.md','MANIFEST.sha256','RECEIPT.json','INTERVAL_PROOF.md'):
        paths.add(str((HERE.parent/'native_current_clock'/name).relative_to(ROOT)))
    for name in ('openwave/row_05_magnetic_moment/THEOREMS.md','openwave/row_06_spin_half/THEOREMS.md',
                 'reverse-engineering-reality/model/observer/native_loop_readout.py',
                 'reverse-engineering-reality/Lean/Screen/SeamU1HolonomyClassification.lean'):
        paths.add(name)
    return {p:digest(ROOT/p) for p in sorted(paths)}

def freeze():
    assert not (HERE/'RECEIPT.json').exists()
    source=pins();(HERE/'SOURCE_PINS.json').write_text(json.dumps(source,sort_keys=True,indent=2)+'\n')
    files={p:digest(HERE/p) for p in ('DESIGN.md','TASK.md','probe.py','SOURCE_PINS.json')}
    (HERE/'PREREG.json').write_text(json.dumps({'schema':'oph.native_field_cycle.prereg.v1',
          'files':files,'evidence_executed':False},sort_keys=True,indent=2)+'\n')
    print(json.dumps({'preregistration_sha256':digest(HERE/'PREREG.json'),'source_pin_count':len(source)}))

def check_pins():
    for p,h in json.loads((HERE/'PREREG.json').read_text())['files'].items():assert digest(HERE/p)==h,p
    for p,h in json.loads((HERE/'SOURCE_PINS.json').read_text()).items():assert digest(ROOT/p)==h,p

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true');parser.add_argument('--record',action='store_true')
    args=parser.parse_args()
    if args.freeze:freeze();return
    check_pins()
    if args.record and (HERE/'RECEIPT.json').exists():raise FileExistsError('Evidence is immutable')
    start=time.perf_counter();frame,_,spins=cc.geometric_gates()
    a,b,c=(frame[x] for x in sg.TRIANGLE);normal=np.cross(b-a,c-a);normal/=np.linalg.norm(normal)
    weights=np.array([np.cross(frame[a],frame[b]) for a,b in sg.EDGES]).T
    arrays={'frame':frame,'normal':normal,'axial_weights':weights};plants=plant_reader(weights,normal)
    cases=[]
    for q in (0,1):
        cases.append(run_case(q,frame,spins,normal,weights,arrays));print('Completed q='+str(q),flush=True)
    receipt=canonical({'schema':'oph.native_field_cycle.v1','candidate':'CAR-current-loop-clock-v1',
                       'class':'N+ with declared coherent services','cases':cases,'plants':plants,
                       'nature_measurements_used':[],'physical_row5_claimed':False,
                       'full_integer_and_added_reflection_native_pulse_compilation':False,
                       'preregistration_sha256':digest(HERE/'PREREG.json')})
    check_pins()
    if args.record:
        (HERE/'RECEIPT.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
        np.savez_compressed(HERE/'RECORDS.npz',**arrays)
    else:
        assert receipt==json.loads((HERE/'RECEIPT.json').read_text())
        with np.load(HERE/'RECORDS.npz') as old:
            assert set(old.files)==set(arrays)
            for name,value in arrays.items():assert np.array_equal(value,old[name]),name
    print(json.dumps({'all_pass':True,'elapsed_seconds':time.perf_counter()-start,
        'native_calls':sum(c['source_native_calls']+c['loop']['native_path_calls'] for c in cases),
        'stationary_current':[c['stationary_diagnostics']['native']['max_realized_charge_current'] for c in cases],
        'returned_field_response':[c['stationary_diagnostics']['native']['max_returned_loop_electric_response'] for c in cases],
        'physical_rows_complete':[]},indent=2))

if __name__=='__main__':main()
