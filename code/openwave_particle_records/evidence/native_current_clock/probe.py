#!/usr/bin/env python3
"""Complete source-current query walk. Replay is read-only."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy import sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RER = ROOT/'reverse-engineering-reality'
sys.path.insert(0, str(RER))
provider_path = HERE.parent/'spin_gauge_join/verify.py'
spec = importlib.util.spec_from_file_location('current_clock_spin_provider', provider_path)
sg = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = sg
spec.loader.exec_module(sg)
from model.observer import native
from model.observer.source_scalar import compiler
from model.native_observer_program import reflection


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value):
    return json.loads(json.dumps(value, sort_keys=True, allow_nan=False))


def maxabs(matrix):
    data = matrix.data if sp.issparse(matrix) else np.asarray(matrix)
    return float(np.max(np.abs(data))) if np.size(data) else 0.


def geometric_gates():
    frame = np.array([[float(a)+float(b)*np.sqrt(5) for a,b in row]
                     for row in sg.slow_geometry.raw_frame()])
    frame /= np.linalg.norm(frame[0])
    steps = np.array([frame[b]-frame[a] for a,b in sg.EDGES])
    lengths = np.linalg.norm(steps, axis=1)
    assert np.max(abs(lengths-lengths[0])) < 1e-13
    axes = steps/lengths[:,None]
    spins = np.einsum('ei,iab->eab', axes, sg.SIGMA)
    return frame, steps, spins


def word_certificates(spins):
    prep_word = reflection.plan(3)['preparation']
    actual_prep = compiler.execute(prep_word,2)
    embedded = compiler.remap(prep_word,[1,2])
    phase = np.full(8,np.pi); phase[0]=0
    ref_word = compiler.inverse(embedded)+compiler.compile_diagonal(phase)+embedded
    actual_ref = compiler.execute(ref_word,3)
    seed = np.zeros(8); seed[[0,2,4]]=1/np.sqrt(3)
    ideal_ref = 2*np.outer(seed,seed)-np.eye(8)
    spin_words=[];actual_spins=[]
    for matrix in spins:
        axis = np.array([matrix[0,1].real,-matrix[0,1].imag,matrix[0,0].real])
        theta=float(np.arccos(np.clip(axis[2],-1,1)))
        phi=float(np.arctan2(axis[1],axis[0]))
        word=[compiler.pulse('Z',phi/2,[0]),compiler.pulse('Y',theta/2,[0]),
              compiler.pulse('I',-np.pi/2,[0]),compiler.pulse('Z',np.pi/2,[0]),
              compiler.pulse('Y',-theta/2,[0]),compiler.pulse('Z',-phi/2,[0])]
        spin_words.append(word);actual_spins.append(compiler.execute(word,1))
    multiplex=[]
    for event,word in enumerate(spin_words):
        for pulse in word:
            # The event projector is the exact four-term Z expansion.
            for mask in range(4):
                label=[pulse['kind']]+['Z' if mask>>bit&1 else 'I' for bit in range(2)]
                coefficient=(-1)**((event&mask).bit_count())/4
                multiplex+=compiler.pauli_rotation(label,pulse['angle']*coefficient)
    actual_multiplex=compiler.execute(multiplex,3)
    ideal_multiplex=np.eye(8,dtype=complex)
    for event in range(3):ideal_multiplex[2*event:2*event+2,2*event:2*event+2]=spins[event]
    errors={'preparation':float(np.linalg.norm(actual_prep[:,0]-np.array([1,1,1,0])/np.sqrt(3))),
            'reflection':float(np.linalg.norm(actual_ref-ideal_ref,'fro')),
            'spin_gates':[float(np.linalg.norm(a-b,'fro')) for a,b in zip(actual_spins,spins)],
            'complete_spin_multiplex':float(np.linalg.norm(actual_multiplex-ideal_multiplex,'fro'))}
    assert max(errors['preparation'],errors['reflection'],errors['complete_spin_multiplex'],*errors['spin_gates'])<1e-9
    count = len(prep_word)+16*(len(ref_word)+len(multiplex))
    bound=reflection.complete_bound(count,Fraction(1,10**12))
    arrays={'source_preparation':actual_prep,'source_reflection':actual_ref,
            'source_spin_gates':np.array(actual_spins),'source_spin_multiplex':actual_multiplex}
    words={'preparation':prep_word,'reflection':ref_word,'spin_gates':spin_words,'spin_multiplex':multiplex}
    return seed,ideal_ref,errors,bound,words,arrays


def source_columns(q):
    states=sg.finite_slice(q); index={s:i for i,s in enumerate(states)};n=len(states)
    targets=np.zeros((3,n),int);currents=np.zeros((3,n),int)
    signs=np.zeros((3,n),int);overflows=np.zeros((3,n),bool);proposed=np.zeros((3,n),int)
    raw=[]
    for e,seam in enumerate(sg.EDGES):
        for col,s in enumerate(states):
            hop,sign,_=sg.joined(s,(0,1),seam,limit=1)
            targets[e,col]=index[hop.state]
            currents[e,col]=0 if hop.overflow else hop.proposed_number_current
            proposed[e,col]=hop.proposed_number_current
            signs[e,col]=sign;overflows[e,col]=hop.overflow
            loads=tuple(int(v in s.pair) for v in range(12))
            for coin in (0,1):
                out,record=native.native(loads,seam,coin)
                mask,_=sg.encoded(s,(0,1))
                label,target,car_sign,j=sg.car.branch_action(mask,'C',coin,sites=12,flavors=2,seam=seam)
                assert j==record[1]
                assert tuple(v for v,x in enumerate(out) if x)==tuple(v for v in range(12) if ((target>>(2*v))&3))
                raw.append((e,col,coin,j,label,target,car_sign))
            assert sg.g.gauss(s)==sg.g.gauss(hop.state)
            assert all(hop.state.flux[k]-s.flux[k]==(-q*currents[e,col] if k==sg.g.EDGE_INDEX[seam] else 0)
                       for k in range(30))
    return states,targets,currents,signs,overflows,proposed,np.array(raw,dtype=np.int64)


def select_matrix(targets,currents,signs,spins,mode='native'):
    n=targets.shape[1];N=n+1;rows=[];cols=[];values=[]
    def address(e,h,s,i):return (((e*2+h)*2+s)*N+i)
    for e in range(4):
        for h in range(2):
            for s in range(2):
                for i in range(N):
                    col=address(e,h,s,i)
                    if e==3 or i==n:
                        rows.append(col);cols.append(col);values.append(1);continue
                    j=int(currents[e,i]) if mode!='zero_current' else 0
                    out=int(targets[e,i]) if j else i
                    hh=h if j else 1-h
                    coefficient=1j*j*int(signs[e,i])*(-1)**h if j else 1
                    for ss in range(2):
                        value=coefficient*spins[e,ss,s]
                        if value:
                            rows.append(address(e,hh,ss,out));cols.append(col);values.append(value)
    return sp.coo_matrix((values,(rows,cols)),shape=(16*N,16*N)).tocsr()


def current_operators(targets,currents,signs):
    n=targets.shape[1];ops=[]
    for e in range(3):
        ops.append(sp.coo_matrix((1j*currents[e]*signs[e],(targets[e],np.arange(n))),shape=(n,n)).tocsr())
    return ops


def reference_read(initial,output,phase):
    """Full reference plus/minus and every source address; no conditioning."""
    plus=(initial+phase*output)/2;minus=(initial-phase*output)/2
    return np.array([abs(plus)**2,abs(minus)**2])


def plant_reader():
    out=[]
    for lam in (0.,.5):
        b=np.sqrt(1-lam*lam);w=np.array([[lam,-b],[b,lam]])
        initial=np.array([1.,0.]);v=initial.copy()
        for t in range(17):
            for phase in (1.,1j):
                law=reference_read(initial,v,phase)
                measured=float(law[0].sum()-law[1].sum())
                expected=float((phase*np.cos(t*np.arccos(lam))).real)
                assert abs(measured-expected)<1e-10
                out.append({'lambda':lam,'time':t,'phase':[complex(phase).real,complex(phase).imag],
                            'plus':float(law[0].sum()),'minus':float(law[1].sum()),'quadrature':measured})
            v=w@v
    assert any(abs(x['quadrature'])>.9 for x in out)
    return out


def physical_signal_basis(states):
    groups={}
    for i,s in enumerate(states):
        if not any(sg.g.gauss(s)):
            key=(s.pair,s.tagged_port,s.flux)
            groups.setdefault(key,[]).append(i)
    groups=[groups[key] for key in sorted(groups)]
    assert all(len(g)==6 for g in groups)
    n=len(states);cols=[]
    for spin in range(2):
        for group in groups:
            v=np.zeros(2*n,complex);v[np.array(group)+spin*n]=1/np.sqrt(6);cols.append(v)
    return np.array(cols).T


def run_case(q,seed,ref,spins,arrays):
    states,target,current,sign,overflow,proposed,raw=source_columns(q)
    n=len(states);N=n+1;dim=16*N
    valid=sp.diags([1]*n+[0],format='csr');invalid=sp.diags([0]*n+[1],format='csr')
    R=sp.kron(sp.csr_matrix(ref),sp.kron(sp.eye(2),valid))+sp.kron(sp.eye(16),invalid)
    S=select_matrix(target,current,sign,spins);W=S@R
    ops=current_operators(target,current,sign)
    A=sum(sp.kron(sp.csr_matrix(spins[e]),ops[e]) for e in range(3))/3
    # Signal embedding uses every valid source state and both spin states.
    erows=[];ecols=[]
    for e in range(3):
        for spin in range(2):
            erows.extend((e*4+spin)*N+np.arange(n));ecols.extend(spin*n+np.arange(n))
    embed=sp.coo_matrix((np.full(len(erows),1/np.sqrt(3)),(erows,ecols)),shape=(dim,2*n)).tocsr()
    compression=embed.getH()@(S@embed)
    errors={'select_hermiticity':maxabs(S-S.getH()),'select_involution':maxabs(S@S-sp.eye(dim)),
            'reflection_hermiticity':maxabs(R-R.getH()),'reflection_involution':maxabs(R@R-sp.eye(dim)),
            'walk_unitarity':maxabs(W.getH()@W-sp.eye(dim)),
            'signal_compression':maxabs(compression-A)}
    electric=[];current_after_ref=[]
    for e,seam in enumerate(sg.EDGES):
        values=np.array([s.flux[sg.g.EDGE_INDEX[seam]] for s in states]+[0])
        E=sp.diags(np.tile(values,16),format='csr');electric.append(E)
        qdiag=np.zeros((4,2,2,N));qdiag[e,:,:,:n]=current[e]
        current_reader=R.getH()@sp.diags(qdiag.ravel(),format='csr')@R
        current_after_ref.append(current_reader)
        errors[f'heisenberg_electric_continuity_{e}']=maxabs(W.getH()@E@W-E+q*current_reader)
    # All source Gauss sectors, including coherences, are preserved.
    for vertex in range(12):
        gvals=[sg.g.gauss(s)[vertex] for s in states]+[0]
        G=sp.diags(np.tile(gvals,16),format='csr')
        errors[f'gauss_commutator_{vertex}']=maxabs(W@G-G@W)
    assert max(errors.values())<1e-10
    arrays.update({f'q{q}_targets':target,f'q{q}_currents':current,f'q{q}_car_signs':sign,
                   f'q{q}_overflow':overflow,f'q{q}_proposed':proposed,f'q{q}_native_records':raw})
    # The emitted pair is on Gauss zero and memory character zero.
    signal=np.zeros(2*n,complex)
    for m in range(6):signal[states.index(sg.g.pair_string((0,1),q,m))]=1/np.sqrt(6)
    initial=embed@signal
    moments=[];laws_by_mode={};source_states=[];chebyshev_error=None
    for mode in ('native','zero_current','omit_reflection'):
        w=S if mode=='omit_reflection' else (select_matrix(target,current,sign,spins,mode)@R)
        v=initial.copy();laws=[];t0=signal.copy();t1=A@signal;recurrence=[]
        for t in range(17):
            pair=np.array([reference_read(initial,v,z) for z in (1,1j)])
            assert abs(pair[0].sum()-1)<1e-10 and abs(pair[1].sum()-1)<1e-10
            laws.append(pair)
            if mode=='native':
                expected=t0 if t==0 else t1
                observed=embed.getH()@v
                recurrence.append(float(np.linalg.norm(observed-expected)))
                correlation=np.vdot(initial,v)
                moments.append({'time':t,'real':float(correlation.real),'imag':float(correlation.imag),
                                'query_one_probability':float(abs(v.reshape(4,2,2,N)[:,1]).ravel()@abs(v.reshape(4,2,2,N)[:,1]).ravel()),
                                'padding_probability':float(np.vdot(v.reshape(4,2,2,N)[3],v.reshape(4,2,2,N)[3]).real)})
                source_states.append(v.copy())
                if t>=1:t0,t1=t1,2*(A@t1)-t0
            v=w@v
        laws_by_mode[mode]=np.array(laws)
        arrays[f'q{q}_{mode}_complete_reference_laws']=laws_by_mode[mode]
        if mode=='native':
            chebyshev_error=max(recurrence)
            assert chebyshev_error<1e-10
    changes={mode:float(np.max(abs(law-laws_by_mode['native']))) for mode,law in laws_by_mode.items() if mode!='native'}
    assert all(v>1e-6 for v in changes.values())
    # Every eigenstate in the fixed Gauss0 / memory0 signal subspace is kept.
    B=physical_signal_basis(states);small=B.conj().T@(A@B)
    assert np.linalg.norm(A@B-B@small)<1e-10
    eigenvalues,V=np.linalg.eigh(small)
    stationary=[];stationary_vectors=[]
    for k,lam in enumerate(eigenvalues):
        v=B@V[:,k];u=embed@v;beta=np.sqrt(max(0.,1-float(lam)**2))
        if beta>1e-10:
            orth=(S@u-lam*u)/beta;psi=(u-1j*orth)/np.sqrt(2);phase=complex(lam,beta)
        else:psi=u;phase=complex(lam)
        eigenerror=float(np.linalg.norm(W@psi-phase*psi))
        realized=[float(np.vdot(psi,K@psi).real) for K in current_after_ref]
        deltas=[float(np.vdot(W@psi,E@(W@psi)).real-np.vdot(psi,E@psi).real) for E in electric]
        coherent=[float(np.vdot(v,sp.kron(sp.eye(2),J)@v).real) for J in ops]
        assert eigenerror<1e-9
        assert max(abs(np.array(deltas)+q*np.array(realized)))<1e-9
        if q:assert max(abs(np.array(realized)))<1e-9
        stationary.append({'index':k,'lambda':float(lam),'phase_radians':float(np.arccos(np.clip(lam,-1,1))),
                           'eigenstate_residual':eigenerror,'realized_number_current':realized,
                           'electric_change':deltas,'signal_coherent_J_quadrature':coherent})
        stationary_vectors.append(psi)
    # Invalid input is a complete, flagged identity sector.
    invalid_input=np.zeros(dim,complex);invalid_input[n]=1
    assert np.array_equal(W@invalid_input,invalid_input)
    # An incorrect overflow current cannot be promoted to a query.
    mutant=select_matrix(target,proposed,sign,spins)
    defect=maxabs(mutant-mutant.getH())+maxabs(mutant@mutant-sp.eye(dim))
    if q:assert defect>1e-6
    arrays[f'q{q}_stationary_eigenvalues']=eigenvalues
    arrays[f'q{q}_stationary_walk_vectors']=np.array(stationary_vectors)
    arrays[f'q{q}_signal_basis']=B
    return {'charge_number':q,'valid_source_columns':n,'full_walk_columns':dim,
            'literal_native_calls':len(raw),'all_matrix_errors':errors,
            'source_current_query_holds':int(np.count_nonzero(current==0)),
            'overflow_columns':int(np.count_nonzero(overflow)),
            'reference_records_per_time':2*2*dim,'times':17,'moments':moments,
            'changed_law_max_record_difference':changes,'chebyshev_max_error':chebyshev_error,
            'stationary_signal_dimension':len(eigenvalues),'stationary_states':stationary,
            'proposed_current_mutant_query_defect':defect,'invalid_failure_probability':1,
            'spectator_spin':'full two-state tensor identity',
            'physical_energy_or_magnetic_moment_claimed':False}


def source_pins():
    paths=set()
    for module in tuple(sys.modules.values()):
        f=getattr(module,'__file__',None)
        if f:
            p=Path(f).resolve()
            if p.is_relative_to(RER) and p.is_file():paths.add(p)
    inherited=json.loads((HERE.parent/'spin_gauge_join/FROZEN_SOURCE_PINS.json').read_text())
    paths.update(ROOT/p for p in inherited)
    paths.update([provider_path,HERE.parent/'spin_gauge_join/MANIFEST.json',
                  HERE.parent/'spin_gauge_join/PROOF.md',HERE.parent/'spin_gauge_join/NEXT_INTERFACE.md'])
    for file in ('Lean/EventAlgebra/PauliMinimalCoupling.lean',
                 'Lean/Screen/SeamCurrentDirichletGenerator.lean',
                 'Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean',
                 'model/observer/source_scalar/compiler.py','model/native_observer_program/reflection.py',
                 'model/native_electric_gauss/PROOF.md','model/native_circulation/README.md',
                 'model/observer/holonomy_reconstruction.py','code/source_selection_model/response.json',
                 'code/m1_source_realization/receipt.json'):
        paths.add(RER/file)
    return {str(p.relative_to(ROOT)):digest(p) for p in sorted(paths)}


def freeze():
    assert not (HERE/'RECEIPT.json').exists()
    pins=source_pins();(HERE/'SOURCE_PINS.json').write_text(json.dumps(pins,sort_keys=True,indent=2)+'\n')
    files={name:digest(HERE/name) for name in ('DESIGN.md','TASK.md','probe.py','SOURCE_PINS.json')}
    (HERE/'PREREG.json').write_text(json.dumps({'schema':'oph.native_current_clock.prereg.v1','files':files,
                                              'evidence_executed':False},sort_keys=True,indent=2)+'\n')
    print('PREREG',digest(HERE/'PREREG.json'),'source_inputs',len(pins))


def check_pins():
    prereg=json.loads((HERE/'PREREG.json').read_text())
    for name,h in prereg['files'].items():assert digest(HERE/name)==h,name
    for name,h in json.loads((HERE/'SOURCE_PINS.json').read_text()).items():assert digest(ROOT/name)==h,name


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true');parser.add_argument('--record',action='store_true')
    args=parser.parse_args()
    if args.freeze:freeze();return
    check_pins()
    if args.record and (HERE/'RECEIPT.json').exists():raise FileExistsError('Existing evidence must not be overwritten')
    start=time.perf_counter();frame,steps,spins=geometric_gates()
    seed,ref,errors,bound,words,arrays=word_certificates(spins)
    print('Complete source-word matrices passed.',flush=True)
    plants=plant_reader();cases=[]
    for q in (0,1):
        cases.append(run_case(q,seed,ref,spins,arrays));print(f'Full q={q} query walk and stationary currents passed.',flush=True)
    receipt=canonical({'schema':'oph.native_current_clock.v1','candidate':'CAR-current-query-clock-v1',
              'class':'N+ observer program; no physical kinetic interpretation','word_errors':errors,
              'source_word_complete_bound':bound,'cases':cases,'plants':plants,
              'nature_measurements_used':[],'new_full_integer_query_M6_compilation':False,
              'stationary_dipole_g_larmor_established':False,'source_native_calls':sum(c['literal_native_calls'] for c in cases),
              'preregistration_sha256':digest(HERE/'PREREG.json')})
    check_pins()
    if args.record:
        (HERE/'RECEIPT.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
        (HERE/'SOURCE_WORDS.json').write_text(json.dumps(words,sort_keys=True,indent=2)+'\n')
        np.savez_compressed(HERE/'RECORDS.npz',**arrays)
    else:
        assert receipt==json.loads((HERE/'RECEIPT.json').read_text())
        stored=np.load(HERE/'RECORDS.npz');assert set(stored.files)==set(arrays)
        for k,value in arrays.items():assert np.array_equal(value,stored[k]),k
        assert canonical(words)==json.loads((HERE/'SOURCE_WORDS.json').read_text())
    print(json.dumps({'all_pass':True,'elapsed_seconds':time.perf_counter()-start,
                      'native_calls':receipt['source_native_calls'],
                      'full_walk_columns':[c['full_walk_columns'] for c in cases],
                      'stationary_signal_dimensions':[c['stationary_signal_dimension'] for c in cases],
                      'physical_rows_complete':[]},indent=2))


if __name__=='__main__':main()
