"""Independent qubit/log and full continued-instrument modular observer audit."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.linalg import expm
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(BASE/'next/C10'))
from charged_reservoir import BoundarySource


def run():
    paths=[Path(__file__),HERE/'interacting_modular_observer.py',HERE/'INTERACTING_MODULAR_OBSERVER.json',HERE/'interacting_modular_observer_states.npz',BASE/'next/C10/canonical_states.npz',BASE/'next/C10/joint_source_states.npz']
    pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    saved=json.loads((HERE/'INTERACTING_MODULAR_OBSERVER.json').read_text());z=np.load(HERE/'interacting_modular_observer_states.npz');canonical=np.load(BASE/'next/C10/canonical_states.npz');joint=np.load(BASE/'next/C10/joint_source_states.npz');h=joint['hamiltonian'];native=canonical['selected_protocol_poststate']
    indices=[row['archive_index'] for row in saved['observer']['configurations']];p=float(np.trace(native[np.ix_(indices,indices)]).real);rho=native[np.ix_(indices,indices)]/p
    tr=float(np.trace(rho).real);gap=float(np.sqrt((rho[0,0].real-rho[1,1].real)**2+4*abs(rho[0,1])**2));lo=(tr-gap)/2;hi=(tr+gap)/2
    if lo<=0:raise ValueError('notfaithful; nofloorallowed')
    frequency=float(np.log(hi/lo));logrho=np.log(lo)*np.eye(2)+(np.log(hi)-np.log(lo))/gap*(rho-lo*np.eye(2));k=-logrho
    root=np.sqrt(lo)*np.eye(2)+(np.sqrt(hi)-np.sqrt(lo))/gap*(rho-lo*np.eye(2))
    a=z['modular_A0'];s=.37;flow=np.cos(frequency*s)*a+1j*np.sin(frequency*s)*(rho@a-a@rho)/gap
    flowerror=float(np.max(abs(flow-z['modular_As'])));aerror=max(float(np.max(abs(a@a-np.eye(2)))),float(np.max(abs(a-a.conj().T))),float(np.max(abs((rho-.5*np.eye(2))@a+a@(rho-.5*np.eye(2))))))
    delta=np.kron(np.linalg.inv(rho).T,rho);omega=root.reshape(-1,order='F');u=expm(1j*s*logrho);du=np.kron(expm(-1j*s*logrho).T,u)
    gns=max(float(np.max(abs(delta-z['GNS_Delta']))),float(np.max(abs(omega-z['GNS_Omega']))),float(np.max(abs(du@np.kron(np.eye(2),a)@du.conj().T-np.kron(np.eye(2),flow)))))
    rng=np.random.default_rng(731);kms=0.;wrong=0.
    for _ in range(5):
        x=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2));y=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2));ys=u@y@u.conj().T
        rhs=np.trace(rho@ys@x);kms=max(kms,float(abs(np.trace(rho@x@rho@ys@np.linalg.inv(rho))-rhs)));wrong=max(wrong,float(abs(np.trace(rho@x@np.linalg.inv(rho)@ys@rho)-rhs)))
    model=BoundarySource(6);allN=np.array([mask.bit_count() for mask,_ in model.states]);ids=np.flatnonzero(allN==3);local=[list(ids).index(i) for i in indices];hb=h[np.ix_(ids,ids)];unitary=expm(-.37j*hb)
    effect=np.zeros((model.n,model.n),complex);posterror=0.;continuationerror=0.;ledgererror=0.;probabilityerror=0.;pureerror=0.;probabilities=[]
    for record in saved['complete_observer_records']:
        if record.get('status')=='outside_corner':
            branch=native.copy();branch[indices,:]=0;branch[:,indices]=0
            q=np.ones(model.n);q[indices]=0;effect+=np.diag(q)
        else:
            first=record['first_A0'];last=record['second_As'];pa=(np.eye(2)+first*a)/2;pb=(np.eye(2)+last*flow)/2;small=pb@pa
            branch=np.zeros_like(native,dtype=complex);branch[np.ix_(indices,indices)]=p*small@rho@small.conj().T;effect[np.ix_(indices,indices)]+=small.conj().T@small
            expected=p*(1+first*last*np.cos(frequency*s))/4;probabilityerror=max(probabilityerror,abs(expected-record['probability']))
            conditional=branch[np.ix_(indices,indices)]/np.trace(branch).real;pureerror=max(pureerror,float(abs(np.linalg.det(conditional))))
        posterror=max(posterror,float(np.max(abs(branch-z[record['poststate_archive_key']]))));probabilities.append(float(np.trace(branch).real))
        coherent=unitary@branch[np.ix_(ids,ids)]@unitary.conj().T;continued=.71*coherent;ledger=[]
        for operations in model.operators['MIX']:
            for _,fullk in operations:
                kk=fullk[np.ix_(ids,ids)];part=.29*kk@coherent@kk.T/3;continued+=part;ledger.append(float(np.trace(part).real))
        ledger.append(.71*float(np.trace(branch).real));ledgererror=max(ledgererror,float(np.max(abs(np.array(ledger)-[row['joint_probability'] for row in record['complete_native_continuation']]))))
        continuationerror=max(continuationerror,float(np.max(abs(continued-z[record['continued_archive_key']][np.ix_(ids,ids)]))))
    # Reproduce the declared after-inspection selection rule independently.
    ancestry=saved['observer']['selection_ancestry'];active=np.flatnonzero(np.diag(native).real>ancestry['active_diagonal_threshold']);selected=None
    for offset,i in enumerate(active):
        for j in active[offset+1:]:
            block=native[np.ix_([i,j],[i,j])];conditional=block/np.trace(block).real
            trace=float(np.trace(conditional).real);diff=float(np.sqrt((conditional[0,0].real-conditional[1,1].real)**2+4*abs(conditional[0,1])**2))
            if (trace-diff)/2>ancestry['conditional_minimum_threshold'] and abs(native[i,j])>ancestry['offdiagonal_threshold']:selected=[int(i),int(j)];break
        if selected is not None:break
    conditioned=np.zeros((len(ids),)*2,complex);conditioned[np.ix_(local,local)]=rho;evolved=unitary@conditioned@unitary.conj().T;leakage=float(np.trace(evolved).real-np.trace(evolved[np.ix_(local,local)]).real)
    hc=h[np.ix_(indices,indices)];hc0=hc-np.trace(hc)*np.eye(2)/2;k0=k-np.trace(k)*np.eye(2)/2;alpha=float(np.vdot(k0,hc0).real/np.vdot(k0,k0).real);fit=float(np.linalg.norm(hc0-alpha*k0))
    correlation=np.cos(frequency*s)+1j*(lo-hi)*np.sin(frequency*s)
    configuration_ok=all(model.states[row['archive_index']]==(row['occupation_mask'],tuple(row['half_record_residues'])) for row in saved['observer']['configurations'])
    eigen,basis=np.linalg.eigh(rho)
    for j in range(2):basis[:,j]*=np.exp(-1j*np.angle(basis[int(np.argmax(abs(basis[:,j]))),j]))
    frame_error=float(np.max(abs(a-basis@np.array([[0,1],[1,0]])@basis.conj().T)))
    checks={'explicit_global_configuration_labels':configuration_ok,'declared_eigenframe_phase':frame_error<1e-12,'native_input_exact_C10_output':bool(np.array_equal(native,z['native_input'])),'qubit_eigenvalues_without_floor':bool(np.allclose([lo,hi],saved['observer']['conditional_eigenvalues'],atol=1e-12)),'frozen_family_qubit_formula':flowerror<1e-12 and aerror<1e-12,'independent_log_GNS_and_KMS':gns<1e-11 and kms<1e-11 and wrong>.1,'every_PQ_sequential_effect_complete':float(np.max(abs(effect-np.eye(model.n))))<1e-12,'every_poststate_and_probability':posterror<1e-12 and probabilityerror<1e-12,'every_subnormalized_native_continuation':continuationerror<1e-12 and ledgererror<1e-12,'pure_branch_has_no_faithful_M2_flow':pureerror<1e-12,'after_inspection_corner_selection_reproduced':selected==indices,'native_H_clock_identification_rejected':leakage>.01 and fit>.01,'full_M384_faithfulness_structurally_impossible':1+sum(len(ops) for ops in model.operators['MIX'])<model.n}
    checks={key:bool(value) for key,value in checks.items()}
    if not all(checks.values()):raise ValueError(checks)
    if pins!={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}:raise ValueError('sourcechanged')
    return {'schema':'oph.c11.interacting_modular_independent_review.v1','checks':checks,'corner_probability':p,'eigenvalues':[lo,hi],'frequency':frequency,'independent_complex_correlation':{'real':float(correlation.real),'imag':float(correlation.imag)},'errors':{'family':flowerror,'GNS':gns,'KMS':kms,'poststates':posterror,'probabilities':probabilityerror,'native_continuation':continuationerror,'joint_native_probabilities':ledgererror},'H_clock_control':{'leakage':leakage,'compressed_fit_residual':fit},'scope':'Postselectedglobaltwo-configurationM2corner, notlocalsubsystem orphysicalobserverclock. Cornerselectedafterknown-outputinspection, frozenbeforemeasurements. Fullcomplement/subnormalizedbranches retained. GNSreference algebraic; modulars injected andnotnativephysicaltime. RealPVMcorrelation only; no subjectiveexperience or physicalrowaward.','source_sha256':pins,'full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('INTERACTING_MODULAR_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'checks':r['checks'],'errors':r['errors'],'frequency':r['frequency']}))
