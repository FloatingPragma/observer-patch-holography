"""Independent grouped energy projection, rotation and native continuation audit."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(HERE))
from charged_reservoir import BoundarySource


def run():
    archive=np.load(HERE/'joint_source_states.npz');h=archive['hamiltonian'];m=BoundarySource(6);ids=np.flatnonzero([mask.bit_count()==3 for mask,_ in m.states]);smallh=h[np.ix_(ids,ids)];values,vectors=np.linalg.eigh(smallh)
    groups=[]
    for j,e in enumerate(values):
        if groups and abs(e-groups[-1][0])<1e-9:groups[-1][1].append(j)
        else:groups.append([float(e),[j]])
    doublet=next(indices for _,indices in groups if len(indices)==2);v=vectors[:,doublet]
    coefficient=np.array([1,1j])/np.sqrt(2);c=np.outer(coefficient,coefficient.conj());state=v@c@v.conj().T
    theta=.43;phi=.71;u=np.array([[np.cos(theta),-np.sin(theta)*np.exp(-1j*phi)],[np.sin(theta)*np.exp(1j*phi),np.cos(theta)]])
    rotated=(v@u)@(u.conj().T@c@u)@(v@u).conj().T;rotationerror=float(np.max(abs(state-rotated)))
    fine=v@np.diag(np.diag(c))@v.conj().T;fine_distance=float(np.linalg.norm(state-fine,ord='nuc')/2)
    native=np.zeros_like(state);nativefine=np.zeros_like(state);probabilities=[]
    for ops in m.operators['MIX']:
        for _,fullk in ops:
            k=fullk[np.ix_(ids,ids)]/np.sqrt(3);branch=k@state@k.T;native+=branch;nativefine+=k@fine@k.T;probabilities.append(float(np.trace(branch).real))
    # Explicit entangled reference with nonzero off-diagonal marginal.
    psi=np.zeros((len(ids),2),complex);psi[:,0]=(vectors[:,0]+v[:,0])/2;psi[:,1]=(1j*v[:,0]+v[:,1])/2
    before=psi.conj().T@psi;after=np.zeros((2,2),complex);total=0.;referencepost=np.zeros((2*len(ids),)*2,complex)
    for _,indices in groups:
        basis=vectors[:,indices];coeff=basis.conj().T@psi;branch=basis@coeff
        total+=float(np.sum(abs(coeff)**2));after+=branch.conj().T@branch
        flat=branch.reshape(-1);referencepost+=np.outer(flat,flat.conj())
    checks={'degenerate_basis_rotation_invariant':rotationerror<1e-12,'fine_basis_forget_counterexample_retained':abs(fine_distance-.5)<1e-12,
      'complete_subsequent_native_outcomes':abs(sum(probabilities)-1)<1e-12,
      'subsequent_native_preserves_distinction':float(np.linalg.norm(native-nativefine,ord='nuc')/2)>.01,
      'all_group_probabilities_complete':abs(total-1)<1e-12,
      'external_reference_not_measured':float(np.max(abs(before-after)))<1e-12 and abs(before[0,1])>.2,
      'external_reference_joint_positive':float(np.linalg.eigvalsh(referencepost).min())>-1e-12}
    checks={key:bool(value) for key,value in checks.items()}
    if not all(checks.values()):raise ValueError(checks)
    producer=json.loads((HERE/'CANONICAL_CONSTRUCTION.json').read_text())['spectral_readout_instrument']
    compact=np.load(HERE/'canonical_states.npz');joint=np.load(HERE/'joint_source_states.npz')
    coefficient_error=0.;probability_error=0.;spectral_error=0.;resolution_error=0.
    sample=next(iter(producer['complete_branch_records'].values()))
    for number in range(7):
        bases=[compact[row['basis_archive_key']] for row in sample if row['source_number']==number]
        fullbasis=np.concatenate(bases,axis=1);resolution_error=max(resolution_error,float(np.max(abs(fullbasis@fullbasis.conj().T-np.eye(len(fullbasis))))))
    for label,rows in producer['complete_branch_records'].items():
        if label in joint.files:inputstate=joint[label]
        elif label in compact.files:inputstate=compact[label]
        else:
            populated=next(row for row in rows if row['probability']>.99);indices=compact[populated['position_indices_archive_key']];basis=compact[populated['basis_archive_key']]
            vector=basis@np.ones(basis.shape[1])/np.sqrt(basis.shape[1]);inputstate=np.zeros_like(h,dtype=complex);inputstate[np.ix_(indices,indices)]=np.outer(vector,vector.conj())
        for row in rows:
            indices=compact[row['position_indices_archive_key']];basis=compact[row['basis_archive_key']];c_saved=compact[row['coefficient_archive_key']]
            c_expected=basis.conj().T@inputstate[np.ix_(indices,indices)]@basis
            coefficient_error=max(coefficient_error,float(np.max(abs(c_saved-c_expected))));probability_error=max(probability_error,abs(np.trace(c_expected).real-row['probability']))
            spectral_error=max(spectral_error,float(np.max(abs(h[np.ix_(indices,indices)]@basis-row['energy']*basis))))
    continuation=producer['actual_populated_degenerate_continuation'];row=continuation['branch'];indices=compact[row['position_indices_archive_key']];basis=compact[row['basis_archive_key']];coefficient=compact[row['coefficient_archive_key']]
    branch=basis@coefficient@basis.conj().T
    from scipy.linalg import expm
    unitary=expm(-.37j*h[np.ix_(indices,indices)]);coherent=unitary@branch@unitary.conj().T;out=.71*coherent;prob=[]
    for ops in m.operators['MIX']:
        for _,k in ops:
            small=k[np.ix_(indices,indices)];part=.29*small@coherent@small.T/3;out+=part;prob.append(float(np.trace(part).real))
    prob.append(.71*float(np.trace(branch).real));target=compact[continuation['poststate_archive_key']][np.ix_(indices,indices)]
    continuation_error=float(np.max(abs(out-target)));joint_probability_error=float(np.max(abs(np.array(prob)-[r['joint_probability'] for r in continuation['complete_joint_native_outcomes']])))
    from spectral_readout import GroupedEnergyInstrument
    from types import SimpleNamespace
    rejected=0
    for tolerance in (float('nan'),-1,0,True):
        try:GroupedEnergyInstrument(SimpleNamespace(spectral={}),tolerance)
        except ValueError:rejected+=1
    checks.update({'producer_PVM_resolution_complete':resolution_error<1e-12,'producer_all_compact_coefficients_reproduced':coefficient_error<1e-12,'producer_all_terminal_probabilities':probability_error<1e-12,'producer_spectral_bases_actual_H':spectral_error<1e-9,'actual_subnormalized_continuation':continuation_error<1e-12 and joint_probability_error<1e-12,'bad_group_tolerances_rejected':rejected==4})
    checks={key:bool(value) for key,value in checks.items()}
    if not all(checks.values()):raise ValueError(checks)
    paths=[Path(__file__),HERE/'spectral_readout.py',HERE/'CANONICAL_CONSTRUCTION.json',HERE/'canonical_states.npz',HERE/'joint_source.py',HERE/'joint_source_states.npz',HERE/'charged_reservoir.py']
    return {'schema':'oph.c10.luders_independent_review.v1','checks':checks,'basis_rotation_error':rotationerror,'fine_forget_distance':fine_distance,'post_native_fine_forget_distance':float(np.linalg.norm(native-nativefine,ord='nuc')/2),'actual_doublet_indices':doublet,'actual_energy':float(values[doublet[0]]),'producer_compared':True,'compact_coefficient_error':coefficient_error,'terminal_probability_error':float(probability_error),'spectral_error':spectral_error,'projector_resolution_error':resolution_error,'actual_continuation_error':continuation_error,'joint_native_probability_error':joint_probability_error,'scope':'IndependentactualgroupedLuders/reference/nativecontinuation controls; producer compactadapter compared on every branch. No physicalapparatussynthesis or physicalrow award.','source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('LUDERS_CONTINUATION_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
