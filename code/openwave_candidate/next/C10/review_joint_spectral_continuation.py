"""Second independent Codex-worker audit using archived H and block projectors.

No producer constructors or spectral instrument implementation are imported.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import eigh
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]

def require(ok,message):
    if not ok:raise ValueError(message)

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def direct_native(masks,records,indices):
    """Direct endpoint balanced occupation + intervening-bit CAR parity."""
    count=len(indices);lookup={(int(masks[i]),tuple(records[i])):j for j,i in enumerate(indices)}
    out={}
    for seam,(a,b) in enumerate(((0,1),(1,2),(2,0))):
        for law in ('C','E'):
            for coin in (0,1):
                alphabet=((0,1) if coin==0 else (-1,0)) if law=='C' else (-1,0,1)
                for value in alphabet:out[(seam,law,coin,value)]=np.zeros((count,count))
                for col,i in enumerate(indices):
                    mask=int(masks[i]);na=(mask>>a)&1;nb=(mask>>b)&1
                    target_a=(na+nb)//2 if coin==0 else (na+nb+1)//2
                    current=na-target_a;target=mask;sign=1
                    if current:
                        target^=(1<<a)|(1<<b)
                        sign=(-1)**sum((mask>>z)&1 for z in range(min(a,b)+1,max(a,b)))
                    f=records[i].copy();f[2*seam]=(f[2*seam]+current)%6;f[2*seam+1]=(f[2*seam+1]-current)%6
                    key=(seam,law,coin,current if law=='C' else na-nb)
                    out[key][lookup[(target,tuple(f))],col]=sign/np.sqrt(12)
    return out

def run():
    names=['JOINT_SOURCE.json','joint_source_states.npz','CANONICAL_CONSTRUCTION.json','canonical_states.npz','CANONICAL_SOURCE_CONTRACT.json','joint_source.py','charged_reservoir.py','modular_field_adapter.py','spectral_readout.py','canonical_construction.py']
    paths=[Path(__file__),*[HERE/name for name in names]];before={str(p.relative_to(BASE)):digest(p) for p in paths}
    joint=json.loads((HERE/'JOINT_SOURCE.json').read_text());canonical=json.loads((HERE/'CANONICAL_CONSTRUCTION.json').read_text());contract=json.loads((HERE/'CANONICAL_SOURCE_CONTRACT.json').read_text())
    require(digest(HERE/joint['state_archive']['filename'])==joint['state_archive']['sha256'],'joint archive custody')
    for receipt in (joint,canonical,contract):
        for path,pin in receipt['source_sha256'].items():require(digest(BASE/path)==pin,'source custody '+path)
    with np.load(HERE/'joint_source_states.npz',allow_pickle=False) as z:inputs={k:z[k] for k in z.files}
    with np.load(HERE/'canonical_states.npz',allow_pickle=False) as z:archive={k:z[k] for k in z.files}
    h=inputs['hamiltonian'];masks=inputs['occupation_masks'];records=inputs['half_records'];numbers=np.array([int(x).bit_count() for x in masks]);n=len(h)
    require(n==384 and np.array_equal(inputs['boundary_B'],3-numbers),'wrong joined source')
    spectral={};ground=(float('inf'),None,None);projector_error=0.;coefficient_error=0.;prob_error=0.;rows_checked=0
    for number in range(7):
        ids=np.flatnonzero(numbers==number);values,vectors=eigh(h[np.ix_(ids,ids)],driver='evr')
        spectral[number]=(ids,values,vectors)
        if values[0]<ground[0]:ground=(float(values[0]),number,vectors[:,0])
    e0,q0,gvec=ground;ids0=spectral[q0][0];rho=inputs['global_allN_reference'];expected=np.outer(gvec,gvec.conj())
    ref_error=float(np.max(abs(rho[np.ix_(ids0,ids0)]-expected)))
    require(q0==3 and abs(e0+2.7262233938630343)<1e-11 and ref_error<1e-11,'interacting reference mismatch')
    instrument=canonical['spectral_readout_instrument'];all_inputs={k:v for k,v in inputs.items() if k!='hamiltonian' and v.shape==h.shape};all_inputs['selected_protocol_poststate']=archive['selected_protocol_poststate']
    # Independent EVR projectors, rather than copying the producer eigenvectors.
    for label,rows in instrument['complete_branch_records'].items():
        if label not in all_inputs:continue
        density=all_inputs[label];total=0.
        for row in rows:
            ids,values,vectors=spectral[row['source_number']];take=abs(values-row['energy'])<1e-8
            require(int(take.sum())==row['degeneracy'],'degeneracy mismatch')
            independent=vectors[:,take];p=independent@independent.conj().T
            basis=archive[row['basis_archive_key']];coefficient=archive[row['coefficient_archive_key']]
            projector_error=max(projector_error,float(np.max(abs(p-basis@basis.conj().T))))
            sub=density[np.ix_(ids,ids)];projected=independent@(independent.conj().T@sub@independent)@independent.conj().T
            compressed=basis@coefficient@basis.conj().T
            coefficient_error=max(coefficient_error,float(np.max(abs(projected-compressed))))
            total+=row['probability'];prob_error=max(prob_error,abs(float(np.trace(projected).real)-row['probability']))
            require(np.linalg.eigvalsh(coefficient).min()>-1e-10,'nonpositive subnormalized coefficient')
            rows_checked+=1
        prob_error=max(prob_error,abs(total-float(np.trace(density).real)))
    # Reproduce canonical native protocol directly from occupation/record tables.
    native=direct_native(masks,records,ids0);rho0=rho[np.ix_(ids0,ids0)]
    commit=sum((k@rho0@k.conj().T for k in native.values()),np.zeros_like(rho0))
    selected=.71*rho0+.29*commit
    selected_error=float(np.max(abs(selected-archive['selected_protocol_poststate'][np.ix_(ids0,ids0)])))
    completeness=float(np.max(abs(sum((k.conj().T@k for k in native.values()),np.zeros_like(rho0))-np.eye(len(ids0)))))
    # Actual populated measured branch: full subnormalization and all native events.
    actual=instrument['actual_populated_degenerate_continuation'];row=actual['branch'];ids,values,vectors=spectral[row['source_number']]
    take=abs(values-row['energy'])<1e-8;p=vectors[:,take]@vectors[:,take].conj().T;density=all_inputs[actual['input_archive_label']][np.ix_(ids,ids)];branch=p@density@p
    unitary=(vectors*np.exp(-.37j*values))@vectors.conj().T;coherent=unitary@branch@unitary.conj().T
    operators=direct_native(masks,records,ids);continued=.71*coherent;ledger={}
    for key,k in operators.items():
        part=.29*k@coherent@k.conj().T;continued+=part;ledger[key]=float(np.trace(part).real)
    ledger[None]=.71*float(np.trace(branch).real)
    saved=archive[actual['poststate_archive_key']][np.ix_(ids,ids)]
    continuation_error=float(np.max(abs(continued-saved)));ledger_error=0.
    for record in actual['complete_joint_native_outcomes']:
        key=None if record.get('status')=='no_commit' else tuple(record[x] for x in ('seam','law','coin','record_value'))
        ledger_error=max(ledger_error,abs(ledger.pop(key)-record['joint_probability']))
    require(not ledger,'incomplete native output alphabet')
    # Complex control in an actual doublet; invariance under a different unitary.
    control=instrument['fine_basis_rejecting_control'];energy=control['energy'];ids,values,vectors=spectral[3];basis=vectors[:,abs(values-energy)<1e-8]
    c=np.array([[.5,.3+.2j],[.3-.2j,.5]]);state=basis@c@basis.conj().T
    rotation=np.array([[1,1],[-1,1]],complex)/np.sqrt(2);rotation=np.diag([np.exp(.37j),np.exp(-.21j)])@rotation
    b=basis@rotation;rotated=b@(rotation.conj().T@c@rotation)@b.conj().T
    rotation_error=float(np.max(abs(rotated-state)));fine=basis@np.diag(np.diag(c))@basis.conj().T
    fine_defect=float(np.linalg.norm(state-fine,ord='nuc')/2)
    native3=direct_native(masks,records,ids)
    channel=lambda x:.71*x+.29*sum((k@x@k.conj().T for k in native3.values()),np.zeros_like(x))
    continued_fine_defect=float(np.linalg.norm(channel(state)-channel(fine),ord='nuc')/2)
    errors={'reference':ref_error,'independent_projectors':projector_error,'full_branch_coefficients':coefficient_error,'complete_probabilities':prob_error,'canonical_native_protocol':selected_error,'native_completeness':completeness,'actual_subnormalized_continuation':continuation_error,'all_native_joint_events':ledger_error,'complex_basis_rotation':rotation_error}
    require(max(errors.values())<1e-9 and fine_defect>.3 and continued_fine_defect>.1,'joint/spectral audit failed')
    require(before=={str(p.relative_to(BASE)):digest(p) for p in paths},'source snapshot changed during review')
    return {'schema':'oph.c10.joint_spectral.second_independent_review.v1','reviewer':'Independent Codex worker distinct from producer and first reviewer; not opposite-team Fable audit','outcome':'PASS_WITH_FINITE_DECLARED_SCOPE','source_dimension':n,'recomputed_interacting_ground':e0,'spectral_rows_checked':rows_checked,'errors':errors,'independent_complex_control':{'preparation':'supplied faithful complex density within an actual N3 doublet; not substituted for canonical reference','fine_forget_trace_distance':fine_defect,'after_native_channel_fine_forget_distance':continued_fine_defect},'scope':'Archived actual interacting H/reference; independent SciPy EVR block projectors and direct occupation/parity native Kraus reconstruction. Complete grouped number/energy PVM, including zeros and unnormalized continuation. No full source reconstruction or producer rerun; floating checks without interval enclosure. No physical species/vacuum/mass/radiation or official criterion award.','source_sha256':before,'full_criterion_met':False,'official_rows_changed':[]}
if __name__=='__main__':
    result=run();Path(__file__).with_name('JOINT_SPECTRAL_CONTINUATION_INDEPENDENT_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
