"""Independent actual-state modular/PVM/continuation review, no old flow output."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.linalg import expm
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1];C10=BASE/'next/C10'
sys.path.insert(0,str(C10))
from joint_source import JointNativeSource


def execute():
    path=HERE/'INTERACTING_MODULAR_OBSERVER.json';r=json.loads(path.read_text());archive=HERE/r['state_archive']['filename']
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==r['state_archive']['sha256']
    for relative,pin in r['source_sha256'].items():assert hashlib.sha256((BASE/relative).read_bytes()).hexdigest()==pin
    with np.load(C10/'canonical_states.npz') as data:native=data['selected_protocol_poststate']
    with np.load(C10/'joint_source_states.npz') as data:h=data['hamiltonian'];masks=data['occupation_masks'];words=data['half_records']
    with np.load(archive) as data:z={key:data[key] for key in data.files}
    configs=[(7,(5,1,5,1,5,1)),(11,(1,0,0,0,1,5))]
    indices=[next(i for i,(mask,f) in enumerate(zip(masks,words)) if int(mask)==a and tuple(f)==b) for a,b in configs]
    p=np.zeros_like(h);p[indices,indices]=1;q=np.eye(len(h))-p;prob=float(np.trace(p@native).real);rho=native[np.ix_(indices,indices)]/prob
    determinant=float(np.linalg.det(rho).real);gap=np.sqrt(float(np.trace(rho).real)**2-4*determinant)
    eigen=np.array([(np.trace(rho).real-gap)/2,(np.trace(rho).real+gap)/2]);assert min(eigen)>0
    low=(rho-eigen[1]*np.eye(2))/(eigen[0]-eigen[1]);projectors=[low,np.eye(2)-low]
    def power(exponent):return sum((np.exp(exponent*np.log(value))*proj for value,proj in zip(eigen,projectors)),np.zeros((2,2),complex))
    def flow(a,t):return power(1j*t)@a@power(-1j*t)
    a=z['modular_A0'];s=r['modular_observable']['s'];as_=flow(a,s)
    as_error=float(np.max(abs(as_-z['modular_As'])));spectrum_error=float(max(abs(eigen-r['observer']['conditional_eigenvalues'])))
    # Independent complex AB and analytic continuation, not producer's KMS sample.
    first=np.array([[.3+.7j,-1+.2j],[.8-.3j,.4]]);second=np.array([[1-.2j,.9j],[.2+.4j,-.6]])
    kms=float(abs(np.trace(rho@first@flow(second,s-1j))-np.trace(rho@flow(second,s)@first)))
    wrong=float(abs(np.trace(rho@first@flow(second,s+1j))-np.trace(rho@flow(second,s)@first)))
    delta=np.kron(np.linalg.inv(rho).T,rho);du=np.kron(power(1j*s).conj(),power(1j*s))
    gns=float(np.max(abs(du@np.kron(np.eye(2),a)@du.conj().T-np.kron(np.eye(2),as_))))
    omega=power(.5).reshape(-1,order='F');cyclic=float(np.linalg.norm(delta@omega-omega))
    numbers=np.array([int(mask).bit_count() for mask in masks]);ids=np.flatnonzero(numbers==3)
    assert np.max(abs(native[numbers!=3]))<1e-12 # exact conserved-number support implies nonfaithfulness on M384.
    model=JointNativeSource();assert np.array_equal(h,model.hamiltonian);u=expm(-.37j*h[np.ix_(ids,ids)])
    expected=[]
    for first_label in (-1,1):
        for last_label in (-1,1):
            small=((np.eye(2)+last_label*as_)/2)@((np.eye(2)+first_label*a)/2)
            op=np.zeros_like(native,complex);op[np.ix_(indices,indices)]=small;expected.append(op@native@op.conj().T)
    expected.append(q@native@q);post_error=0.;continuation_error=0.;joint_error=0.;signed=0.
    premeasurement_energy=float(np.trace(h@native).real);postmeasurement_energy=0.;continued_energy=0.
    for record,branch in zip(r['complete_observer_records'],expected):
        saved=z[record['poststate_archive_key']];post_error=max(post_error,float(np.max(abs(saved-branch))))
        postmeasurement_energy+=float(np.trace(h@branch).real)
        b=branch[np.ix_(ids,ids)];coherent=u@b@u.conj().T;out=.71*coherent
        probabilities={('no_commit',):.71*float(np.trace(branch).real)}
        for ops in model.field.operators_by_law['MIX']:
            for meta,k in ops:
                k=k[np.ix_(ids,ids)];out+=.29*k@coherent@k.conj().T/3
                probabilities[tuple(meta.values())]=float(.29*np.trace(k@coherent@k.conj().T).real/3)
        continued=z[record['continued_archive_key']];continuation_error=max(continuation_error,float(np.max(abs(out-continued[np.ix_(ids,ids)]))),float(np.max(abs(continued[numbers!=3,:]))))
        continued_energy+=float(np.trace(h@continued).real)
        for entry in record['complete_native_continuation']:
            key=('no_commit',) if entry.get('status')=='no_commit' else tuple(entry[k] for k in ('seam','law','coin','record_value'))
            joint_error=max(joint_error,abs(probabilities[key]-entry['joint_probability']))
        if record.get('status')!='outside_corner':signed+=record['first_A0']*record['second_As']*float(np.trace(branch).real)
    correlation=np.trace(rho@a@as_);correlation_error=float(abs(signed/prob-correlation.real))
    conditional=(p@native@p)[np.ix_(ids,ids)]/prob;evolved=u@conditional@u.conj().T;local=[list(ids).index(i) for i in indices]
    leakage=float(np.trace(evolved).real-np.trace(evolved[np.ix_(local,local)]).real)
    maximum=max(as_error,spectrum_error,kms,gns,cyclic,post_error,continuation_error,joint_error,correlation_error)
    assert maximum<1e-10 and wrong>.1 and leakage>.01 and prob<.1
    files=[Path(__file__),HERE/'interacting_modular_observer.py',path,archive,C10/'CANONICAL_CONSTRUCTION.json',C10/'canonical_states.npz',C10/'JOINT_SOURCE.json',C10/'joint_source_states.npz']
    return {'schema':'oph.c11.interacting_modular_reference_review.v1','reviewer':'Independent Codex worker; no producer/provider math import','outcome':'PASS_WITH_DECLARED_CONDITIONAL_OBSERVER_SCOPE','checks':{'actual_C10_poststate_corner_reconstructed':True,'faithful_M2_without_floor':bool(min(eigen)>0),'full_M384_nonfaithful_by_number_support':True,'polynomial_modular_flow_and_GNS':bool(max(as_error,gns,cyclic)<1e-10),'independent_KMS_minus_i':kms<1e-10,'wrong_plus_i_rejected':wrong>.1,'allfive_actual_PVM_poststates':post_error<1e-10,'allnative_continued_poststates':continuation_error<1e-10,'alljoint_record_probabilities':joint_error<1e-10,'real_sequential_not_complex_probability':correlation_error<1e-10,'autonomous_native_corner_flow_rejected':leakage>.01},'independent_method':'Reconstruct source configurations from actual archive; closed2x2trace/determinant eigenvalues and polynomial projectors; independent complexAB KMS; direct five Kraus maps; exact invariantN3 block of archivedH and explicit every nativeK for unnormalized continuation.','corner_probability':prob,'conditional_eigenvalues':eigen.tolist(),'dimensionless_modular_frequency':float(np.log(eigen[1]/eigen[0])),'errors':{'modular_observable':as_error,'conditional_eigenvalues':spectrum_error,'KMS_minus_i':kms,'KMS_wrong_plus_i':wrong,'GNS':gns,'GNS_cyclic':cyclic,'five_poststates':post_error,'native_continuation':continuation_error,'complete_joint_records':joint_error,'real_correlation':float(correlation_error)},'native_clock_leakage':leakage,'additional_work_accounting':{'pre_observer_energy':premeasurement_energy,'post_observer_energy':postmeasurement_energy,'observer_measurement_mean_work':postmeasurement_energy-premeasurement_energy,'post_native_continuation_energy':continued_energy,'subsequent_native_mean_work':continued_energy-postmeasurement_energy,'scope':'Additional ideal observer measurement can do work; no energy conservation or physicalclock interpretation asserted.'},'scope':['Known-output selected coherent two-configuration corner is added apparatus, not unique/native observer selection.','Fullsource is nonfaithful; conditionalM2faithfulness cannot be silently reset after purePVM outcomes.','Modulars is injected dimensionless observableparameter, not nativeHtime or subjectiveexperience.','Complex algebraic correlation is distinct from real sequentialPVM outcome correlation.','Fullcomplement and unnormalized continuation retained; no oldone-load source result is transferred.'],'full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
if __name__=='__main__':
 r=execute();HERE.joinpath('INTERACTING_MODULAR_REFERENCE_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
