"""Conditional modular observer on the actual interacting C10 native output.

Reuses existing finite modular-flow/GNS conventions and native continuation.
"""
from pathlib import Path
import argparse,json,hashlib,importlib.util,sys
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1];C10=BASE/'next/C10'
sys.path.insert(0,str(C10))
from joint_source import JointNativeSource
from spectral_readout import GroupedEnergyInstrument
PROVIDER=BASE/'validation/observables/modular_time/run_modular_time.py'
spec=importlib.util.spec_from_file_location('existing_conditional_modular_math',PROVIDER)
modular=importlib.util.module_from_spec(spec);spec.loader.exec_module(modular)


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def require(ok,message):
    if not ok:raise ValueError(message)

def execute():
    paths=[Path(__file__),PROVIDER,C10/'CANONICAL_CONSTRUCTION.json',C10/'CANONICAL_SOURCE_CONTRACT.json',C10/'canonical_states.npz',C10/'JOINT_SOURCE.json',C10/'joint_source_states.npz',C10/'joint_source.py',C10/'charged_reservoir.py',C10/'modular_field_adapter.py',C10/'spectral_readout.py']
    paths += [BASE/rel for rel in ['next/C9/fermionic_shared_records.py','candidate/oph_openwave/construction.py','candidate/oph_openwave/fiber_preparation.py','candidate/oph_openwave/uniform_native.py','candidate/oph_openwave/modular_native.py','validation/bridges/providers/fock.py']]
    before={str(p.relative_to(BASE)):digest(p) for p in paths}
    contract=json.loads((C10/'CANONICAL_SOURCE_CONTRACT.json').read_text())
    for rel,pin in contract['source_sha256'].items():require(digest(BASE/rel)==pin,'canonical source custody '+rel)
    with np.load(C10/'canonical_states.npz',allow_pickle=False) as z:native=z['selected_protocol_poststate']
    with np.load(C10/'joint_source_states.npz',allow_pickle=False) as z:h=z['hamiltonian'];masks=z['occupation_masks'];records=z['half_records'];reference=z['global_allN_reference']
    model=JointNativeSource();require(np.array_equal(h,model.hamiltonian),'not the selected interacting H')
    configurations=[(7,(5,1,5,1,5,1)),(11,(1,0,0,0,1,5))]
    indices=[next(i for i,(mask,f) in enumerate(zip(masks,records)) if int(mask)==m and tuple(f)==word) for m,word in configurations]
    n=model.n;P=np.zeros((n,n));P[indices,indices]=1;Q=np.eye(n)-P
    inside=P@native@P;outside=Q@native@Q;prob=float(np.trace(inside).real);rho=inside[np.ix_(indices,indices)]/prob
    values,vectors=modular.faithful_spectrum(rho)
    for j in range(2):
        pivot=int(np.argmax(abs(vectors[:,j])));vectors[:,j]*=np.exp(-1j*np.angle(vectors[pivot,j]))
    A=vectors@np.array([[0.,1.],[1.,0.]])@vectors.conj().T;s=.37
    As=modular.modular_flow(rho,A,s);frequency=float(np.log(values[1]/values[0]))
    K=-(vectors*np.log(values))@vectors.conj().T
    delta=np.kron(np.linalg.inv(rho).T,rho);dv,db=np.linalg.eigh(delta);du=(db*np.exp(1j*s*np.log(dv)))@db.conj().T
    omega=((vectors*np.sqrt(values))@vectors.conj().T).reshape(-1,order='F')
    gns_error=float(np.max(abs(du@np.kron(np.eye(2),A)@du.conj().T-np.kron(np.eye(2),As))))
    cyclic_error=float(np.linalg.norm(delta@omega-omega))
    a=np.array([[1+.2j,2-1j],[.4+.7j,-.6j]]);b=np.array([[.2,1+.8j],[-.9+.1j,1.3]])
    bs=modular.modular_flow(rho,b,s);rhs=np.trace(rho@bs@a)
    kms=float(abs(np.trace(rho@a@rho@bs@np.linalg.inv(rho))-rhs))
    wrong_kms=float(abs(np.trace(rho@a@np.linalg.inv(rho)@bs@rho)-rhs))
    # Conditioning and a separate two-time projective observable apparatus.
    # The modular family is frozen from PREmeasurement rho, not reset to a pure branch.
    arrays={'native_input':native,'corner_projector':P,'complement_projector':Q,'inside_unnormalized':inside,'outside_unnormalized':outside,'conditional_M2_density':rho,'modular_A0':A,'modular_As':As,'GNS_Delta':delta,'GNS_Omega':omega}
    apparatus=[];effects=[];continuation=GroupedEnergyInstrument(model)
    signed=0.;number=0;pure_rejections=0;complete_out=np.zeros_like(native)
    for first in (-1,1):
        pa=(np.eye(2)+first*A)/2
        for last in (-1,1):
            pb=(np.eye(2)+last*As)/2;small=pb@pa;op=np.zeros((n,n),complex);op[np.ix_(indices,indices)]=small
            branch=op@native@op.conj().T;continued,ledger=continuation.continue_branch(branch)
            p=float(np.trace(branch).real);label=f'corner_{first}_{last}';arrays[label+'_poststate']=branch;arrays[label+'_continued']=continued
            effects.append(op.conj().T@op);complete_out+=continued;signed+=first*last*p
            if p>1e-12:
                try:modular.faithful_spectrum(branch[np.ix_(indices,indices)]/p)
                except ValueError:pure_rejections+=1
            apparatus.append({'first_A0':first,'second_As':last,'probability':p,'poststate_archive_key':label+'_poststate','continued_archive_key':label+'_continued','complete_native_continuation':ledger})
            require(abs(sum(r['joint_probability'] for r in ledger)-p)<1e-11,'lost sequential native outcome')
            number+=1
    qout,qledger=continuation.continue_branch(outside);arrays['complement_continued']=qout;effects.append(Q);complete_out+=qout
    apparatus.append({'status':'outside_corner','probability':float(np.trace(outside).real),'poststate_archive_key':'outside_unnormalized','continued_archive_key':'complement_continued','complete_native_continuation':qledger})
    require(abs(sum(r['joint_probability'] for r in qledger)-np.trace(outside).real)<1e-11,'lost complement continuation')
    effect_error=float(np.max(abs(sum(effects)-np.eye(n))));trace_error=float(abs(np.trace(complete_out)-1))
    real_correlation=float((np.trace(rho@A@As)).real);sequential_correlation=signed/prob
    # Compare to the supplied native H clock without identifying t and s.
    ids,e,v=model.spectral[3];u=(v*np.exp(-.37j*e))@v.conj().T;conditioned=inside[np.ix_(ids,ids)]/prob
    evolved=u@conditioned@u.conj().T;local=[list(ids).index(i) for i in indices]
    leakage=float(np.trace(evolved).real-np.trace(evolved[np.ix_(local,local)]).real)
    hc=h[np.ix_(indices,indices)];hc0=hc-np.trace(hc)*np.eye(2)/2;k0=K-np.trace(K)*np.eye(2)/2
    alpha=float(np.vdot(k0,hc0).real/np.vdot(k0,k0).real);fit=float(np.linalg.norm(hc0-alpha*k0))
    noninvariance=float(np.linalg.norm(h@P-P@h));full_eigen=np.linalg.eigvalsh(native);rank=int(np.count_nonzero(full_eigen>1e-10))
    full_rejected=False
    try:modular.faithful_spectrum(native)
    except ValueError:full_rejected=True
    checks={'same_interacting_H':True,'corner_numerically_faithful':bool(values[0]>1e-8),'complete_PQ_conditioning':bool(abs(np.trace(inside+outside)-1)<1e-11),'complete_sequential_observer_instrument':bool(effect_error<1e-11),'all_native_continuations_complete':bool(trace_error<1e-11),'GNS_representation':bool(max(gns_error,cyclic_error)<1e-10),'KMS_minus_i':bool(kms<1e-10),'real_sequential_correlation':bool(abs(sequential_correlation-real_correlation)<1e-11)}
    mutations={'full_M384_faithfulness_rejected':full_rejected,'pure_postmeasurement_M2_faithfulness_rejected':pure_rejections==4,'wrong_KMS_plus_i_rejected':wrong_kms>.1,'discarded_complement_rejected':bool(np.trace(outside).real>.9),'autonomous_native_corner_clock_rejected':leakage>.01 and fit>.01}
    require(all(checks.values()) and all(mutations.values()),'interacting modular observer failed')
    require(before=={str(p.relative_to(BASE)):digest(p) for p in paths},'source changed during experiment')
    result={'schema':'oph.c11.interacting_modular_observer.v1','construction_id':contract['construction_id'],'source':'actual C10 selected_protocol_poststate after the interacting global reference native .37/.29 protocol; no old one-load modular output transferred','reuse':{'modular_math':str(PROVIDER.relative_to(BASE))+' faithful_spectrum/modular_flow; finite standard GNS convention','continuation':'next/C10/spectral_readout.py GroupedEnergyInstrument.continue_branch, actual H and complete MIX with native time .37/probability .29'},'checks':checks,'mutations':mutations,'observer':{'algebra':'declared two-configuration P M384 P = M2; global occupation/record apparatus, not a local port or a partial trace','selection_ancestry':{'status':'Apparatus chosen after finite inspection of the known actual C10 output; not an a-priori prediction or uniquely derived observer','rule':'Lexicographically first basis pair with both diagonals above active_diagonal_threshold, conditional minimum eigenvalue above conditional_minimum_threshold, and raw offdiagonal magnitude above offdiagonal_threshold','active_diagonal_threshold':1e-8,'conditional_minimum_threshold':1e-4,'offdiagonal_threshold':1e-4,'execution':'The two disclosed configurations are fixed during execution; no fresh adaptive search or empirical physical target is used'},'configurations':[{'occupation_mask':m,'half_record_residues':list(f),'archive_index':i} for (m,f),i in zip(configurations,indices)],'corner_probability':prob,'complement_probability':float(np.trace(outside).real),'conditional_density':modular.matrix_data(rho),'conditional_eigenvalues':values.tolist(),'dimensionless_modular_frequency':frequency,'full_source_numerical_rank':rank,'full_source_rank_threshold':1e-10,'apparatus_eigenframe_phase':'Each ordered density eigenvector has its largest-magnitude component positive real','full_source_faithful_on_M384':False,'GNS_reference':'vec(sqrt(rho)) on the four-dimensional Hilbert-Schmidt representation; no physical reference copy or subjective clock','no_floor_or_pseudoinverse':True},'modular_observable':{'s':s,'observable_A0':modular.matrix_data(A),'observable_As':modular.matrix_data(As),'algebraic_complex_correlation':{'real':float(np.trace(rho@A@As).real),'imag':float(np.trace(rho@A@As).imag)},'executed_real_two_PVM_correlation':sequential_correlation,'scope':'Only the real correlation is realized by the declared sequential PVM; the complex algebraic correlation is not a classical joint probability.'},'complete_observer_records':apparatus,'native_clock_comparison':{'supplied_native_duration':.37,'corner_leakage_under_H':leakage,'corner_H_commutator_norm':noninvariance,'best_compressed_H_equals_alpha_K_plus_scalar_alpha':alpha,'compressed_generator_fit_residual':fit,'identification':'No autonomous corner H flow, physical-time conversion or subjective/experienced-time identification.'},'errors':{'GNS':gns_error,'GNS_cyclic':cyclic_error,'KMS':kms,'wrong_KMS_shift':wrong_kms,'conditional_density_hermiticity':float(np.max(abs(rho-rho.conj().T))),'observer_effect_completeness':effect_error,'continued_state_trace':trace_error},'resources':['actual canonical interacting state and ideal preparation inherited from C10','explicit fixed two-configuration coherent corner apparatus, chosen using the known simulated input density; tomography and sample costs not derived','frozen premeasurement state-dependent observable family and injected modular parameter','full complement and all subnormalized poststates retained','additional sequential projective apparatus, not synthesized by native C/E','native continuation schedule supplied; no autonomous clock or consciousness map'], 'scope':'Finite source-conditioned modular algebra and complete continued apparatus. Full source and pure observed branches are not faithful M2/M384 states. Modular flow is not inherited as the postmeasurement state flow. Numerical matrix propagation has no certified enclosure. No physical species, vacuum, mass, gravity or official row award.','source_sha256':before,'full_criterion_met':False,'official_rows_changed':[]}
    return result,arrays

def compare_scientific(a,b,key=''):
    if isinstance(a,dict):return isinstance(b,dict) and a.keys()==b.keys() and all(compare_scientific(a[k],b[k],k) for k in a)
    if isinstance(a,list):return isinstance(b,list) and len(a)==len(b) and all(compare_scientific(x,y,key) for x,y in zip(a,b))
    if type(a) is float:
        if key in {'s','supplied_native_duration','full_source_rank_threshold','active_diagonal_threshold','conditional_minimum_threshold','offdiagonal_threshold'}:return type(b) is float and a==b
        return type(b) in (float,int) and bool(np.isclose(a,b,atol=1e-10,rtol=1e-10))
    return type(a) is type(b) and a==b

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True);mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    result,arrays=execute();path=HERE/'INTERACTING_MODULAR_OBSERVER.json';archive=HERE/'interacting_modular_observer_states.npz'
    if args.write:
        np.savez_compressed(archive,**arrays);result['state_archive']={'filename':archive.name,'sha256':digest(archive)};path.write_text(json.dumps(result,indent=2)+'\n')
    else:
        saved=json.loads(path.read_text());require(digest(archive)==saved['state_archive']['sha256'],'state archive custody')
        saved.pop('state_archive');require(compare_scientific(result,saved),'scientific receipt changed')
        with np.load(archive,allow_pickle=False) as z:require(z.files==list(arrays) and all(np.allclose(value,z[key],atol=1e-10,rtol=1e-10) for key,value in arrays.items()),'continued state changed')
    print(json.dumps({'checks':result['checks'],'mutations':result['mutations'],'frequency':result['observer']['dimensionless_modular_frequency'],'corner_probability':result['observer']['corner_probability']},indent=2))
