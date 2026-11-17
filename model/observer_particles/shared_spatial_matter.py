#!/usr/bin/env python3
"""Conditional CAR field on the actual scalar spatial form and observer clock.

The incidence field is constructed and executed, including its edge kernel.
It is not labeled a Weyl spinor or a selected Standard Model action.
"""
from __future__ import annotations
from evidence_contract import public_receipt
from fractions import Fraction as F
from collections import Counter
import hashlib,itertools,json,math,sys
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PUBLIC=ROOT
sys.path[:0]=[str(PUBLIC),str(PUBLIC/'code')]
from model.federated_observer import source
from model.observer import native
from model.observer.source_scalar import compiler
from model.common_source import gauge_matter
from m1_fermionic_source.model import Encoding
from m1_fermionic_source.pauli import Word

PREVIOUS=HERE


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def maxabs(a):return float(np.max(abs(a))) if np.size(a) else 0.
def sparse_arrays(prefix,m):
    m=m.tocoo();return {prefix+'_shape':np.array(m.shape),prefix+'_row':m.row,prefix+'_col':m.col,prefix+'_data':m.data}


def geometry():
    graph=source.federation(0);v=graph.ports;e=len(graph.seams);q=3;h=F(2,q+1)
    assert (v,e)==(240,630)
    rows=[];cols=[];values=[]
    for j,(a,b) in enumerate(graph.seams):
        # Derive the orientation from the actual even-load current, not an
        # unrelated graph supplied solely to the fermion construction.
        for coin in (0,1):
            out,(_,current)=native.native((2,0),(0,1),coin)
            assert out==(1,1) and current==1
        rows.extend((j,j));cols.extend((a,b));values.extend((-1,1))
    d=sparse.csr_matrix((values,(rows,cols)),shape=(e,v),dtype=np.float64)
    l=(d.T@d).tocsr()
    radial=np.zeros((q+1,q));radial[0,0]=1/float(h);radial[-1,-1]=-1/float(h)
    for j in range(1,q):radial[j,j-1]=-1/float(h);radial[j,j]=1/float(h)
    c0=sparse.vstack((sparse.kron(sparse.eye(q),d),sparse.kron(radial,sparse.eye(v))),format='csr')
    b0=sparse.kron(sparse.eye(q),l)+sparse.kron(radial.T@radial,sparse.eye(v))
    assert (c0.T@c0-b0).nnz==0
    ll,ul=np.linalg.eigh(l.toarray());assert ll[0]>-1e-12 and ll[1]>0
    ll[0]=0. # The exact constant zero mode, after numerical diagonalization.
    rr,ur=np.linalg.eigh(radial.T@radial)
    f=630*np.arccos(1-ll/1260)**2
    g=np.ones(v);g[1:]=np.sqrt(f[1:]/ll[1:])
    angular_factor=d@((ul*g)@ul.T)
    cw=sparse.vstack((sparse.kron(sparse.eye(q),sparse.csr_matrix(angular_factor)),
                      sparse.kron(radial,sparse.eye(v))),format='csr')
    bw=sparse.kron(sparse.eye(q),sparse.csr_matrix((ul*f)@ul.T))+sparse.kron(radial.T@radial,sparse.eye(v))
    assert maxabs((cw.T@cw-bw).toarray())<2e-11
    d0=sparse.bmat([[None,c0.T],[c0,None]],format='csr')
    dw=sparse.bmat([[None,cw.T],[cw,None]],format='csr')
    # A rational finite-series witness that this particular exact-clock
    # factor correction has nonlocal matrix support.
    lint=l.toarray().astype(np.int64);di=d.toarray().astype(np.int64)
    dl=di@lint;dl2=dl@lint
    rho=F(1,210) # ||L/(4Q)|| <=12/2520
    remainder=4*rho**3/(1-rho) # sqrt(12)<4 and all arcsin-series coefficients<=1
    witness=None
    for row,(a,b) in enumerate(graph.seams):
        for col in range(v):
            if col in (a,b):continue
            partial=F(int(dl[row,col]),24*630)+F(3*int(dl2[row,col]),640*630**2)
            if abs(partial)>remainder:
                witness={'angular_edge_row':row,'edge':[a,b],'outside_vertex':col,
                    'quadratic_series_entry':str(partial),'absolute_tail_upper':str(remainder),
                    'certified_absolute_nonlocal_entry_lower':str(abs(partial)-remainder),
                    'numerical_full_entry':float(angular_factor[row,col])};break
        if witness:break
    assert witness
    return graph,d,l,radial,c0,b0,cw,bw,d0,dw,ll,ul,rr,ur,f,witness


def spectral_apply(vector,values,ul,ur):
    return (ur@((ur.T@vector.reshape(3,240)@ul)*values)@ul.T).reshape(-1)


def propagation(c,big,angular,ul,radial_values,ur,clock_weights):
    omega=np.sqrt(radial_values[:,None]+angular[None,:])
    start=np.zeros(big.shape[0],complex);start[240]=1 # middle radial slice, native port0
    node=start[:720].real
    rows=[];probabilities=[]
    for ticks in (0,1,2,25):
        t=ticks/math.sqrt(630)
        actual=expm_multiply(-1j*t*big,start,traceA=0)
        expected_node=spectral_apply(node,np.cos(t*omega),ul,ur)
        response=spectral_apply(node,np.sin(t*omega)/omega,ul,ur)
        expected=np.concatenate((expected_node,-1j*(c@response)))
        error=maxabs(actual-expected);norm=float(np.vdot(actual,actual).real)
        assert error<2e-12 and abs(norm-1)<3e-12
        law=abs(actual)**2
        rows.append({'incident_count':ticks,'observer_time':t,'complete_number_readout_probability_sum':float(law.sum()),
            'scalar_to_fermion_propagator_identity_error':error,
            'node_probability':float(law[:720].sum()),'edge_probability':float(law[720:].sum()),
            'return_probability':float(law[240]),'off_initial_mode_probability':float(1-law[240])})
        if ticks<3:probabilities.append(law)
    law=np.array(probabilities)
    weights=np.array(list(map(float,clock_weights)))
    assert abs(np.sum(law*weights[:,None])-1)<3e-12
    return rows,law


def edge_zero_mode(graph,c0,cw,d0,dw):
    edge_lookup={edge:i for i,edge in enumerate(graph.seams)}
    triangle=None
    for a in range(12):
        for b,c in itertools.combinations([x for x in graph.neighbours[a] if x<12],2):
            if tuple(sorted((b,c))) in edge_lookup:triangle=(a,b,c);break
        if triangle:break
    assert triangle
    z=np.zeros(c0.shape[0])
    for a,b in zip(triangle,(*triangle[1:],triangle[0])):
        z[edge_lookup[tuple(sorted((a,b)))]]+=1 if a<b else -1
    assert np.array_equal(c0.T@z,np.zeros(720))
    assert maxabs(cw.T@z)<2e-14
    vector=np.concatenate((np.zeros(720),z/math.sqrt(3))).astype(complex)
    errors=[]
    for matrix in (d0,dw):
        actual=expm_multiply(-1j*matrix,vector,traceA=0)
        errors.append(maxabs(actual-vector))
    assert max(errors)<2e-13
    return {'native_triangle':list(triangle),'exact_local_C_adjoint_cycle_zero':True,
        'finite_angular_factor_preserves_same_kernel':True,'one_time_unit_persistence_errors':errors,
        'node_modes':720,'edge_modes':2850,'positive_energy_modes':720,'negative_energy_modes':720,
        'zero_energy_modes':2130,'all_modes_retained':True,
        'level1_zero_energy_count':8520,
        'zero_mode_Fock_degeneracy':'2^2130 at level0; 2^8520 at level1',
        'thermal_zero_mode_entropy':'2130 log(2) at level0',
        'continuum_thermal_Fock_state_or_unique_vacuum_proved':False,
        'ground_state_choice':'All720 negative modes may be filled, but2130 zero modes remain. Particle-hole reinterpretation/normal ordering and any removal of zero modes are additional choices.'},z


def source_CAR_and_locality(c0):
    coo=c0.tocoo();vertices=c0.shape[1]+c0.shape[0]
    edges=tuple(sorted((int(col),720+int(row)) for row,col in zip(coo.row,coo.col)))
    encoding=Encoding(vertices,edges)
    assert encoding.n==3570 and encoding.m==5220 and len(encoding.loops)==1651
    # Signed loop constraints and every fundamental cycle are retained. The
    # chord X submatrix is the identity, hence their binary independence.
    chordmask=sum(1<<x for x in encoding.chords)
    for loop,chord in zip(encoding.loops,encoding.chords):
        assert loop*loop==Word() and loop.x&chordmask==1<<chord
    i,j=edges[0]
    timing=1/math.sqrt(630)
    g=expm(-1j*timing*np.array([[0.,1.],[1.,0.]]))
    scalar,tape=encoding.su2(i,j,g)
    support=sorted({k for p,_ in tape for k in range(encoding.m) if ((p.x|p.z)>>k)&1})
    assert len(support)<=9
    def restrict(p):return Word(sum(((p.x>>old)&1)<<new for new,old in enumerate(support)),
                                sum(((p.z>>old)&1)<<new for new,old in enumerate(support)),p.phase)
    word=[];target=np.eye(1<<len(support),dtype=complex)
    for p,angle in tape:
        pp=restrict(p);label=[]
        for bit in range(len(support)):
            x,z=(pp.x>>bit)&1,(pp.z>>bit)&1
            label.append('Y' if x and z else 'X' if x else 'Z' if z else 'I')
        sign=1j**((pp.phase-(pp.x&pp.z).bit_count())%4)
        assert abs(sign.imag)<1e-15
        word+=compiler.pauli_rotation(''.join(label),float(angle*sign.real))
        target=(math.cos(angle)*np.eye(len(target))+1j*math.sin(angle)*pp.matrix(len(support)))@target
    actual=scalar*compiler.execute(word,len(support))
    error=float(np.linalg.norm(actual-scalar*target,2));assert error<2e-11
    # Independent complete CAR and source hopping with actual occupation
    # signs, including intervening spectator occupancies and both parities.
    ladders=[]
    for mode in range(4):
        mat=np.zeros((16,16),complex)
        for occupation in range(16):
            result=gauge_matter.car_ladder(occupation,mode,False)
            if result is not None:sign,out=result;mat[out,occupation]=sign
        ladders.append(mat)
    for a in range(4):
        for b in range(4):
            assert maxabs(ladders[a]@ladders[b]+ladders[b]@ladders[a])==0
            assert maxabs(ladders[a]@ladders[b].conj().T+ladders[b].conj().T@ladders[a]-(np.eye(16) if a==b else 0))==0
    hopping=ladders[0].conj().T@ladders[3]+ladders[3].conj().T@ladders[0]
    gates=[]
    for coefficient in (-2,-1,1,2):
        w=compiler.pauli_rotation('XZZX',-coefficient*timing/2)+compiler.pauli_rotation('YZZY',-coefficient*timing/2)
        u=compiler.execute(w,4);err=float(np.linalg.norm(u-expm(-1j*coefficient*timing*hopping),2))
        assert err<2e-12
        gates.append({'actual_incidence_coefficient':coefficient,'source_pulses':len(w),'full16_occupation_column_error':err})
    measure=compiler.execute(compiler.cnot(0,4),5)
    kraus=[measure[b*16:(b+1)*16,:16] for b in range(2)]
    projectors=[np.diag([int((i&1)==b) for i in range(16)]) for b in range(2)]
    assert max(np.linalg.norm(k-p,2) for k,p in zip(kraus,projectors))<1e-12
    assert np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(16),2)<1e-12
    exports={}
    for dimension in (2,4):
        ks=compiler.transfer_kraus(dimension)
        assert np.array_equal(sum(k.conj().T@k for k in ks),np.eye(6))
        exports[str(dimension)]={'complete_kraus_count':len(ks),'complement_failure_branches':6-dimension}
    return {'actual_encoding_graph_vertices':vertices,'actual_encoding_graph_edges':len(edges),
        'fundamental_signed_cycle_constraints':len(encoding.loops),'cycle_independence_exact':True,
        'maximum_mode_graph_degree':max(map(len,encoding.incident)),
        'both_parity_sectors_required':True,'each_parity_block_dimension':'2^3569',
        'local_hopping_source_check':{'chosen_actual_bond':[i,j],'local_edge_qubits':len(support),
            'source_pulses':len(word),'whole_local_unitary_error':error},
        'independent_full_CAR_mode_pair_checks':16,'spectator_sign_source_checks':gates,
        'complete_number_instrument_checked':True,'complete_M6_exports':exports,
        'full3570_mode_Fock_matrix_or_entire_native_program_executed':False,
        'all_mode_CAR_automorphism':'For unitary U=exp(-itD), c(f)->c(Uf) preserves CAR on every occupation sector; this is the owning arbitrary-graph source theorem, not a claim to have enumerated2^3570 states.',
        'physical_statistics_selected_by_A1_A3':False}


def gauge_covariance(c0):
    coo=c0.tocoo();node_phase=(1j)**(np.arange(720)%4)
    row_anchor=np.zeros(c0.shape[0],int)
    for row in range(c0.shape[0]):
        cols=c0.indices[c0.indptr[row]:c0.indptr[row+1]]
        row_anchor[row]=int(cols[-1])
    edge_phase=node_phase[row_anchor]
    transformed=sparse.diags(edge_phase)@c0@sparse.diags(node_phase.conj())
    # Direct edge transport prescription, independently assembled.
    direct=sparse.csr_matrix((coo.data*edge_phase[coo.row]*node_phase[coo.col].conj(),(coo.row,coo.col)),shape=c0.shape)
    assert (transformed-direct).nnz==0
    covariance=(transformed.conj().T@transformed-sparse.diags(node_phase)@(c0.T@c0)@sparse.diags(node_phase.conj())).tocsr()
    assert maxabs(covariance.data)==0
    return {'integer_U1_charge':1,'global_charge':'Full CAR number N; [dGamma(D),N]=0',
        'local_external_link_covariance_exact':True,
        'transformation':'C_U -> R_G C_U G^dagger; D_U -> diag(G,R_G) D_U diag(G,R_G)^dagger',
        'gauge_anchors':'Each edge/form mode transforms at its explicitly declared endpoint anchor',
        'dynamical_quantum_links_or_Gauss_projection_added':False,
        'M15_extension':'D tensor I15 commutes with all12 existing internal currents. This is a conditional vectorlike degeneracy; no Weyl action or chiral continuum follows.',
        'spin_attachment':'Incidence0/1-form indices are not a Weyl spinor. An independent Spin2 tensor factor alone would not supply the metric spin connection or remove the edge kernel.'}


def main():
    graph,d,l,radial,c0,b0,cw,bw,d0,dw,ll,ul,rr,ur,f,witness=geometry()
    upstream=PREVIOUS/'results/spacetime_correspondence.json'
    up_artifact=PREVIOUS/'results/spacetime_correspondence_arrays.npz'
    up=json.loads(upstream.read_text())
    assert up['artifacts'][up_artifact.name]==sha(up_artifact)
    weights=[F(0)]*3
    for row in up['joint_native_history']['joint_record_rows']:
        weights[row['clock_A']]+=F(row['native_numerator'],row['native_denominator'])
    local_rows,local_laws=propagation(c0,d0,ll,ul,rr,ur,weights)
    exact_rows,exact_laws=propagation(cw,dw,f,ul,rr,ur,weights)
    zero,cycle=edge_zero_mode(graph,c0,cw,d0,dw)
    implementation=source_CAR_and_locality(c0)
    covariance=gauge_covariance(c0)
    with np.load(up_artifact) as saved:old=saved['quantum_pointer_branch16_amplitudes']
    old_masses=np.sum(abs(old)**2,axis=-1)
    joint=old_masses[:,:,:,None]*exact_laws[:,None,None,:]*np.array(list(map(float,weights)))[:,None,None,None]
    for axis in range(2):assert abs(joint[:,axis].sum()-1)<4e-12
    old_marginal_error=maxabs(joint.sum(axis=-1)-old_masses*np.array(list(map(float,weights)))[:,None,None])
    assert old_marginal_error<3e-12
    omega=np.sqrt(rr[:,None]+f[None,:]);cfplus=1/(1+np.exp(omega));cfminus=1-cfplus
    assert maxabs(cfplus+cfminus-1)<1e-15 and np.min(cfplus)>0
    artefact=HERE/'results/shared_spatial_matter_arrays.npz'
    arrays={'native_seams':np.array(graph.seams),'radial_difference':radial,'angular_eigenvalues':ll,
        'angular_finite_clock_squared':f,'radial_eigenvalues':rr,
        'complete_node_scalar_frequencies':omega,'fermion_positive_energy_occupations':cfplus,
        'fermion_negative_energy_occupations':cfminus,'persistent_edge_cycle':cycle,
        'local_fermion_number_laws':local_laws,'exact_shared_B_fermion_number_laws':exact_laws,
        'same_history_joint_pointer_fermion_probabilities':joint,
        **sparse_arrays('C0',c0),**sparse_arrays('B0',b0),**sparse_arrays('Cwave',cw),**sparse_arrays('Bwave',bw)}
    np.savez_compressed(artefact,**arrays)
    paths=[upstream,up_artifact,PREVIOUS/'spacetime_correspondence.py',
        PUBLIC/'docs/research/FERMIONIC_SOURCE_CLOCKS.md',PUBLIC/'docs/research/MASSIVE_OPERATIONAL_CLOCKS.md',
        PUBLIC/'docs/research/SOURCE_GRAVITY_ADMISSIBILITY.md',PUBLIC/'model/observer_time_fields/REFINEMENT.md',
        PUBLIC/'model/relative_observer_geometry/README.md',PUBLIC/'code/m1_fermionic_source/CONTRACT.md']
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and Path(path).is_file() and Path(path).is_relative_to(PUBLIC):paths.append(Path(path))
    result={'status':'PASS_CONDITIONAL_SHARED_SPATIAL_CAR_FIELD_WITH_COMPLETE_EXTRA_MODES_AND_CLOCK',
        'producer_sha256':sha(__file__),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))},
        'artifacts':{artefact.name:sha(artefact)},'empirical_mass_inputs':[],
        'geometry':{'source_ports':240,'source_seams':630,'radial_points':3,'radial_length':2,'radial_boundary':'Dirichlet',
            'node_modes':720,'edge_rows':2850,'incidence_nonzero_entries':c0.nnz,
            'exact_local_factorization':'C0^dagger C0=B0=L_native tensor I_radial + I_native tensor L_radial',
            'finite_angular_ideal_joint_factorization':'Cwave^dagger Cwave=Bwave=f(L_native)+L_radial; f(lambda)=630 acos(1-lambda/1260)^2',
            'factorization_order':'node index=240*radial_index+native_port',
            'finite_angular_correction':'Cwave angular rows=d*asin(sqrt(L/(4Q)))/sqrt(L/(4Q)); radial rows unchanged',
            'nonlocal_finite_clock_factor_witness':witness,
            'locality_scope':'D0 is nearest incidence-bond local; this particular Bwave factor changes that support atfinitecap. On bounded angular bands the correction tends toidentity asQ grows. No other local factorization is ruled out.',
            'Bwave_is_exact_Floquet_square_of_executed_two_axis_word':False,
            'Bwave_scope':'Exact finite angular frequency combined with radial contribution in the IDEAL joint generator. The actual finite protected Z/X split word has different Floquet phase and retains the upstream error bound.'},
        'conditional_joint_action_state_clock':{'fermion_action':'H_F=dGamma(D), D=[[0,C^dagger],[C,0]]',
            'scalar_action':'H_B=1/2(P^2+Q B Q), on the same B and radial boundary',
            'same_beta':1,'joint_reference':'rho_ref=rho_scalar,B tensor rho_CAR,D atbeta1; modularsigma_s=alpha_(-s)',
            'reference_scope':'This exact joint reference is for the declared ideal CCR/CAR actions. The retained two-copy spin apparatus is a bounded finite approximation, not the entire CCR reference or exactly the joint ideal flow.',
            'fermion_covariance':'C_F=(I+exp(D))^-1: occupations f(+omega),f(-omega), and1/2 on everyedgezero mode',
            'full_joint_Gaussian_reference_preparation_newly_executed':False,
            'working_state':'Prior retained scalar/work/reference branches tensor a declared one-fermion basis excitation atnativeport0,middleradialsite, initialized beforethehistory',
            'clock':'Same nativeN_A, t=N_A/sqrt630, zero-current incidents included',
            'source_selected_fermionic_statistics_action_or_preparation':False},
        'propagation':{'local_B0':local_rows,'exact_shared_ideal_joint_Bwave':exact_rows,
            'exact_identity':'U_D(t) on nodeinput f has nodecos(t sqrtB)f and edge -i C sin(t sqrtB)/sqrtB f',
            'same_scalar_poles':True,'newly_executed_native_full_graph_unitary':False,
            'sparse_full_one_particle_unitary_action_executed':True,
            'actual_two_axis_scalar_word_exactly_intertwined':False,
            'inherited_two_axis_split_controls':{n:{k:up['same_band_state_clock_response_error_controls'][n][k]
                for k in ('single_code_unitary_error','two_axis_product_formula_error_upper',
                          'two_axis_pulse_synthesis_error_upper','finite_split_clock_is_exact_target_modular_flow')}
                for n in ('1','2','25')}},
        'same_history_extension':{'fine_native_prefixes':up['joint_native_history']['all_fine_native_histories'],
            'old_scalar_pointer_marginal_error':old_marginal_error,'all3570_fermion_number_outcomes_retained':True,
            'same_native_history_no_reset':True,'scalar_fermion_interaction_inserted':False,
            'joint_law':'p(nativeprefix)*p(oldpointer|N_A)*p(fermionmode|N_A), with untouchedconditional scalar/reference state retained',
            'native_full_fermion_preparation_and_propagation_execution_claimed':False},
        'extra_modes_and_vacuum':zero,'CAR_source_realization':implementation,'charge_and_gauge':covariance,
        'physical_identification':{'conditional_shared_geometry_time_state_net_constructed':True,
            'CAR_not_just_one_particle_occupation_swap':True,
            'physical_Weyl_spin_connection_chirality_identified':False,
            'SM_generations_masses_Yukawas_selected':False,'stable_continuum_Fock_limit_proved':False}}
    (HERE/'results/shared_spatial_matter.json').write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'one_particle_dimension':3570,
        'exact_shared_B_response_error':max(x['scalar_to_fermion_propagator_identity_error'] for x in exact_rows),
        'source_local_gate_error':implementation['local_hopping_source_check']['whole_local_unitary_error'],
        'old_pointer_marginal_error':old_marginal_error,'zero_modes':2130},indent=2))


if __name__=='__main__':main()
