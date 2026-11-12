"""Volume-linear coherent error with every CAR factor and reference retained."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.linalg import expm
from model.lepton_mass_compiler import analyzer,compiler as native
from model.openwave_kinetic_phase_cycle.compiler import formal_order_certificate,rounds_for_error as old_rounds,trace_bound as old_bound
from model.openwave_localized_rest_packet.packet import certify,sqrt_bounds


def block_radius(spacing=F(1)):
    if type(spacing)not in(int,F)or spacing<=0:raise ValueError('positive rational spacing required')
    # |sin(p)| <= sqrt(3)/a < 2/a; centered four-mode kinetic norm2|kappa|.
    return F(7,3)*(4/F(spacing)+1)


def trace_bound(rounds,time=F(1),*,sites=1,spacing=F(1)):
    if type(rounds)is not int or rounds<1 or type(sites)is not int or sites<1:
        raise ValueError('positive integer round and site counts required')
    if type(time)not in(int,F)or time<0:raise ValueError('nonnegative rational time required')
    r=block_radius(spacing)
    if rounds<r*time:raise ValueError('block Taylor radius exceeds one')
    return min(F(2),F(6*sites,10)*r**5*F(time)**5/rounds**4)


def rounds_for_error(error,time=F(1),*,sites=1,spacing=F(1)):
    error=F(error);time=F(time)
    if error<=0 or time<0:raise ValueError('positive error and nonnegative time required')
    rt=block_radius(spacing)*time;lo=max(1,(rt.numerator+rt.denominator-1)//rt.denominator);hi=lo
    while trace_bound(hi,time,sites=sites,spacing=spacing)>error:hi*=2
    while lo<hi:
        mid=(lo+hi)//2
        if trace_bound(mid,time,sites=sites,spacing=spacing)<=error:hi=mid
        else:lo=mid+1
    return lo


def block(effect,kappa):
    if not 0<effect<1:raise ValueError('interior source effect required')
    c=[native._annihilator(4,i)for i in range(4)];w=np.array([[effect,np.sqrt(effect*(1-effect))],[np.sqrt(effect*(1-effect)),1-effect]])
    V=sum((w[i,j]*c[i].conj().T@c[2+j].conj().T for i in range(2)for j in range(2)),np.zeros((16,16),complex));V+=V.conj().T
    number=sum((a.conj().T@a for a in c),np.zeros((16,16),complex));H=kappa*number+V
    counts=np.array([((mask>>1)&1)+((mask>>3)&1)for mask in range(16)])
    phases=[np.diag((1j)**(r*counts))for r in range(4)]
    letters=[R.conj().T@H@R for R in phases]
    assert np.linalg.norm(V,2)<1+1e-13
    assert np.linalg.norm(kappa*(number-2*np.eye(16)),2)<=2*abs(kappa)+1e-13
    return letters,sum(letters)/4


def step(letters,duration):
    p=1/(4-4**(1/3));order=[(0,1/8),(1,1/8),(2,1/8),(3,1/4),(2,1/8),(1,1/8),(0,1/8)]
    def second(t):
        U=np.eye(16,dtype=complex)
        for r,w in order:U=expm(-1j*t*w*letters[r])@U
        return U
    A=second(p*duration);B=second((1-4*p)*duration)
    return A@A@B@A@A


def complete_control():
    effects=[float(sum(map(F,x),F(0))/2)for x in analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']]
    pairs=[(effects[1],.7),(effects[2],-1.2)];t=.8;n=16;actual=[];ideal=[];errors=[];marginals=[]
    J=np.eye(16)[:,[0,1,4,5]]
    for e,k in pairs:
        letters,mean=block(e,k);U=np.linalg.matrix_power(step(letters,t/n),n);V=expm(-1j*t*mean)
        actual.append(U);ideal.append(V);errors.append(float(np.linalg.norm(U-V,2)))
        out=U@J;probs=[float(np.linalg.norm(out[[word|mask for word in(0,1,4,5)]],'fro')**2/4)for mask in(0,2,8,10)]
        assert abs(sum(probs)-1)<2e-13;marginals.append(probs)
    full=np.kron(actual[0],actual[1]);target=np.kron(ideal[0],ideal[1]);full_error=float(np.linalg.norm(full-target,2))
    assert full_error<=sum(errors)+2e-13
    # Complete maximally entangled physical/reference input, all16 ancillary branches.
    J2=np.kron(J,J);A=full@J2;B=target@J2;overlap=np.vdot(A.ravel(),B.ravel())/16
    trace=2*np.sqrt(max(0.,1-abs(overlap)**2));bound=float(trace_bound(n,F(4,5),sites=1))
    assert trace<=bound+1e-12
    probabilities=np.outer(marginals[0],marginals[1]);assert abs(probabilities.sum()-1)<4e-13
    # Invalid noncommuting inter-block coupling prevents exact pulse factorization.
    Z=np.diag([(-1)**(mask&1)for mask in range(16)]);X=native._annihilator(4,0);X+=X.conj().T
    coupled=np.kron(block(*pairs[0])[0][0],np.eye(16))+np.kron(np.eye(16),block(*pairs[1])[0][0])+.2*np.kron(Z,X)
    defect=float(np.linalg.norm(expm(-.01j*coupled)-np.kron(expm(-.01j*block(*pairs[0])[0][0]),expm(-.01j*block(*pairs[1])[0][0])),2))
    assert defect>1e-4
    return {'full_two_block_Fock_rows':256,'physical_input_and_reference_dimension':16,
        'unhalved_reference_trace_distance_diagnostic':float(trace),'complete_trace_upper':bound,
        'full_operator_error':full_error,'sum_of_block_operator_errors':sum(errors),
        'all_original_ancillary_records':probabilities.tolist(),
        'coupled_block_mutation_factorization_defect':defect,
        'numeric_controls_are_not_interval_enclosures':True}


def generic_momentum_factorization():
    sigma=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.]))
    eps=np.array([[0,1],[-1,0]],complex);effect=.27;w=np.array([[effect,np.sqrt(effect*(1-effect))],[np.sqrt(effect*(1-effect)),1-effect]])
    def idx(p,leg,spin,flag):return 8*p+4*leg+2*spin+flag
    groups=[[idx(p,0,spin,f)for f in range(2)]+[idx(1-p,1,1-spin,f)for f in range(2)]for p in range(2)for spin in range(2)]
    output=[]
    for vector in((.4,-1/3,1/7),(0,0,1.),(0,0,0)):
        S=sum((v*a for v,a in zip(vector,sigma)),np.zeros((2,2),complex));length=float(np.linalg.norm(vector))
        if length:
            _,U=np.linalg.eigh(S);U=U[:,::-1];U[:,1]/=np.linalg.det(U)
        else:U=np.eye(2,dtype=complex)
        assert np.linalg.norm(U.conj().T@eps@U.conj()-eps)<1e-13
        h=np.zeros((16,16),complex);pair=np.zeros_like(h);rotation=np.zeros_like(h)
        for p in range(2):
            for leg in range(2):
                for flag in range(2):
                    indices=[idx(p,leg,t,flag)for t in range(2)]
                    h[np.ix_(indices,indices)]=(-1 if p==0 else 1)*S
                    rotation[np.ix_(indices,indices)]=U
            for a in range(2):
                for b in range(2):
                    for f in range(2):
                        for g in range(2):
                            i,j=idx(p,0,a,f),idx(1-p,1,b,g);pair[i,j]=eps[a,b]*w[f,g];pair[j,i]=-pair[i,j]
        transformed_h=rotation.conj().T@h@rotation;transformed_pair=rotation.conj().T@pair@rotation.conj()
        H=np.block([[transformed_h,transformed_pair],[transformed_pair.conj().T,-transformed_h.T]])
        errors=[]
        for number,group in enumerate(groups):
            outside=[j for j in range(16)if j not in group]
            errors.extend([np.linalg.norm(transformed_h[np.ix_(group,outside)]),np.linalg.norm(transformed_pair[np.ix_(group,outside)])])
            p,spin=divmod(number,2);kappa=(-1 if p==0 else 1)*(1 if spin==0 else -1)*length
            errors.append(np.linalg.norm(transformed_h[np.ix_(group,group)]-kappa*np.eye(4)))
            block_pair=transformed_pair[np.ix_(group,group)];expected=np.block([[np.zeros((2,2)),(-1)**spin*w],[-(-1)**spin*w.T,np.zeros((2,2))]])
            errors.append(np.linalg.norm(block_pair-expected))
        assert max(errors)<2e-13
        output.append({'momentum':list(vector),'full_CAR_modes':16,'complete_Nambu_rows':32,'four_complete_CAR_factors':groups,'maximum_intertwiner_error':float(max(errors))})
    return output


def state_energy_error(delta,energy_range,target_mean_upper):
    delta=F(delta);energy_range=F(energy_range);target_mean_upper=F(target_mean_upper)
    if min(delta,energy_range,target_mean_upper)<0:raise ValueError('nonnegative norm and energy bounds required')
    root=sqrt_bounds(energy_range*target_mean_upper)[1]
    return 2*root*delta+energy_range*delta**2


def state_sensitive_rounds(budget,time,sites,energy_range,target_mean_upper):
    budget=F(budget);time=F(time)
    if budget<=0:raise ValueError('positive physical energy budget required')
    rt=block_radius()*time;lo=max(1,(rt.numerator+rt.denominator-1)//rt.denominator);hi=lo
    def acceptable(n):
        # In the admissible branch epsilon<2, the unitary bound is epsilon/2.
        eps=trace_bound(n,time,sites=sites)
        return eps<=F(1,10**12) and state_energy_error(eps/2,energy_range,target_mean_upper)<=budget
    while not acceptable(hi):hi*=2
    while lo<hi:
        mid=(lo+hi)//2
        if acceptable(mid):hi=mid
        else:lo=mid+1
    return lo


def resources():
    rows=[]
    for side in(3,4,16):
        s=side**3;eps=F(1,10**12);n=rounds_for_error(eps,sites=s);old=old_rounds(eps,sites=s)
        rows.append({'side':side,'sites':s,'target_trace_error':str(eps),'cycles':n,'prior_global_norm_cycles':old,
            'cycle_reduction_factor':float(F(old,n)),'complete_trace_bound':str(trace_bound(n,sites=s))})
    packet=certify();s=packet['spatial_sites'];time=F(packet['source_time_horizon']);d=F(packet['mean_source_energy_excess_upper']);energy=F(282*s)
    eps=min(F(1,10**12),d/(10**6*energy));n=rounds_for_error(eps,time,sites=s);bound=trace_bound(n,time,sites=s)
    assert energy*bound<=d/10**6
    old=packet['finite_full_Fock_compiler']['cycles'];assert old//n>10**28
    mean=F(packet['electron_mass_interval'][1])+d
    refined=state_sensitive_rounds(d/10**6,time,s,energy,mean)
    refined_trace=trace_bound(refined,time,sites=s)
    refined_energy=state_energy_error(refined_trace/2,energy,mean)
    assert refined_energy<=d/10**6 and refined<n
    return {'finite_examples':rows,'same_localized_packet':{'spatial_sites':s,'source_time':str(time),
        'cycles':n,'prior_global_norm_cycles':old,'cycle_reduction_factor':float(F(old,n)),
        'complete_trace_bound':str(bound),'physical_energy_mean_error_upper':str(energy*bound),
        'physical_energy_error_budget':str(d/10**6),'all_ancilla_outputs_retained':True,
        'same_target_state_and_coefficients':True,'full_large_carrier_executed':False},
        'state_sensitive_same_packet':{'cycles':refined,'prior_global_norm_cycles':old,
            'cycle_reduction_factor':float(F(old,refined)),'further_reduction_from_state_energy_bound':float(F(n,refined)),
            'complete_trace_bound':str(refined_trace),'physical_energy_mean_error_upper':str(refined_energy),
            'physical_energy_error_budget':str(d/10**6),'target_mean_energy_upper':str(mean),
            'same_physical_record_and_energy_guarantees':True,'all_ancilla_outputs_retained':True,
            'bound_domain':'prepared packet with target physical excitation mean at most m_upper+d; arbitrary purifying references',
            'ancillary_energy_itself_not_bounded':True,'full_large_carrier_executed':False}}


def certificate():
    return {'schema':'oph.factorized_native_compiler.v1','formal_fourth_order':formal_order_certificate(),
        'exact_factor_count':'6*sites even disjoint four-CAR blocks; all other kinetic factors exact',
        'complete_bound':'min(2,6*sites*[(7/3)(4/a+1)]^5*t^5/(10*n^4))',
        'domain':'homogeneous periodic kinetic source, four-phase native analyzer; all CAR states and external quantum references',
        'complete_small_carrier_control':complete_control(),'generic_three_axis_momentum_factors':generic_momentum_factorization(),'resources':resources(),
        'scope':'Same controlled compiler and complete records; global parity/inversion services and empirical ancestry retained. No passive microscopic clock or efficient local hardware implementation.',
        'verdict':'PASS_VOLUME_LINEAR_COMPLETE_COMPILER_BOUND'}
