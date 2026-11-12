"""Complete quantum angular, energy, field and localization observations."""
from fractions import Fraction as F
import numpy as np
import sympy as sp
from scipy.linalg import expm,eigh
from .source import source,photon,ladders,matter_exact,native_reduction,cross,radial_effect,MASS,ETA
from .enclosure import certify,rational_bounds


def exact_rotation():
    m=matter_exact();zero=sp.zeros(20);errors=[]
    def vanishes(matrix):return all(sp.simplify(x)==0 for x in matrix)
    assert vanishes(m['B']**2*(m['B']**2-sp.Rational(5,2)*sp.eye(20)))
    probe_square=m['external'][2]**2
    assert vanishes(probe_square*(probe_square-sp.eye(20)/4)*(probe_square-sp.eye(20)/2))
    for i in range(3):
        assert vanishes(m['J'][i]*m['B']-m['B']*m['J'][i])
        for j in range(3):
            expected=sp.I*sum((cross(i,j,k)*m['derivative'][k]for k in range(3)),zero)
            assert vanishes(m['J'][i]*m['derivative'][j]-m['derivative'][j]*m['J'][i]-expected)
    for cap in (1,2,3):
        basis,a=ladders(cap,True);n=len(basis);zero=sp.zeros(n)
        J=tuple(sum((-sp.I*cross(i,j,k)*a[j].H*a[k]for j in range(3)for k in range(3)),zero)for i in range(3))
        for i in range(3):
            for j in range(3):
                assert vanishes(J[i]*J[j]-J[j]*J[i]-sp.I*sum((cross(i,j,k)*J[k]for k in range(3)),zero))
                assert vanishes(J[i]*(a[j]+a[j].H)-(a[j]+a[j].H)*J[i]-sp.I*sum((cross(i,j,k)*(a[k]+a[k].H)for k in range(3)),zero))
    return {'complete_matter_vector_identities':12,'complete_photon_caps_checked':[1,2,3],
      'method':'exact radicals for component covariance; scalar joint coupling and functional calculus give full-H rotational invariance',
      'full_native_intertwiner':native_reduction()}


def full_fock_order_control():
    # A positive operator-valued2x2block matrix. Transposing only its matter
    # indices loses positivity, so the order in the CAR lift is essential.
    v=np.array([1,0,0,1],complex)/np.sqrt(2);k=np.outer(v,v.conj())
    c=[]
    for i in range(2):
        a=np.zeros((4,4),complex)
        for state in range(4):
            if state>>i&1:a[state^(1<<i),state]=(-1)**((state&((1<<i)-1)).bit_count())
        c.append(a)
    good=np.zeros((8,8),complex);bad=good.copy()
    for i in range(2):
        for j in range(2):
            good+=np.kron(c[i].conj().T@c[j],k[2*i:2*i+2,2*j:2*j+2])
            bad+=np.kron(c[i].conj().T@c[j],k[2*j:2*j+2,2*i:2*i+2])
    low=float(np.linalg.eigvalsh(good)[0]);wrong=float(np.linalg.eigvalsh(bad)[0]);assert low>=-1e-14 and wrong<-.49
    return {'full_two_mode_Fock_dimension':8,'correct_order_minimum':low,'transposed_field_block_mutation_minimum':wrong}


def complete_records(cap):
    d=source(cap);v=d['psi'];nf=len(d['photon']['basis']);mat=matter_exact();ball=radial_effect(5)
    prob=[]
    for n in range(cap+1):
        photon=np.diag([int(sum(x)==n)for x in d['photon']['basis']])
        inside=np.kron(np.kron(np.eye(2),ball),photon)
        number=np.kron(np.eye(20),photon)
        for label,effect in [('inside_radius_five',inside),('outside_radius_five',number-inside)]:
            p=float(np.vdot(v,effect@v).real);assert p>=-1e-13
            prob.append({'photons':n,'position_record':label,'conditional_probability':max(0,p)})
    success=abs(np.vdot(d['seed'],v))**2
    for row in prob:row['original_input_probability']=success*row['conditional_probability']
    original=prob+[{'position_record':'ground_filter_failure','original_input_probability':float(1-success)}]
    assert abs(sum(x['original_input_probability']for x in original)-1)<1e-12
    photon_prob=[float(np.vdot(v,np.kron(np.eye(20),np.diag([int(sum(x)==n)for x in d['photon']['basis']]))@v).real)for n in range(cap+1)]
    jmf=[float(np.vdot(v,j@v).real)for j in d['Jm']];jff=[float(np.vdot(v,j@v).real)for j in d['Jf']]
    r2=np.kron(np.kron(np.eye(2),np.array(mat['r2'],complex)),np.eye(nf))
    current=float(np.vdot(v,d['current_derivative']@v).real)
    fieldmeans=[float(np.vdot(v,np.kron(np.eye(20),q)@v).real)for q in d['photon']['q']]
    # Complete photon number readout changes interaction energy, even on a
    # stationary ground. Its work is accounted for separately.
    post_energy=0.
    for n in range(cap+1):
        projection=np.kron(np.eye(20),np.diag([int(sum(x)==n)for x in d['photon']['basis']]))
        part=projection@v;post_energy+=float(np.vdot(part,d['H']@part).real)
    return {'cap':cap,'evaluated_reduced_rows':len(v),'complete_native_one_particle_rows':8*len(v),
      'photon_basis_dimension':nf,'fermion_modes':160,'full_fermion_Fock_dimension':'2^160',
      'stationary_ground_energy_native':d['E'],'J_squared_eigenvalue':.75,'Jz_eigenvalue':.5,
      'Casimir_residual':float(np.linalg.norm(d['J2']@v-.75*v)),'Jz_residual':float(np.linalg.norm(d['J'][2]@v-.5*v)),
      'matter_angular_expectation':jmf,'Maxwell_angular_expectation':jff,'Maxwell_amplitude_expectation':fieldmeans,
      'mean_radius_squared':float(np.vdot(v,r2@v).real),'mean_photon_number':float(np.vdot(v,d['Nf']@v).real),
      'photon_number_probabilities':photon_prob,'source_current_probe_derivative':current,
      'static_g_with_same_rest_energy_and_total_J':4*d['E']*abs(current),
      'complete_preparation_position_photon_records':original,'ground_filter_success':float(success),
      'complete_angular_records':[{'J_squared':.75,'Jz':.5,'conditional_probability':1.,'original_input_probability':float(success)},
                                  {'other_joint_angular_outcomes':True,'conditional_probability':0.,'original_input_probability':0.},
                                  {'ground_filter_failure':True,'original_input_probability':float(1-success)}],
      'energy_parts':{'positive_Dirac_kernel':float(np.vdot(v,(d['K_minus_mass']+MASS*np.eye(len(v)))@v).real),
          'oscillator_trap':float(np.vdot(v,d['Nm']@v).real),'Maxwell_above_vacuum':float(np.vdot(v,d['Nf']@v).real),
          'idle_compensator':0.,'photon_number_measurement_work':post_energy-d['E']},
      'numerical_scope':'Floating evaluations of the exact isolated line; the separate interval receipt bounds energy, static current and vector error.'}


def quantum_dynamics_control():
    d=source(1);h=d['H'];n=len(h);rng=np.random.default_rng(511042)
    psi=rng.normal(size=(n,2))+1j*rng.normal(size=(n,2));psi/=np.linalg.norm(psi)
    u=expm(-.7j*h);after=u@psi
    torque_m=np.array([sum((1j*np.vdot(psi[:,k],(h@j-j@h)@psi[:,k])).real for k in range(2))for j in d['Jm']])
    torque_f=np.array([sum((1j*np.vdot(psi[:,k],(h@j-j@h)@psi[:,k])).real for k in range(2))for j in d['Jf']])
    drift=max(abs(sum(np.vdot(after[:,k],j@after[:,k])-np.vdot(psi[:,k],j@psi[:,k])for k in range(2)))for j in d['J'])
    algebra=max(np.linalg.norm(h@j-j@h)for j in d['J']);assert algebra<1e-11 and drift<1e-12
    assert np.linalg.norm(torque_m)>1e-5 and np.linalg.norm(torque_m+torque_f)<1e-13
    f=d['photon'];eta=float(ETA);raw=sum((eta**2*p@p+q@q/(4*eta**2)for p,q in zip(f['pi'],f['q'])),np.zeros_like(f['N']))-1.5*np.eye(len(f['basis']))
    boundary=(1+3)*f['top']/2
    assert np.linalg.norm(raw+boundary-f['N'])<1e-14 and np.linalg.eigvalsh(raw)[0]<0
    return {'reference_dimension':2,'full_entangled_input_evolved':True,'unitary_norm_error':abs(float(np.linalg.norm(after))-1),
       'total_generator_commutator_residual':float(algebra),'total_angular_drift':float(drift),
       'matter_torque':torque_m.tolist(),'field_torque':torque_f.tolist(),
       'frozen_field_torque_mutation':float(np.linalg.norm(torque_m)),
       'omitted_Maxwell_boundary_completion_minimum':float(np.linalg.eigvalsh(raw)[0])}


def run():
    rotations=exact_rotation();rows=[complete_records(c)for c in(1,2,3)];bounds=[certify(c)for c in(1,2,3)]
    anchor=0.00051099895069;elo,ehi=bounds[1]['actual_ground_energy_interval'];exact_anchor=F('0.00051099895069');unit=[float(np.nextafter(float(exact_anchor/F.from_float(ehi)),-np.inf)),float(np.nextafter(float(exact_anchor/F.from_float(elo)),np.inf))]
    comparisons=[]
    for a,b in zip(bounds,bounds[1:]):
        e0,e1=a['actual_ground_energy_interval'],b['actual_ground_energy_interval'];g0,g1=a['rest_energy_normalized_static_g_interval'],b['rest_energy_normalized_static_g_interval']
        comparisons.append({'caps':[a['cap'],b['cap']],'ground_energy_change_interval':[e1[0]-e0[1],e1[1]-e0[0]],
            'static_g_change_interval':[g1[0]-g0[1],g1[1]-g0[0]],'scope':'finite complete carriers; neither an infinite-photon convergence theorem nor physical uncertainty'})
    for row,bound in zip(rows,bounds):
        lo,hi=bound['actual_ground_energy_interval'];row['common_cap_two_electron_unit_GeV_interval']=[float(np.nextafter(lo*unit[0],-np.inf)),float(np.nextafter(hi*unit[1],np.inf))]
    return {'schema':'oph.quantum_rotational_particle.v1','verdict':'PASS_QUANTUM_CASIMIR_CALIBRATED_RECORDS_AND_FINITE_CAP_CONTROLS',
      'declared_source':{'mass':16,'charge':-1,'spatial_total_shell':2,'trap':1,'Maxwell_frequency':1,'field_radius_squared':'5/2','field_amplitude_scale':'1/8',
        'Maxwell_kappa':'c0/32, with c0=pi^(3/2)R^5/4','photon_total_occupation_caps':[1,2,3],
        'native_blocks':4,'signed_Dirac_branches_retained':True,'all_CAR_sectors_retained':True,'physical_local_Gauss_law_or_continuum_QED':False},
      'exact_rotation_algebra':rotations,'operator_valued_Fock_positivity_control':full_fock_order_control(),
      'linear_reference_uniform_dynamics':quantum_dynamics_control(),'stationary_records':rows,'interval_certificates':bounds,
      'finite_cap_comparisons':comparisons,'calibration':{'electron_anchor_GeV':anchor,'anchor_member_photon_cap':2,
        'GeV_per_native_energy_interval':unit,'other_caps_reuse_same_unit':True,'physical_electron_identity_proved':False,
        'muon_tau_spectrum_transferred':False},
      'scope':'Complete finite quantum matter/Maxwell construction, actual dressed half-spin Casimir and source energy/current records. The rotational regulator, positive kernel, trap, charge assignment, coupling, preparation/filter services and one empirical electron unit are declared.'}
