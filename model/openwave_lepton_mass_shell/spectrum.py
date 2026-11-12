"""Literal native Dirac frames, all valleys and complete probe outcomes."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
import sympy as s
from scipy.linalg import block_diag,eigh
from model.lepton_mass_compiler import analyzer,compiler as native
from model.common_source import native_matter_attachment as matter
from model.common_source.chiral_candidate import kinetic
from model.openwave_kinetic_phase_cycle.compiler import rounds_for_error,trace_bound


@lru_cache(maxsize=1)
def native_frames():
    rep=matter.representation();indices=list(rep['selected'])
    weak=rep['operators']['weak_C2'].extract(indices,indices)
    left=s.diag(*(int(weak[i,i]!=0)for i in range(15)))
    row=next(row for row in matter.invariant_higgs_lines()if row['neutral_rank']==2)
    D=left*row['components'][1]*(s.eye(15)-left)
    assert D.rank()==1 and D*D==s.zeros(15) and s.trace(D.H*D)==1
    column=next(j for j in range(15)if D[:,j]!=s.zeros(15,1))
    ell=D[:,column]/s.sqrt((D[:,column].H*D[:,column])[0]);right=D.H*ell
    assert D==ell*right.H and (right.H*right)[0]==1
    basis=s.Matrix.hstack(s.conjugate(ell),right)
    assert basis.H*basis==s.eye(2)
    assert basis.H*(s.conjugate(D)+D.H)*s.conjugate(basis)==s.Matrix([[0,1],[1,0]])
    Q=(rep['operators']['T3']+rep['operators']['q']/6).extract(indices,indices)
    assert Q*basis==basis*s.diag(-1,1)
    px,py,pz=s.symbols('px py pz',real=True);masses=s.symbols('m0:3',positive=True)
    sigma=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1));epsilon=s.Matrix([[0,1],[-1,0]])
    K=s.kronecker_product(D,s.diag(*masses));K+=K.T
    delta=s.kronecker_product(epsilon,s.conjugate(K))
    h=-s.kronecker_product(sum((p*t for p,t in zip((px,py,pz),sigma)),s.zeros(2)),s.eye(45))
    H=s.Matrix.vstack(s.Matrix.hstack(h,delta),s.Matrix.hstack(delta.H,h.T))
    # lower block is -h(-p)^T=h(p)^T, retaining opposite momentum.
    charge=s.diag(s.kronecker_product(s.eye(2),s.kronecker_product(Q,s.eye(3))),
                  -s.kronecker_product(s.eye(2),s.kronecker_product(Q.T,s.eye(3))))
    small=s.zeros(8,4);small[0,0]=small[2,1]=small[7,2]=1;small[5,3]=-1
    frames=[]
    for k in range(3):
        u=s.kronecker_product(basis,s.eye(3)[:,k]);spin_u=s.kronecker_product(s.eye(2),u)
        frame=s.diag(spin_u,s.conjugate(spin_u))*small
        assert frame.H*frame==s.eye(4) and charge*frame==-frame
        expected=-s.kronecker_product(s.diag(1,-1),sum((p*t for p,t in zip((px,py,pz),sigma)),s.zeros(2)))+masses[k]*s.kronecker_product(s.Matrix([[0,1],[1,0]]),s.eye(2))
        assert (H*frame-frame*expected).applyfunc(s.simplify)==s.zeros(180,4)
        for spin in sigma:
            original=s.diag(s.kronecker_product(spin/2,s.eye(45)),-s.kronecker_product(spin.T/2,s.eye(45)))
            assert frame.H*original*frame==s.kronecker_product(s.eye(2),spin/2)
        frames.append(frame)
    return {'D':D,'Q':Q,'frames':frames,'mass_symbols':masses,'momenta':(px,py,pz),'H':H,'representation':rep}


def dirac_hamiltonian(momentum,mass):
    p=np.asarray(momentum,float)
    if p.shape!=(3,)or not np.isfinite(p).all()or not np.isfinite(mass)or mass<=0:raise ValueError('finite momentum and positive mass required')
    sig=kinetic.pauli_matrices();tauz=np.diag([1.,-1.]);taux=np.array([[0.,1.],[1.,0.]])
    return -np.kron(tauz,sum((x*matrix for x,matrix in zip(p,sig)),np.zeros((2,2),complex)))+mass*np.kron(taux,np.eye(2))


def exact_massless_charge_control():
    data=native_frames();rep=data['representation'];indices=list(rep['selected']);Q=data['Q']
    weak=rep['operators']['weak_C2'].extract(indices,indices)
    color=rep['operators']['color_C2'].extract(indices,indices)
    candidates=[i for i in range(15)if Q[i,i]==s.Rational(1,3)and weak[i,i]==0 and color[i,i]==s.Rational(4,3)]
    assert len(candidates)==3
    assert all(data['D'][:,i]==data['D'][i,:].T==s.zeros(15,1)for i in candidates)
    # The exterior cube of this complete colour module is the determinant
    # line. Every traceless colour generator acts by its trace, exactly0.
    for current in rep['matter_currents']:
        block=current.extract(candidates,candidates)
        outside=[i for i in range(15)if i not in candidates]
        assert current.extract(outside,candidates)==s.zeros(12,3)
    for generator in matter.canonical_currents()[:8]:
        current=rep['encoder'].H*matter.exterior_derivative(matter.spin_lift(generator))*rep['encoder']
        assert s.trace(current.extract(candidates,candidates))==0
    charge=sum(Q[i,i]for i in candidates);assert charge==1
    return {'actual_native_indices':candidates,'native_electric_charge':'1',
        'source_energy_at_zero_lattice_momentum':'0','fermion_number':3,
        'colour_singlet_exterior_determinant_line':True,'exact_native_colour_generator_trace_checks':8,
        'needs_three_physical_zero_modes_not_a_new_particle_coordinate':True,
        'total_charge_one_sector_has_positive_mass_floor':False,
        'interpretation':'A three-quark zero-mode state of the declared charged-lepton-only free action; no bound baryon interpretation.'}


def all_lattice_shells(side=4):
    if side!=4:raise ValueError('the complete exact sine census uses periodic side4')
    boxes=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']
    masses=np.array([float(sum(map(F,box),F(0))/2)for box in boxes])
    data=native_frames();D=np.asarray(data['D'],complex);K=np.kron(D,np.diag(masses));K+=K.T
    exact_sin=(0,1,0,-1);rows=[];ground=0.;zero_count=0
    for address in product(range(4),repeat=3):
        momentum=tuple(exact_sin[n]for n in address);norm2=sum(x*x for x in momentum)
        H=kinetic.bdg_symbol(momentum,K,speed=1);values=np.linalg.eigvalsh(H)
        positive=np.r_[np.repeat(np.sqrt(norm2+masses*masses),4),np.repeat(np.sqrt(norm2),78)]
        expected=np.sort(np.r_[-positive,positive])
        assert len(values)==180 and np.max(abs(values-expected))<2e-13
        if norm2==0:zero_count+=78
        ground-=sum(positive)/2
        for k,frame in enumerate(data['frames']):
            U=np.asarray(frame,complex);block=dirac_hamiltonian(momentum,masses[k])
            assert np.linalg.norm(H@U-U@block)<1e-13
        rows.append({'momentum_grid_address':list(address),'sin_momentum':list(momentum),
            'kinetic_norm_squared':norm2,'complete_180_Nambu_energies':values.tolist(),
            'charged_branch_positive_energies':np.sqrt(norm2+masses*masses).tolist(),
            'native_frame_intertwiner_maximum_error':float(max(np.linalg.norm(H@np.asarray(frame,complex)-np.asarray(frame,complex)@dirac_hamiltonian(momentum,masses[k]))for k,frame in enumerate(data['frames'])))})
    assert zero_count==624
    return {'periodic_side':4,'spatial_sites':64,'physical_CAR_modes':5760,'complete_Nambu_rows':11520,
            'three_distinct_positive_branch_minima_k0_k1_k2':masses.tolist(),
            'rest_valley_count':8,'charged_Dirac_fields_per_valley':3,'charged_family_valley_copies':24,
            'massless_internal_coefficients_retained':39,'zero_energy_physical_CAR_modes':zero_count,
            'full_ground_reference_degeneracy':'2^624, from retained massless zero modes',
            'chosen_reference':'all negative modes occupied; zero modes vacant; exact same quadratic Hamiltonian',
            'ground_energy_in_unshifted_pairing_convention_diagnostic':float(ground),
            'same_source_single_excitation_energy':'sqrt(sum_i sin(ap_i)^2/a^2+e_k^2)',
            'rest_minimum_Hessian_in_each_valley':'delta_ij/e_k at spacing1 and speed1',
            'rows':rows,'exactly_three_physical_species_inferred':False}


def spatial_records(mass,time,side=4):
    # Positive-band creation at one actual spatial source coordinate.
    points=list(product(range(side),repeat=3));volume=side**3
    matrices=[];projectors=[];momenta=[]
    for address in points:
        p=2*np.pi*np.array(address)/side;momenta.append(p)
        sinp=np.sin(p);sinp[np.abs(sinp)<1e-14]=0
        h=dirac_hamiltonian(sinp,mass);energy=np.sqrt(sinp@sinp+mass*mass)
        matrices.append(h);projectors.append((np.eye(4)+h/energy)/2)
    f=np.array([1.,0,0,0],complex)
    initial=np.array([P@f/np.sqrt(volume)for P in projectors]);weight=float(np.vdot(initial,initial).real)
    assert abs(weight-.5)<3e-15
    evolved=[]
    for h,phi in zip(matrices,initial):
        energy=np.sqrt(float(np.trace(h@h).real/4))
        evolved.append((np.cos(energy*time)*np.eye(4)-1j*np.sin(energy*time)*h/energy)@phi)
    real=np.zeros((volume,4),complex)
    for x,point in enumerate(points):
        real[x]=sum((np.exp(1j*np.dot(p,point))*v for p,v in zip(momenta,evolved)),np.zeros(4,complex))/np.sqrt(volume)
    original=np.abs(real)**2
    assert abs(original.sum()-weight)<3e-15
    return {'source_time':time,'reservoir_failure_probability':1-weight,
        'unconditional_all_site_spinor_arrival_probabilities':original.tolist(),
        'original_input_probability_sum':float(1-weight+original.sum()),
        'arrival_effects':'positive-band-filtered source CAR effects, summed over all256 site/spinor outcomes; full source complement retained',
        'reference':'filled negative bands of this source, explicit vacant zero spectators; charged reservoir and opposite idle reference retained',
        'raw_source_local_position_measurement_claimed':False,
        'quantum_failure_output_retained':True}


def complete_rest_records():
    instrument=native.exact_rest_instrument();rows=[]
    for character in range(3):
        for valley in product((0,1),repeat=3):
            for outcome in instrument['complete_records']:
                row=dict(outcome);row.update({'character':character,'valley':list(valley),
                    'probability_exact':str(s.Rational(outcome['probability_exact'])/24)})
                rows.append(row)
    assert sum(s.Rational(row['probability_exact'])for row in rows)==1
    return {'complete_setting_reservoir_energy_charge_records':rows,'record_count':864,
            'zero_probability_records':sum(row['probability_exact']=='0'for row in rows),
            'all_other_family_valley_and_spectator_modes_retained':True,
            'source_energy_unit':'energy_in_units_of_same_m times native e_character',
            'declared_spectral_measurement_service':True}


def certificate():
    frames=native_frames();shell=all_lattice_shells();control=exact_massless_charge_control()
    masses=shell['three_distinct_positive_branch_minima_k0_k1_k2'];spatial=[]
    for k,mass in enumerate(masses):
        records=[spatial_records(mass,t)for t in (0,1,3)]
        changes=[]
        initial=np.asarray(records[0]['unconditional_all_site_spinor_arrival_probabilities'])
        for row in records[1:]:
            change=float(np.abs(np.asarray(row['unconditional_all_site_spinor_arrival_probabilities'])-initial).sum()/2)
            assert change>1e-4;changes.append(change)
        spatial.append({'character':k,'records':records,'arrival_total_variation_from_time_zero':changes})
    n=rounds_for_error(F(1,10**12),sites=64)
    return {'schema':'oph.native_lepton_mass_shell.v1',
        'actual_native_Dirac_frame':{'complete_native_Nambu_rows':180,'families':3,
            'exact_Hamiltonian_intertwiners':3,'exact_native_spin_intertwiners':9,
            'same_native_electric_charge_on_each_particle_frame':'-1',
            'other_signed_charge_frames_and_spectators_retained_in_full_spectrum':True},
        'complete_lattice_spectrum':shell,'complete_rest_instrument':complete_rest_records(),
        'spatial_complete_record_experiments':spatial,'massless_charge_countercontrol':control,
        'source_calibration':native.electron_calibration(),
        'finite_compiler_transfer':{'sites':64,'cycles':n,'unhalved_trace_error':str(trace_bound(n,sites=64)),
            'same_full_source_state_reference_and_measurements':True,'entire_output_probability_TV_bound':'half the unhalved trace bound',
            'all_nonvacuum_ancillary_outcomes_retained':True,'whole_enlarged_volume_execution_claimed':False},
        'scope':'Three distinct free charged-lepton branch masses with exact native frame/charge/spin attachment, all eight valleys and39 massless internal spectators. The coefficients and electron anchor are inherited; phase2/9 remains a retrospective input. No exactly-three observed species, interacting vacuum, autonomous mass minima across all charged sectors, canonical physical Yukawa normalization or measured physical clock is inferred.',
        'verdict':'PASS_NATIVE_MASSIVE_BRANCHES_AND_COMPLETE_SPATIAL_RECORDS'}
