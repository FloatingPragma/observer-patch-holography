"""Compact smooth native packets, rigorous kinetic moments and full records."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import numpy as np
from model.lepton_mass_compiler import analyzer,compiler as native
from model.openwave_kinetic_phase_cycle.compiler import rounds_for_error,trace_bound
from model.openwave_lepton_mass_shell.spectrum import native_frames,dirac_hamiltonian


def pi_interval(terms=60):
    def atan(q):
        partial=sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(terms)),F())
        next_term=F((-1)**terms,(2*terms+1)*q**(2*terms+1))
        return min(partial,partial+next_term),max(partial,partial+next_term)
    a,b=atan(5),atan(239)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]


def sqrt_bounds(value,digits=60):
    value=F(value)
    if value<0:raise ValueError('nonnegative radicand required')
    scale=10**digits;n=isqrt(value.numerator*scale**2//value.denominator)
    low=F(n,scale);high=F(n+1,scale)
    assert low*low<=value<=high*high
    return low,high


def one_axis_moments(length):
    """Exact expressions evaluated numerically, with rational global bounds below."""
    if type(length)is not int or length<4:raise ValueError('compact envelope length at least4')
    q=np.pi/(length+1);A2=8/(3*(length+1));sq=np.sin(q);s2=np.sin(2*q)
    second=s2*s2/3+4*sq**4/(3*(length+1))
    fourth=A2*(length*s2**4/8-s2*s2*np.cos(2*q)*sq*sq/2+sq**4/4)
    return second,fourth


def moment_upper_bounds(length):
    if type(length)is not int or length<4:raise ValueError('compact envelope length at least4')
    q=pi_interval()[1]/(length+1)
    # cos(2q)>=0 for length>=4, so the exact negative cross term can be dropped.
    d2=4*q*q/3+4*q**4/(3*(length+1))
    d4=F(8,3*(length+1))*(2*length*q**4+q**4/4)
    return {'one_axis_second':d2,'one_axis_fourth':d4,
            'kinetic_second':3*d2,'kinetic_fourth':3*d4+6*d2*d2}


def electron_mass_interval():
    return tuple(map(F,analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2'][1]))


def certify(length=2**30,side=None,time=10**9):
    side=4*(length+1)if side is None else side
    if type(side)is not int or side<length+5 or type(time)is not int or time<0:
        raise ValueError('periodic box with nonoverlapping two-step collar and nonnegative integral time required')
    lo,hi=electron_mass_interval();mom=moment_upper_bounds(length);k2=mom['kinetic_second'];k4=mom['kinetic_fourth']
    eta=k2/(4*lo*lo)
    if eta>=1:raise ValueError('chosen packet does not certify positive source preparation')
    wmin=1-eta;d=k2/(2*lo*wmin);b2=k4/(4*lo*lo*wmin);b=sqrt_bounds(b2)[1]
    outside0=sqrt_bounds(eta/wmin)[1];outside=min(F(1),(outside0+time*b)**2)
    phase_error=time*d
    if phase_error>=1:raise ValueError('requested phase window has no coherent lower bound')
    frequency_error=d/(1-phase_error)
    anchor=F('0.00051099895069');relative=d/lo
    # Complete bin detector E <= m(1+10^-6); Markov on its nonnegative excess.
    high_energy_probability=min(F(1),relative*10**6)
    energy_range=F(282*side**3)
    compiler_error=min(F(1,10**12),d/(10**6*energy_range))
    resources=rounds_for_error(compiler_error,F(time),sites=side**3)
    actual_error=trace_bound(resources,F(time),sites=side**3)
    compiled_energy_error=energy_range*actual_error
    assert compiled_energy_error<=d/10**6
    qlo,qhi=(p/(length+1)for p in pi_interval())
    k2lower=(2*qlo-(2*qhi)**3/6)**2
    assert hi<1 and k2lower/4>compiled_energy_error
    # E(p)<2, w<=1 and a2>=sin^2(2q) keep the compiled mean above m.
    mean_excess_lower=k2lower/4
    return {'length':length,'periodic_side':side,'spatial_sites':side**3,
        'physical_CAR_modes':90*side**3,'enlarged_CAR_modes':324000*side**3,
        'electron_mass_interval':list(map(str,(lo,hi))),
        'exact_compact_envelope':'sqrt(8/[3(L+1)]) sin^2(pi*j/(L+1)), j=1..L, zero elsewhere; tensor product in three dimensions',
        'moment_upper_bounds':{k:str(v)for k,v in mom.items()},
        'preparation_success_interval':[str(wmin),'1'],
        'positive_projection_failure_upper':str(eta),
        'mean_source_energy_excess_upper':str(d),
        'mean_source_energy_excess_lower':str(mean_excess_lower),
        'mean_relative_rest_energy_excess_upper':str(relative),
        'calibrated_rest_energy_GeV':str(anchor),
        'calibrated_mean_energy_GeV_interval':[str(anchor),str(anchor*(1+relative))],
        'calibrated_mean_energy_GeV_interval_decimal':[float(anchor),float(anchor*(1+relative))],
        'original_probability_inside_compact_cube_initial_lower':str(1-2*eta),
        'conditional_initial_outside_cube_upper':str(eta/wmin),
        'source_time_horizon':time,'rms_energy_excess_upper':str(b),
        'conditional_outside_cube_uniform_time_upper':str(outside),
        'original_probability_inside_cube_uniform_time_lower':str(wmin*(1-outside)),
        'conditional_probability_above_rest_plus_one_ppm_upper':str(high_energy_probability),
        'autonomous_clock':{'reference_Hamiltonian':'0 on retained idle charge-reference carrier, inherited from the lepton source',
            'complex_return_error_from_exp_minus_i_mt_upper':str(phase_error),
            'visibility_uniform_lower':str(1-phase_error),
            'instantaneous_angular_frequency_absolute_error_upper':str(frequency_error),
            'instantaneous_angular_frequency_relative_error_upper':str(frequency_error/lo),
            'source_rest_cycles_at_horizon_interval':[str(lo*time/(2*pi_interval()[1])),str(hi*time/(2*pi_interval()[0]))],
            'autonomous_law':'declared quadratic physical target Hamiltonian; the finite analyzer phase-cycle implementation uses continuing controls',
            'instantaneous_frequency_derivative_transferred_to_controlled_pulses':False,
            'coherent_preparation_required':True,'passive_self_start_proved':False,'rigidity_to_differential_native_scalar_drive':False},
        'finite_full_Fock_compiler':{'cycles':resources,'full_source_unhalved_trace_bound':str(actual_error),
            'physical_excitation_energy_operator_norm_upper':str(energy_range),
            'conditional_on_initial_success_physical_energy_mean_transfer_error_upper':str(compiled_energy_error),
            'compiled_physical_mean_energy_GeV_interval':[str(anchor),str(anchor*(1+(d+compiled_energy_error)/lo))],
            'energy_observable_scope':'physical target Hamiltonian minus its chosen physical ground, tensored with identity on all ancillary outputs; no final ancillary postselection or enlarged-source ground claim',
            'all_ancillary_outputs_and_arbitrary_retained_references':True,
            'full_enlarged_volume_time_evolution_executed':False},
        'complete_original_coarse_instrument':{'preparation_failure':'1-w','inside_cube':'w p_C(t)','outside_cube':'w (1-p_C(t))','sum':'1'},
        'complete_energy_instrument':{'preparation_failure':'1-w','near_rest_bin':'w p_[m,m(1+1e-6)]','all_other_energies':'w (1-p_[m,m(1+1e-6)])','sum':'1'},
        'all_position_records':'(1-w) plus |N^(-3/2) sum_p exp(ip.x-itE_p) P_+(p) fhat(p)|_alpha^2 for every x and all four native spinor coordinates; complete source complement retained',
        'full_finite_record_law_has_been_enumerated':False}


def finite_execution(length=6,side=16,time=3,mass=None):
    if side<length+5:raise ValueError('nonoverlapping support collar required')
    mass=float(sum(electron_mass_interval(),F())/2)if mass is None else mass
    f=np.zeros(side);j=np.arange(1,length+1);f[j]=np.sqrt(8/(3*(length+1)))*np.sin(np.pi*j/(length+1))**2
    assert abs(f@f-1)<2e-15
    scalar=f[:,None,None]*f[None,:,None]*f[None,None,:]
    fhat=np.fft.fftn(scalar,norm='ortho')
    grid=np.indices((side,side,side));p=np.sin(2*np.pi*grid/side);p[np.abs(p)<1e-14]=0
    p2=np.sum(p*p,axis=0);E=np.sqrt(p2+mass*mass)
    seed=np.array([1.,0,1.,0],complex)/np.sqrt(2)
    # Exact native four-row frame; vectorized H action on this beta+ seed.
    Hseed=np.empty((side,side,side,4),complex)
    Hseed[...,0]=(mass-p[2])/np.sqrt(2);Hseed[...,1]=-(p[0]+1j*p[1])/np.sqrt(2)
    Hseed[...,2]=(mass+p[2])/np.sqrt(2);Hseed[...,3]=(p[0]+1j*p[1])/np.sqrt(2)
    positive=fhat[...,None]*(seed+Hseed/E[...,None])/2
    w=float(np.sum(abs(positive)**2));weights=np.sum(abs(positive)**2,axis=-1)
    evolved=np.exp(-1j*time*E)[...,None]*positive
    real=np.fft.ifftn(evolved,axes=(0,1,2),norm='ortho');prob=np.abs(real)**2
    assert abs(prob.sum()-w)<1e-14
    second=float(np.sum(abs(fhat)**2*p2));fourth=float(np.sum(abs(fhat)**2*p2*p2))
    d2,d4=one_axis_moments(length)
    assert abs(second-3*d2)<3e-14 and abs(fourth-(3*d4+6*d2*d2))<3e-14
    C=complex(np.sum(weights*np.exp(-1j*E*time))/w)
    clock=[]
    for quadrature in ('X','Y'):
        contrast=C.real if quadrature=='X'else-C.imag
        for sign in (-1,1):clock.append({'quadrature':quadrature,'outcome':sign,'original_probability':w*(1+abs(C)**2+2*sign*contrast)/8})
        clock.append({'quadrature':quadrature,'outcome':'orthogonal_return_complement','original_probability':w*(1-abs(C)**2)/4})
    clock.append({'outcome':'preparation_failure','original_probability':1-w})
    assert abs(sum(r['original_probability']for r in clock)-1)<1e-14
    inside=float(prob[1:length+1,1:length+1,1:length+1].sum())
    return {'length':length,'periodic_side':side,'source_mass':mass,'source_time':time,
        'original_preparation_failure':1-w,'original_probability_inside_cube':inside,
        'original_probability_outside_cube':w-inside,'complete_position_spinor_probabilities':prob.reshape((-1,4)).tolist(),
        'complete_energy_momentum_probabilities':weights.reshape(-1).tolist(),'energy_momentum_values':E.reshape(-1).tolist(),
        'kinetic_second':second,'kinetic_fourth':fourth,'conditional_mean_energy':float(np.sum(E*weights)/w),
        'conditional_complex_return':[C.real,C.imag],'complete_clock_records':clock,
        'positive_projection_failure_finite_identity':float((1-mass*np.sum(abs(fhat)**2/E))/2)}


def perturbation_controls(certificate):
    t=F(certificate['source_time_horizon']);b=F(certificate['rms_energy_excess_upper']);eta=F(certificate['positive_projection_failure_upper'])
    w=1-eta;initial=sqrt_bounds(eta/w)[1]
    kick=F(1,10**5);delta=F(1,10**14)
    outside=min(F(1),(initial+t*b+kick+t*delta)**2)
    return {'preparation_vector_kick_upper':str(kick),'additional_complete_source_Hamiltonian_operator_norm_upper':str(delta),
        'same_initial_reference_and_time':True,'conditional_outside_cube_upper':str(outside),
        'source_perturbation_scope':'A bounded Hermitian perturbation of the complete source with this total operator norm; no per-site bound is substituted. State-dependent mean-field or extensive local perturbations need their own bound.',
        'reference_trace_error_after_perturbation_upper':str(2*(kick+t*delta)),
        'no_shrinkage_or_shape_minimization_inferred':True}


def clock_controls():
    # Complete source/reference instrument on a vacuum plus two particle lines.
    E=np.array([.1,.1001]);weights=np.array([.4,.6]);times=(0.,10.,1000.)
    rows=[]
    for t in times:
        C=np.sum(weights*np.exp(-1j*E*t));drive=.003;reference=.007
        initial=np.r_[1.,np.sqrt(weights)]/np.sqrt(2)
        U=np.diag(np.exp(-1j*np.r_[0.,E]*t));rho=np.outer(initial,initial.conj())
        dephased=rho.copy();dephased[0,1:]=0;dephased[1:,0]=0
        evolved=U@dephased@U.conj().T
        X=np.zeros((3,3),complex);X[0,1:]=np.sqrt(weights);X[1:,0]=np.sqrt(weights)
        Y=np.zeros((3,3),complex);Y[0,1:]=1j*np.sqrt(weights);Y[1:,0]=-1j*np.sqrt(weights)
        dephased_contrast=[float(np.trace(evolved@A).real)for A in (X,Y)]
        assert max(abs(x)for x in dephased_contrast)<1e-15
        echo=complex(np.vdot(initial,U.conj().T@U@initial))
        assert abs(echo-1)<1e-15
        rows.append({'time':t,'free_return':[C.real,C.imag],
            'differential_native_charge_drive_return':[(C*np.exp(1j*drive*t)).real,(C*np.exp(1j*drive*t)).imag],
            'nonzero_reference_energy_return':[(C*np.exp(-1j*reference*t)).real,(C*np.exp(-1j*reference*t)).imag],
            'same_source_then_actual_inverse_return':[echo.real,echo.imag],
            'dephased_preparation_XY_contrast':dephased_contrast})
    return {'independent_finite_two_line_controls':rows,
        'differential_drive_identity':'H_drive(t)=v(t) Q_source; packet charge-1 gives C_drive=exp(i integral v) C',
        'reference_identity':'reference of opposite charge has cost mu when source particle is present; C_reference=exp(-i mu t) C',
        'inverse_identity':'the compiler source inversion plus charged-leg parity realizes the negative full generator; U_inverse U=I',
        'coherent_herald':'K_success=|0><0|_pointer sqrt(w) I + |1><1|_pointer c^dagger(f) d^dagger_reference; its full-carrier failure complement is retained',
        'no_unconditional_clock_from_diagonal_preparation':True}


def certificate():
    native_frames();base=certify();volume=certify(side=8*(2**30+1));coarser=certify(length=2**29,time=10**8)
    assert base['moment_upper_bounds']==volume['moment_upper_bounds']
    records=[finite_execution(time=t)for t in (0,1,3)]
    control=finite_execution(length=6,side=16,time=3,mass=1.)
    return {'schema':'oph.localized_native_rest_packet.v1','near_rest_localized_source':base,
        'doubled_volume_control':volume,'changed_width_control':coarser,
        'complete_executed_small_box_native_records':records,'independent_numerical_mass_control':control,
        'perturbed_persistence':perturbation_controls(base),'clock_reference_drive_inverse_controls':clock_controls(),
        'calibration':native.electron_calibration(),
        'scope':'Compact preparation and rigorous positive-band/energy/phase bounds on one declared periodic native kinetic source. All native fields, valleys, massless spectators and reference outputs are retained. The large carrier is certified by finite analytic identities, not enumerated. No stationary soliton, interacting electron, passive clock onset or scalar-drive rigidity is inferred.',
        'verdict':'PASS_LOCALIZED_NEAR_REST_PERSISTENCE_AND_PHASE_CLOCK'}
