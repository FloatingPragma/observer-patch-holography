"""Complete finite scalar-wave records from the same declared coefficient density."""
from fractions import Fraction as F
from math import factorial
import numpy as np
from scipy.linalg import eigh
from .geometry import prepared


def wave_matrix(radii, outside=1.):
    r = np.asarray(radii,float); n=len(r)
    if n < 3 or np.any(r<=0) or outside<=0:
        raise ValueError('positive finite density and at least three cells required')
    incidence=np.zeros((n+1,n))
    for edge in range(n+1):
        if edge<n:incidence[edge,edge]=1
        if edge>0:incidence[edge,edge-1]=-1
    extended=np.concatenate(([outside],r,[outside]))
    stiffness=(extended[:-1]**2+extended[1:]**2)/2
    raw=incidence.T@(stiffness[:,None]*incidence)
    return np.sqrt(r)[:,None]*raw*np.sqrt(r)[None,:],raw,incidence,stiffness


def taylor_step(hamiltonian,state,duration,degree=36):
    result=state.copy()
    for k in reversed(range(1,degree+1)):
        result=state+(-1j*duration/k)*(hamiltonian@result)
    return result


def operator_tail(norm_bound,duration,steps,degree=36):
    norm_bound=F(norm_bound); duration=F(duration); x=norm_bound*abs(duration)/steps
    if steps<1 or degree<1 or x>=degree+2:raise ValueError('convergent scaled Taylor certificate required')
    single=x**(degree+1)/F(factorial(degree+1))/(1-x/F(degree+2))
    return steps*single*(1+single)**(steps-1)


def execute(radii,duration,steps=16):
    matrix,raw,incidence,stiffness=wave_matrix(radii)
    values,modes=eigh(matrix); assert values.min()>0
    frequencies=np.sqrt(values); h=(modes*frequencies)@modes.T
    n=len(radii); coordinate=np.arange(n,dtype=float)
    state=np.exp(-(coordinate-7)**2/8+1j*.7*coordinate)
    state/=np.linalg.norm(state)
    evolved=state.copy()
    for _ in range(steps):evolved=taylor_step(h,evolved,duration/steps)
    direct=modes@(np.exp(-1j*duration*frequencies)*(modes.T@state))
    assert np.linalg.norm(evolved-direct)<2e-13
    normerror=abs(float(np.vdot(evolved,evolved).real)-1)
    assert normerror<2e-13
    probabilities=np.abs(evolved)**2
    # Vacuum plus every one-quantum mode is the complete declared finite field carrier.
    records=[.5,*list(.5*probabilities)]
    assert abs(sum(records)-1)<1e-13
    bound=operator_tail(3,F(duration),steps)
    assert np.linalg.norm(h,2)<3
    return {'normal_frequencies':list(map(float,frequencies)),
            'complete_vacuum_and_position_probabilities':list(map(float,records)),
            'one_quantum_mean_position':float(coordinate@probabilities),
            'one_quantum_right_region_probability':float(probabilities[16:].sum()),
            'source_time':duration,'scaled_Taylor_steps':steps,'Taylor_degree':36,
            'analytic_evolution_operator_error_bound':str(bound),
            'independent_spectral_evolution_vector_difference':float(np.linalg.norm(evolved-direct)),
            'raw_norm_error':normerror,
            'symmetric_square_root_residual_Frobenius_diagnostic':float(np.linalg.norm(h@h-matrix))}


def certificate():
    _,radii,_=prepared(); runs=[]
    for duration in (8,12,16):
        curved=execute(radii,duration);flat=execute([F(1)]*len(radii),duration)
        difference=np.asarray(curved['complete_vacuum_and_position_probabilities'])-flat['complete_vacuum_and_position_probabilities']
        total_variation=float(np.abs(difference).sum()/2)
        assert total_variation>1e-7
        runs.append({'curved':curved,'flat_control':flat,
                     'complete_output_total_variation':total_variation,
                     'conditional_mean_position_shift':curved['one_quantum_mean_position']-flat['one_quantum_mean_position']})
    # The constant-density finite dispersion follows from the same action.
    a=F(9,10);n=7;k,_,_,_=wave_matrix([a]*n,outside=float(a))
    flat,_,_,_=wave_matrix([F(1)]*n)
    assert np.max(np.abs(k-float(a**3)*flat))<1e-14
    return {'declared_spatial_action':'1/2 sum (1/r_x) phidot_x^2 -1/2 sum_edges mean(r^2)(D phi)^2',
            'canonical_mass_weighted_frequency_square':'K=R^(1/2) D^T A D R^(1/2)',
            'density_and_stiffness_come_from_C':'C_tt=-1/r, C_transverse=r^2',
            'boundary':'zero Dirichlet scalar field; outside geometry r=1',
            'spatial_sites':31,'normal_modes_retained':31,'finite_quantum_carrier_dimension':32,
            'carrier':'vacuum and all single-quanta normal modes, with exact normal-mode Hamiltonian',
            'initial_state':'equal coherent vacuum and normalized right-moving Gaussian one-quantum packet',
            'detector':'vacuum outcome plus all 31 orthogonal site-frame one-quantum outcomes',
            'uniform_geometry_frequency_square_scaling':'r^3',
            'uniform_geometry_dispersion':'omega^2=4 r^3 sin(k/2)^2',
            'runs':runs,
            'all_modes_and_measurement_complements_retained':True,
            'analytic_Taylor_bound_includes_floating_roundoff':False,
            'field_one_quantum_carrier_equal_to_native_eight_level_carrier':False,
            'finite_group_delay_identified_with_metric_ray_time':False,
            'physical_light_or_measured_gravitational_delay_identified':False}
