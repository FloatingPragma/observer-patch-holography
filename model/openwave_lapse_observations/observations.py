from fractions import Fraction as F
from functools import lru_cache
import numpy as np
import sympy as s
from scipy.linalg import eigh
from model.finite_area_response.response import Interval
from model.openwave_energy_lapse.source import prepare,native,exp_interval,KAPPA,BITS,serializable
from model.openwave_native_force.force import laplacian
from model.openwave_metric_observables.clocks import sine
from model.openwave_metric_observables.waves import taylor_step,operator_tail

@lru_cache(maxsize=1)
def preparation():
    p=prepare();q=tuple(F(int(v),2**BITS)for v in p['integer_geometry'][15,:,15])
    return p,q


def clocks():
    p,q=preparation();near,far=15,23;ratio=exp_interval(q[near]-q[far]);gaps=native()[2]
    rows=[]
    for k in range(1,8):
        phases=[exp_interval(q[j])*gaps[k]for j in(near,far)]
        probabilities=[(Interval.point(1)-sine(phi))/2 for phi in phases]
        gap_ratio=phases[0]/phases[1]
        assert gap_ratio.lo<=ratio.hi and ratio.lo<=gap_ratio.hi
        assert all(0<=r.lo<=r.hi<=1 for r in probabilities)
        local=[(p,Interval.point(1)-p,*([Interval.point(0)]*6))for p in probabilities]
        joint=[{'near_outcome':a,'far_outcome':b,'probability_interval':(local[0][a]*local[1][b]).encoded()}for a in range(8)for b in range(8)]
        rows.append({'level':k,'local_gap_intervals':[v.encoded()for v in phases],
          'complete64_joint_records':joint,
          'near_to_far_frequency_ratio_interval':ratio.encoded(),
          'direct_gap_ratio_interval':gap_ratio.encoded(),
          'plus_Y_probability_intervals':[v.encoded()for v in probabilities],
          'minus_Y_probability_rule':'1-minus plus-Y','other_six_level_probabilities':['0']*6,
          'preparation_work_intervals':[ (v/2).encoded()for v in phases]})
    assert ratio.hi<1
    return {'same_complete_eight_level_K_at_both_clocks':True,'rows':rows,
      'source_time':'1 in M0 energy units; one common clock parameter',
      'prepared_q_near':str(q[near]),'prepared_q_far':str(q[far]),
      'universal_frequency_ratio_interval':ratio.encoded(),
      'frequency_ratio_decimal':[float(ratio.lo),float(ratio.hi)],
      'relative_frequency_shift_ppm_interval':[float((ratio.lo-1)*10**6),float((ratio.hi-1)*10**6)],
      'zero_ground_reference_in_K_is_part_of_declared_architecture':True,
      'all64_joint_clock_outcomes_retained':True,
      'geometry_prepared_before_clock_intervention':True,
      'clock_preparation_energy_is_recorded':True,
      'autonomous_geometry_backreaction_claimed':False,
      'laboratory_clock_or_proper_rod_calibrated':False}


def wave_matrix(q):
    q=np.asarray(q,float);n=len(q)
    if n<3 or np.max(abs(q))>1e-3:raise ValueError('finite weak geometry domain required')
    d=np.zeros((n+1,n))
    for j in range(n+1):
        if j<n:d[j,j]=1
        if j>0:d[j,j-1]=-1
    scale=np.exp(2*q)
    return scale[:,None]*(d.T@d)*scale[None,:]


def wave(q,duration=12):
    k=wave_matrix(q);e,u=eigh(k);assert e.min()>0
    f=np.sqrt(e);h=(u*f)@u.T;x=np.arange(len(q));psi=np.exp(-(x-7)**2/8+1j*.7*x);psi/=np.linalg.norm(psi)
    evolved=psi.copy()
    for _ in range(16):evolved=taylor_step(h,evolved,duration/16)
    direct=u@(np.exp(-1j*duration*f)*(u.T@psi));difference=float(np.linalg.norm(evolved-direct));assert difference<2e-13
    probabilities=.5*abs(evolved)**2
    records=[.5,*map(float,probabilities)];assert abs(sum(records)-1)<1e-13
    energy=.5*float(np.vdot(psi,h@psi).real)
    finalenergy=.5*float(np.vdot(evolved,h@evolved).real)
    measured=float(probabilities@np.diag(h))
    return {'duration':duration,'all_normal_frequencies':list(map(float,f)),
       'complete_vacuum_and_position_probabilities':records,
       'conditional_mean_position':float(2*probabilities@x),
       'source_energy_before':energy,'source_energy_after':finalenergy,
       'mean_detector_work':measured-finalenergy,
       'Taylor_operator_error_bound':str(operator_tail(3,F(duration),16)),
       'independent_spectral_vector_discrepancy':difference,
       'norm_defect':abs(sum(records)-1)}


def metric_and_curvature():
    q=s.symbols('q',real=True);density=s.diag(-s.exp(-4*q),1,1,1)
    g=(s.sqrt(-density.det())*density.inv()).applyfunc(s.simplify)
    assert g==s.diag(-s.exp(2*q),s.exp(-2*q),s.exp(-2*q),s.exp(-2*q))
    t,x,y,z=s.symbols('t x y z',real=True);coords=(t,x,y,z);field=s.Function('q')(x,y,z);eta=s.diag(-1,1,1,1)
    tensors={}
    for gamma in(0,1):
        h=s.diag(-2*field,*([-2*gamma*field]*3));trace=s.trace(eta*h);ricci=s.zeros(4)
        for a in range(4):
            for b in range(4):
                ricci[a,b]=s.expand(sum(eta[c,c]*(s.diff(h[c,b],coords[c],coords[a])+s.diff(h[c,a],coords[c],coords[b])-s.diff(h[a,b],coords[c],coords[c]))for c in range(4))/2-s.diff(trace,coords[a],coords[b])/2)
        tensor=(ricci-eta*s.trace(eta*ricci)/2).applyfunc(s.simplify);tensors[gamma]=tensor
    lap=sum(s.diff(field,v,2)for v in(x,y,z));assert s.simplify(tensors[1][0,0]-2*lap)==0
    assert tensors[1][1:,1:]==s.zeros(3)and all(tensors[1][0,i]==0 for i in(1,2,3))
    bad=tensors[0][1,1].subs(field,-1/s.sqrt(x*x+y*y+z*z)).doit().subs({x:3,y:0,z:0});assert bad!=0
    p,_=preparation();integer=p['integer_geometry'];L=laplacian(integer);gap=native()[2][1]
    rows=[]
    for site in((15,15,15),(21,15,15),(15,21,15)):
        linearG=F(-2*int(L[site]),2**BITS)
        rho=gap*exp_interval(F(int(integer[site]),2**BITS))if site==(15,15,15)else Interval.point(0)
        residual=Interval.point(linearG)-rho*(2/KAPPA)
        bound=2*F(p['residual_infinity_upper'])/KAPPA
        assert max(abs(residual.lo),abs(residual.hi))<=bound+F(1,10**28)
        rows.append({'site':list(site),'discrete_linearized_G00':str(linearG),
          'independently_defined_local_coordinate_energy_interval':rho.encoded(),
          'G00_minus_2rho_over_kappa_interval':residual.encoded()})
    return {'coefficient_density':['-exp(-4q)','1','1','1'],
       'metric':[str(v)for v in g.diagonal()],
       'lapse':'exp(q)','proper_spatial_length_per_coordinate':'exp(-q)',
       'null_coordinate_speed':'exp(2q)','linearized_tensor_00':'2*Laplacian(q)',
       'linearized_spatial_Einstein_tensor':'zero',
       'field_equation_from_energy_stationarity':'kappa*Lq+rho_q=0',
       'source_units_coupling':'8*pi*G_model=2/kappa',
       'finite_residual_rows':rows,
       'gamma_zero_countercontrol_G11_at_3_0_0':str(bad),
       'isotropic_gamma_one_wave_gluing_is_declared':True,
       'gamma_selected_by_clock_and_force_alone':False,
       'full_covariant_stress_or_nonlinear_Einstein_equations_proved':False,
       'regional_entropy_stationarity_premises_discharged':False,
       'physical_Newton_strength_or_length_selected':False}


def run():
    p,q=preparation();time=sum((exp_interval(-2*v)for v in q),Interval.point(0));delay=time-len(q);assert delay.lo>0
    rows=[]
    for t in(8,12,16):
        curved=wave(q,t);flat=wave([F(0)]*31,t)
        tv=float(np.abs(np.array(curved['complete_vacuum_and_position_probabilities'])-flat['complete_vacuum_and_position_probabilities']).sum()/2)
        assert tv>1e-7
        rows.append({'curved':curved,'flat_control':flat,'complete_probability_total_variation':tv,'conditional_mean_position_shift':curved['conditional_mean_position']-flat['conditional_mean_position']})
    return {'geometry':serializable(p),'clock_records':clocks(),'metric_and_curvature':metric_and_curvature(),
       'isotropic_piecewise_metric_ray_time_interval':time.encoded(),
       'isotropic_piecewise_metric_delay_interval':delay.encoded(),
       'isotropic_piecewise_metric_delay_decimal':[float(delay.lo),float(delay.hi)],
       'finite_wave_carrier':'vacuum plus all31 one-quantum modes; same q register; scalar field extension',
       'wave_action':'one half sum_x exp(-4q_x) phidot_x^2 minus one half sum_edges (Dphi)^2',
       'wave_frequency_square':'K=exp(2q) D^T D exp(2q)',
       'full_source':'geometry energy + all controlled native K_x + the controlled finite wave Hamiltonian',
       'wave_runs':rows,'metric_ray_equals_lattice_group_delay_claimed':False,
       'clock_and_field_services_are_declared':True,
       'source_geometry_is_constant_during_probes':True,
       'physical_gravitational_redshift_or_Einstein_completion_claimed':False,
       'verdict':'PASS_COMMON_LAPSE_CLOCK_WAVE_RECORDS'}
