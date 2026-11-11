"""Original Wilson open-box family: exact geometry, state and curvature controls."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import hashlib,json,sys
import sympy as s
HERE=Path(__file__).resolve().parent
RER = HERE.parents[2]
from model.native_wilson_maxwell import field
from model.native_wilson_maxwell.current import matrices
from model.native_wilson_maxwell_refinement import bounds

PARENT='model/native_wilson_maxwell_ground_refinement/results.json'
PARENT_SHA='e5a40a0d4e8d183c933232758f211bcd44bf81944429867da4048c184fc7b96a'


def admit():
    p=RER/PARENT;assert hashlib.sha256(p.read_bytes()).hexdigest()==PARENT_SHA
    parent=json.loads(p.read_text())
    for path,h in parent['input_sha256'].items():assert hashlib.sha256((RER/path).read_bytes()).hexdigest()==h,path
    return parent,json.loads((RER/'model/native_wilson_maxwell/results.json').read_text())


def ceil(x):return -(-x.numerator//x.denominator)


def literal_bulk():
    alpha,beta=matrices();cs=s.symbols('c1:4',real=True);ss=s.symbols('s1:4',real=True)
    C0=s.I*beta*alpha[1]
    assert C0*beta.conjugate()*C0.H==-beta
    H=4*beta
    for a,c,t in zip(alpha,cs,ss):
        B=-(beta+s.I*a)/2
        assert C0*B.conjugate()*C0.H==-B
        H+=B*(c+s.I*t)+B.H*(c-s.I*t)
    expected=sum((a*t for a,t in zip(alpha,ss)),s.zeros(4))+beta*(4-sum(cs))
    assert (H-expected).applyfunc(s.expand)==s.zeros(4)
    scalar=sum(t*t for t in ss)+(4-sum(cs))**2
    assert (H*H-scalar*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
    w=[1-c for c in cs];photon=2*sum(w)
    on_shell=3-sum(c*c for c in cs)+(4-sum(cs))**2
    assert s.expand(on_shell-1-2*photon-2*sum(w[i]*w[j] for i in range(3) for j in range(i+1,3)))==0
    d=s.Matrix([c+s.I*t-1 for c,t in zip(cs,ss)])
    C=s.Matrix([[-d[1],d[0],0],[-d[2],0,d[0]],[0,-d[2],d[1]]])
    assert (C*d).applyfunc(s.expand)==s.zeros(3,1)
    lam=(d.H*d)[0]
    assert (C.H*C-(lam*s.eye(3)-d*d.H)).applyfunc(s.expand)==s.zeros(3)
    c=s.Symbol('c',real=True);velocity_squared=4*(1-c*c)/(5-4*c)
    assert s.simplify(1-velocity_squared-(2*c-1)**2/(5-4*c))==0
    return {'literal_source_spin_matrices_and_oriented_links_replayed':True,
        'bulk_h_over_m':'alpha.sin(q)+beta*(1+sum_i(1-cos(q_i)))',
        'photon_squared_over_m_squared':'2 sum_i(1-cos(q_i))',
        'exact_squared_dispersion_identity':'EF^2-m^2=2 omega_gamma^2+2 m^2 sum_i<j(1-cos(q_i))(1-cos(q_j))',
        'axis_squared_dispersion_ratio':'2',
        'axis_fermion_max_group_velocity_at_a_m_1':'1, attained at cos(q)=1/2',
        'axis_photon_max_group_velocity':'1',
        'curvature_mismatch_is_a_maximal_signal_speed_theorem':False,
        'common_clock_rescaling_changes_curvature_ratio':False,
        'volume_refinement_equals_a_to_zero':False}


def hodge(n):
    vertices,edges,faces,D,C=field.cochains(n)
    H=D*D.T+C.T*C;expected=s.zeros(len(edges))
    def lap(size,x,y,edge):
        if x==y:return 2 if edge else int(x>0)+int(x<size-1)
        return -int(abs(x-y)==1)
    for i,(tail,head,a) in enumerate(edges):
        x=vertices[tail]
        for j,(tail2,head2,b) in enumerate(edges):
            if a!=b:continue
            y=vertices[tail2]
            expected[i,j]=sum(lap(n if axis==a else n+1,x[axis],y[axis],axis==a)
                for axis in range(3) if all(x[k]==y[k] for k in range(3) if k!=axis))
    assert H==expected
    V=(n+1)**3;E=3*n*(n+1)**2;T=E-V+1
    assert D.rank()==V-1 and C.rank()==T
    assert max(sum(abs(v) for v in (C.T*C).row(i)) for i in range(E))<=16
    return {'cells_per_axis':n,'vertices':V,'edges':E,'transverse_modes':T,
        'exact_separable_one_form_Hodge_identity':True,'exact_transverse_rank':True}


def preparation(n,p0,target=F(1,2**40)):
    V=(n+1)**3;p0=bounds.rounded(p0,False,60);f=1-p0**4
    # Vanishing failure schedule makes the local success/failure marginal
    # itself converge as boxes grow; all fine flags are retained separately.
    epsilon=target/F((V+1)**2)
    R=1
    while V*f**R>epsilon:R+=1
    return {'cells_per_axis':n,'vertices':V,'parallel_full_box_banks':R,
        'native_constituents_retained':4*V*R,'complete_binary_code_flags':V*R,
        'single_cell_failure_upper':str(bounds.rounded(f)),
        'complete_assembly_failure_upper':str(bounds.rounded(V*f**R)),
        'chosen_failure_target':str(epsilon),
        'conditional_active_state':'exact source free Wilson/code vacuum with photon vacuum and its actual Gauss reference',
        'exact_code_projector_service_required':True,'coupled_ground_prepared':False,
        'all_banks_flags_and_original_spin_references_retained':True}


def cutoff(n,error_bits=40,tau=F(1,64),apparatus_radius=F(1,16)):
    V=(n+1)**3;E=3*n*(n+1)**2;T=E-V+1;g=F(1,64)
    v=F(n+1,4)/bounds.root(2,False);r=bounds.root(g*g*v)
    # The extra V-dependent accuracy makes a local record family Cauchy.
    k=error_bits+2*(V+1).bit_length()
    x=4*E*tau
    J=max(ceil(6*x),k+3)-1
    R=J*r+apparatus_radius
    y=x+R*R
    N=4*max(8,ceil(3*R*R),k+1+2*ceil(y))
    assert N>=32 and N>=12*R*R
    assert J+1>=6*x and J+1>=k+3
    assert 2*ceil(y)-3*N//4<=-k-1
    return {'cells_per_axis':n,'vertices':V,'edges':E,'transverse_modes':T,
        'photon_cap':N-1,'photon_dimension_formula':f'binom({N-1+T},{T})',
        'single_link_coordinate_variance_upper':str(bounds.rounded(v)),
        'Dyson_weight':str(x),'Dyson_order':J,'total_word_radius_upper':str(bounds.rounded(R)),
        'complete_history_vector_error_upper':f'2^-{k}',
        'all_finite_operator_matrices_executed':False,
        'volume_independent_accuracy_after_declared_cap_choice':True}


def square_root_tail(order):
    return F(comb(2*order,order),4**order)


def infrared_constants():
    pi_upper=F(22,7)
    coefficient=pi_upper*pi_upper*bounds.root(pi_upper)/bounds.root(2,False)*(F(1,2)+pi_upper/16)
    assert coefficient<32 and bounds.root(F(32))<6
    for j in range(1,128):
        b=square_root_tail(j)
        a=F(comb(2*j,j),4**j*(2*j-1))
        assert square_root_tail(j-1)-b==a and b*b*(j+1)<=1
    return {'Neumann_dipole_heat_squared_norm_coefficient_upper':str(bounds.rounded(coefficient)),
        'uniform_sufficient_coefficient':'32','integrated_heat_tail_coefficient':'24',
        'square_root_polynomial_exact_telescope_checks':127}


def scientific():
    from . import dynamics
    parent,source=admit();p0=F(source['actual_preparation']['complete_trial_ground_overlap_lower_per_constituent'])
    g=F(1,64);row=16+8*g+4*g*g
    return {'schema':'oph.original_wilson_spatial_family.v1','parent_ground_refinement_sha256':PARENT_SHA,
        'literal_bulk_curvature_control':literal_bulk(),
        'infrared_and_covariance_constants':infrared_constants(),
        'executed_open_box_Hodge_controls':[hodge(n) for n in (1,2,3)],
        'uniform_local_majorant_row_bound':str(row),
        'uniform_path_bound':'C_A C_B exp(kappa exp(mu) tau-mu distance), kappa=16+8g+4g^2',
        'same_source_parameters':'native h, m, g, a=1/m, four constituents per vertex, c(g), full references; only box size and photon cap vary',
        'source_native_code_bank_preparation':[preparation(n,p0) for n in (1,2,3,7,15,31)],
        'constructive_photon_budgets':[cutoff(n) for n in (1,2,3,7,15,31)],
        'free_transverse_photon_gap':'(2 sqrt(2)/a) sin(pi/(2(n+1)))',
        'uniform_interacting_ground_or_ground_gap_transfer':False,
        'equal_time_local_field_record_limit':{
            'magnetic_covariance':'(C C^T)^(1/2)/2',
            'electric_transverse_covariance':'(C^T C)^(1/2)/2',
            'sqrt_polynomial_uniform_error':'sqrt(12)*binom(2J,J)/4^J <= sqrt(12)/sqrt(J+1)',
            'dipole_Green_heat_tail_norm_upper':'24 T^(-1/4), T>=1, uniformly in open 3D boxes',
            'free_Wilson_projector_local_polynomial_error_upper':'100*(99/100)^(J+1)',
            'nonzero_autonomous_time_coupled_volume_limit_proved':True},
        'coupled_positive_time_record_limit':dynamics.scientific(),
        'global_native_guard_removed':False,'microscopic_native_Lieb_Robinson_theorem':False,
        'spatial_continuum_or_Lorentz_limit_proved':False}


if __name__=='__main__':print(json.dumps(scientific(),indent=2,sort_keys=True))
