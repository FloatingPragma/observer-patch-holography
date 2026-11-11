"""Charged local histories and receding-reference controls, at fixed spacing."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys
import sympy as s
HERE=Path(__file__).resolve().parent;WORK=HERE.parent
RER=HERE.parents[2]
from model.native_wilson_maxwell import field,matter
from model.native_wilson_maxwell_refinement import bounds
from model.particle_common_contract.provenance import import_closure
PARENTS={
 'native_orthogonal_wilson_spatial_join':'44ca365fd77c85e1d820d2c2563f42187576c025e551d186c9998ad2428c7bab',
 'native_orthogonal_wilson_charged_band':'310d8d6b5a7491fc0088e011e45dffe02efee90c9838f220f6c99a13e58e411a'}

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def enc(v,upper=True):return str(bounds.rounded(v,upper))
def interval(lo,hi):return [enc(lo,False),enc(hi)]

@lru_cache(None)
def admit():
    vals={};public={};research={}
    for name,sha in PARENTS.items():
        base=WORK/name;assert digest(base/'research_reference.json')==sha,name
        v=json.loads((base/'research_reference.json').read_text());vals[name]=v
        research[name+'/results.json']=sha
        for f,h in v['owned_source_sha256'].items():
            assert digest(base/'retained'/(f+'.raw'))==h;research[name+'/'+f]=h
        for f,h in v['research_input_sha256'].items():
            folder,filename=f.split('/',1)
            retained=WORK/folder/('research_reference.json' if filename=='results.json' else 'retained/'+filename+'.raw')
            assert digest(retained)==h
            assert f not in research or research[f]==h;research[f]=h
        for f,h in v['input_sha256'].items():
            assert digest(RER/f)==h
            assert f not in public or public[f]==h;public[f]=h
    for f in import_closure(set(public)|{'model/native_wilson_maxwell/field.py',
            'model/native_wilson_maxwell/matter.py','model/native_wilson_maxwell_refinement/bounds.py',
            'model/particle_common_contract/provenance.py'}):
        h=digest(RER/f);assert f not in public or public[f]==h;public[f]=h
    return vals,public,research

@lru_cache(None)
def helpers():
    admit()
    from ..native_orthogonal_wilson_spatial_join import model as join
    from ..native_orthogonal_wilson_charged_band import model as band
    return join,band

@lru_cache(None)
def heat_and_reference_fixture():
    vertices,edges,faces,D,C=field.cochains(2);V=len(vertices);L=D.T*D
    G=(L+s.ones(V)).inv()-s.ones(V)/V**2
    assert L*G==s.eye(V)-s.ones(V)/V and G==G.T
    assert G*s.ones(V,1)==s.zeros(V,1)
    assert max(sum(abs(x) for x in L.row(j)) for j in range(V))<=12
    # Actual Neumann tensor factorization, including endpoint degree1.
    T=s.Matrix([[1,-1,0],[-1,2,-1],[0,-1,1]])
    expected=sum((s.kronecker_product(*[T if j==i else s.eye(3) for j in range(3)])for i in range(3)),s.zeros(V))
    assert L==expected
    roots=[vertices.index(x) for x in ((0,0,0),(0,1,1),(1,1,1))]
    rows=[{'root':list(vertices[r]),'G_root_root':str(G[r,r]),
        'reference_self_energy_over_m_g_squared':str(G[r,r]/2)} for r in roots]
    assert len({row['G_root_root'] for row in rows})==3
    assert F(328,21)<16
    # Explicit conservative point-kernel budget; no huge box is executed.
    T_heat=2**20;b=2**27;N=2*b+1;volume=N**3
    assert b>=72*T_heat
    # e<3 and d!>=(d/e)^d give exp(12T)(12T)^d/d!
    # <=2^(24T-d). Avoid materializing a hundred-million-bit Fraction.
    assert T_heat.bit_length()-1+24*T_heat-b < -40
    short=F(1,2**40);root_error=short+F(T_heat,volume)+F(32,2**10)
    local_error=2*short+F(T_heat,volume)+F(64,2**10)
    assert root_error<F(1,16) and local_error<F(1,8)
    return {'literal_27_vertex_Neumann_and_Green_identity':True,'finite_reference_self_energy_fixture':rows,
        'analytic_heat_coefficient_upper':'328/21','centered_heat_bound':'16*t^(-3/2)',
        'point_Green_diagonal_upper':'33','point_Green_tail_after_T_upper':'32/sqrt(T)',
        'monopole_Green_vector_uniform_l2_claimed':False,
        'explicit_point_accuracy_example':{'heat_time':T_heat,'boundary_and_root_distance_lower':b,
            'vertices_per_axis':N,'vertices':volume,'root_local_Green_absolute_upper':enc(root_error),
            'local_Green_comparison_absolute_upper':enc(local_error),'large_box_matrices_executed':False}}

@lru_cache(None)
def actual_packet_controls():
    join,band=helpers();_,_,_,_,H,_=join.literal_box(1)
    alpha,beta,Gamma=join.source_matrices();grading=s.kronecker_product(s.eye(8),beta);A=H-grading
    assert grading*A+A*grading==s.zeros(32) and H*H==s.eye(32)+A*A
    C0=s.kronecker_product(s.eye(8),s.I*beta*alpha[1]);assert C0*H.conjugate()*C0.H==-H
    assert beta[0,0]==1 and beta[2,2]==-1
    native=band.native_algebra();H2,_,minimal,powers,_=native.cube();assert H2==H
    coeff,error=matter.sign_polynomial(minimal,512)
    diagonal=[]
    for index in (0,2):
        value=F(1,2)-sum((c*F(powers[j][index,index])/2 for j,c in enumerate(coeff)),F(0))
        diagonal.append(value)
    eps=error/2;wp=1-diagonal[0];wm=diagonal[1]
    assert wp-eps>F(11,20) and wm-eps>F(11,20)
    return {'actual_32_mode_Gamma_gap_square_and_C_identity':True,
        'beta_positive_creation_and_beta_negative_removal':True,
        'uniform_pump_success_lower':'11/20','uniform_source_free_gap_over_m_lower':'1',
        'uniform_one_body_norm_over_m_upper':'10','uniform_local_mean_density':'0 by onsite C, not total neutrality alone',
        'actual_cube_pump_success_intervals':[interval(wp-eps,wp+eps),interval(wm-eps,wm+eps)],
        'actual_cube_pump_failure_intervals':[interval(1-wp-eps,1-wp+eps),interval(1-wm-eps,1-wm+eps)],
        'sign_polynomial_terms':512,'sign_operator_error_upper':enc(error),
        'normalized_packet_tail_at_range_2Jplus1_upper':'2*100*(99/100)^(J+1)',
        'existing_complete_even_pointer_check':band.pump_controls()}

@lru_cache(None)
def wick_fixture():
    # Eight canonical modes: two vertices with all four components each.
    # This is an exact algebraic CAR countercontrol, not a native-ground claim.
    U=s.Matrix([[(-1)**((i&j).bit_count())/s.sqrt(8) for j in range(8)] for i in range(8)])
    assert U.H*U==s.eye(8)
    occupied=[0,1,2,3];P=U[:,occupied]*U[:,occupied].H
    assert s.trace(P[:4,:4])==s.trace(P[4:,4:])==2
    def state(columns):
        v=s.zeros(256,1)
        for rows in combinations(range(8),len(columns)):
            v[sum(1<<j for j in rows)]=s.det(U[list(rows),columns])
        assert (v.H*v)[0]==1;return v
    omega=state(occupied)
    charges=[[sum((word>>j)&1 for j in range(4*x,4*x+4))-2 for x in range(2)] for word in range(256)]
    def moments(v):
        means=[sum((s.conjugate(v[k])*v[k]*charges[k][x] for k in range(256)),s.S(0)) for x in range(2)]
        cov=s.Matrix(2,2,lambda x,y:sum((s.conjugate(v[k])*v[k]*charges[k][x]*charges[k][y] for k in range(256)),s.S(0)))
        return means,cov
    mean0,C= moments(omega);assert mean0==[0,0]
    G=s.Matrix([[1,-1],[-1,1]])/4;rows=[]
    for q,col,cols in ((1,4,[0,1,2,3,4]),(-1,0,[1,2,3])):
        u=U[:,col];v=state(cols);mu,Cq=moments(v)
        p=[(u[4*x:4*x+4,:].H*u[4*x:4*x+4,:])[0] for x in range(2)]
        assert mu==[q*x for x in p]
        formula=s.Matrix(2,2,lambda x,y:q*(int(x==y)*p[x]-2*s.re((u[4*x:4*x+4,:].H*P[4*x:4*x+4,4*y:4*y+4]*u[4*y:4*y+4,:])[0])))
        assert Cq-C==formula
        D=sum(G[x,x]*p[x] for x in range(2))
        exchange=sum((G[x,y]*(u[4*x:4*x+4,:].H*P[4*x:4*x+4,4*y:4*y+4]*u[4*y:4*y+4,:])[0] for x in range(2)for y in range(2)),s.S(0))
        assert 0<=exchange<=D
        for root in range(2):
            energy0=sum(G[x,y]*C[x,y] for x in range(2) for y in range(2))/2
            energy=sum(G[x,y]*(Cq[x,y]-q*int(x==root)*mu[y]-q*int(y==root)*mu[x]+int(x==root)*int(y==root))for x in range(2)for y in range(2))/2
            correction=q*(D-2*exchange)/2-sum(G[root,x]*p[x] for x in range(2))
            assert s.simplify(energy-energy0-G[root,root]/2-correction)==0
            rows.append({'q':q,'root':root,'neutral_Coulomb':str(energy0),'charged_Coulomb':str(energy),
                'root_self':str(G[root,root]/2),'intrinsic_minus_root_potential':str(correction)})
    return {'full_CAR_dimension':256,'components_per_vertex':4,'source_fixture_identified_as_actual_native_ground':False,
        'both_rank_one_Wick_and_reference_subtraction_identities':True,'all_density_pairs_checked':8,'rows':rows}

@lru_cache(None)
def finite_photon_controls():
    data=field.cube();v=data['Q_edge_vacuum_variance'];cut=field.cutoff(3);Qs=field.link_coordinates()
    vacuum=cut['states'].index((0,0,0,0,0));parity=s.diag(*[(-1)**sum(n) for n in cut['states']])
    for Q in Qs:
        assert parity*Q*parity==-Q
        assert s.simplify((Q*Q)[vacuum,vacuum]-v)==0
    r=F(1,64)*bounds.root(13);cap_error=6*bounds.coherent_word_error(3,r)
    assert cap_error<F(1,70000)
    return {'actual_cube_dimension':56,'all_twelve_literal_Q_parities_and_vacuum_second_moments':True,
        'uniform_spatial_Q_variance_upper':'13','finite_link_bound':'0<=1-<cos(gQ_K)><=g^2*v/2',
        'finite_Gaussian_or_CCR_assumption':False,'energy_cutoff_error_over_m':'6*b_K(g*sqrt(13))',
        'uniform_K3_energy_cutoff_error_over_m_upper':enc(cap_error),
        'K3_energy_error_diagnostic':float(cap_error)}

@lru_cache(None)
def calculate():
    parents,_,_=admit();join,band=helpers();family=join.family
    source=json.loads((RER/'model/native_wilson_maxwell/results.json').read_text())
    p0=F(source['actual_preparation']['complete_trial_ground_overlap_lower_per_constituent'])
    preparations=[]
    for n in (1,2,3,7,15):
        p=family.preparation(n,p0);f=F(p['complete_assembly_failure_upper'])
        p['conditional_active_state']='perpendicular free sea, then literal compensated local charged pump'
        p['original_charged_union_success_lower']=enc((1-f)*F(11,20),False)
        p['original_each_signed_charged_success_lower']=enc((1-f)*F(11,40),False)
        p['original_local_density_imbalance_lower']=enc((1-f)*F(121,1600),False)
        p['pump_failures_retained_and_not_in_assembly_failure']=True
        p['density_imbalance_dimensionless_autonomous_wait']='0'
        p['separate_positive_time_history_example_tau']='1/64'
        p['cutoff']=family.cutoff(n)
        if n<=2:
            V=p['vertices'];R=p['parallel_full_box_banks'];K=p['cutoff']['photon_cap'];d=p['cutoff']['transverse_modes']
            coefficients={2:R*4*V*256+V*R+5,5:R*4*V*8,119:R*4*V*24,
                2**31+1:R*4*V,9:R*(V-1),8*V+1:R,comb(K+d,d):R}
            capacity=(sum(k*band.log_interval(x)[0] for x,k in coefficients.items()),
                sum(k*band.log_interval(x)[1] for x,k in coefficients.items()))
            p['one_detector_full_capacity_nats']=interval(*capacity)
            p['capacity_counts_all_banks_references_flags_pointer_and_setting']=True
        preparations.append(p)
    g=F(1,64);lower=1-F(111,2)*g*g;upper=10+F(33,2)*g*g
    assert lower==F(8081,8192)>0
    return {'schema':'oph.research.orthogonal_wilson_charged_spatial_records.v1',
        'source':{'g':'1/64','a':'1/m','same_original_Psi_rho_and_normalized_Y_rule':True,
            'same_Gamma_spatial_H_and_old_raw_sea_offset':True,'reference_root_static_in_each_finite_program':True,
            'reference_root_recedes_from_fixed_observer':True,'all_six_boundaries_recede_from_observer':True,
            'positive_complement_c':'21/2; optional separate c(V) completion changes no admitted records',
            'accepted_input':'nonstationary localized charged excitation of prepared perpendicular free sea and photon vacuum',
            'coupled_ground_input_or_spectral_threshold_claimed':False},
        'packet_controls':actual_packet_controls(),'heat_reference_controls':heat_and_reference_fixture(),
        'exact_Wick_controls':wick_fixture(),'finite_photon_controls':finite_photon_controls(),
        'complete_local_law':{'root_independent_fixed_spacing_bounded_time_limit':True,
            'gauge_invariant_E_B_current_density_net_only':True,
            'charged_endpoint_is_local_physical_observable':False,
            'literal_photon_cap_and_all_assembly_pump_detector_failures_retained':True,
            'remote_fine_flag_registers_retained_but_not_identified_across_volumes':True,
            'conditional_pump_success_lower':'11/20','limiting_each_fair_signed_success_lower':'11/40',
            'conditional_density_outcome_imbalance_lower':'11/80',
            'limiting_complete_local_density_imbalance_lower':'121/1600',
            'density_instrument_effects':'(I +/- q*rho_x0/4)/2',
            'density_positive_signal_dimensionless_autonomous_wait':'0',
            'density_signal_bound_claimed_at_arbitrary_positive_wait':False,
            'local_density_contrast_claimed_as_interaction_witness':False,
            'charged_Dyson_root_phase':'-q*g^2*t*(D G e_root)_edge',
            'electric_root_phase':'-q*g*(G D^T s)_root'},
        'prepared_energy':{'baseline':'neutral prepared perpendicular free sea under the SAME H',
            'displayed_subtraction':'4+(m*g^2/2)*G(root,root)',
            'subtraction_changes_H':False,'intrinsic_Coulomb_absolute_bound_over_m':'(33/2)*g^2',
            'kinetic_photon_compression_difference_over_m_upper':'39*g^2',
            'reference_subtracted_limit_over_m':[str(lower),str(upper)],
            'lower_diagnostic':float(lower),'exact_energy_effects':'(I +/- H/Braw)/2',
            'all_energy_outcomes_and_preparation_flags_retained':True,
            'energy_compression_equals_self_energy_or_rest_mass':False,
            'raw_finite_charged_spectra_reference_independent':False},
        'finite_program_examples':preparations,
        'clock_capacity':{'duration_one_preparation_and_pump_slots':4,'total_source_time':'4+r+tau/m',
            'all_banks_continue_under_actual_H_and_all_driven_slots_cancel_background':True,
            'additional_pointer_bits_above_neutral_parent':3,'local_observer_record_alphabet':'2^(r+4)',
            'capacity_nats':'R*[4V*(256log2+8log5+24log119+log(2^31+1))+(V-1)log9+log(8V+1)+log binom(K+d,d)]+VRlog2+(r+4)log2',
            'control_norm_upper':'R*Braw_new+pi','gross_switching_work_upper':'2*(4+r)*(R*Braw_new+pi)',
            'exact_finite_matrix_services_required':True,'efficient_local_native_compilation':False},
        'scope':{'finite_program_family_not_infinite_physical_storage':True,'full_native_operator_locality':False,
            'charged_ground_or_sharp_pole_or_relativistic_mass':False,'SM_or_infinite_volume_ground_or_UV_closure':False}}


def certificate():
    v=calculate();_,public,research=admit();v['input_sha256']=public;v['research_input_sha256']=research
    v['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in('.py','.md')}
    return v

if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    value=certificate();target=HERE/'results.json'
    if args.write:target.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:assert value==json.loads(target.read_text())
    print('PASS: charged local records, receding reference, complete failures and prepared-energy limit')
