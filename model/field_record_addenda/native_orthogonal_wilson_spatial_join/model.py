"""Exact spatial/state transfer controls for the declared Gamma regulator."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
import sympy as s

HERE=Path(__file__).resolve().parent;WORK=HERE.parent
RER=HERE.parents[2]
from ..native_wilson_maxwell_spatial_family import family
bounds=family.bounds
PARENTS={
 'native_wilson_maxwell_spatial_family':'3e413171f6ad1c0df6cc4db522cb67bc4d9f5f2ae1cfc80d3f4b92fac03538ed',
 'native_wilson_maxwell_interaction_witness':'b8519b3f2464810d07b78ca90cad8f918007cdee6339f3f710ccefa4b26a8e87',
 'native_orthogonal_wilson_bridge':'75590a8e62090717d5b91ada6ac472498952c3b4ace19b15de815419ecdadc6c',
 'native_wilson_cone_audit':'5a24c444fec3ca0a03f59dde8fd52614669b9d9c730dcd4887563727c4cfa62b'}


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def admit():
    values={};public={};research={}
    for name,sha in PARENTS.items():
        base=WORK/name;assert digest(base/'research_reference.json')==sha,name
        v=json.loads((base/'research_reference.json').read_text());values[name]=v
        research[name+'/results.json']=sha
        for f,h in v['owned_source_sha256'].items():
            assert digest(base/'retained'/(f+'.raw'))==h,(name,f);research[name+'/'+f]=h
        for f,h in v.get('input_sha256',{}).items():
            assert digest(RER/f)==h,f
            assert f not in public or public[f]==h
            public[f]=h
    return values,public,research


def source_matrices():
    alpha,beta=family.matrices();gamma=alpha[0]*alpha[1]*alpha[2]*beta
    return alpha,beta,gamma


def literal_box(n,cosine=s.Rational(0),sine=s.Rational(1)):
    vertices,edges,faces,D,C=family.field.cochains(n)
    alpha,beta,gamma=source_matrices();V=len(vertices);R=cosine*beta+sine*gamma
    assert cosine*cosine+sine*sine==1
    W=D.T*D/2;H=s.MutableSparseMatrix(4*V,4*V,{})
    degree=[0]*V
    for x,y,j in edges:degree[x]+=1;degree[y]+=1
    currents=[]
    for x in range(V):H[4*x:4*x+4,4*x:4*x+4]=beta+degree[x]*R/2
    for x,y,j in edges:
        B=-(R+s.I*alpha[j])/2
        H[4*x:4*x+4,4*y:4*y+4]=B;H[4*y:4*y+4,4*x:4*x+4]=B.H
        J=s.MutableSparseMatrix(4*V,4*V,{})
        J[4*x:4*x+4,4*y:4*y+4]=-s.I*B;J[4*y:4*y+4,4*x:4*x+4]=s.I*B.H
        currents.append(J)
    grading=s.SparseMatrix(s.kronecker_product(s.eye(V),beta));W4=s.SparseMatrix(s.kronecker_product(W,s.eye(4)))
    A=H-grading
    assert grading*A+A*grading==2*cosine*W4
    assert H*H-s.eye(4*V)==2*cosine*W4+A*A
    C0=s.I*beta*alpha[1];charge=s.SparseMatrix(s.kronecker_product(s.eye(V),C0))
    assert charge*H.conjugate()*charge.H==-H
    assert all(charge*J.conjugate()*charge.H==J and s.trace(J)==0 for J in currents)
    return vertices,edges,D,C,H,currents


def finite_controls():
    vertices,edges,D,C,H,currents=literal_box(1)
    alpha,beta,gamma=source_matrices();V=len(vertices)
    sy=s.Matrix([[0,-s.I],[s.I,0]]);T0=s.diag(s.I*sy,s.I*sy)
    inv=s.zeros(V);lookup={tuple(x):j for j,x in enumerate(vertices)}
    for j,x in enumerate(vertices):inv[lookup[tuple(1-k for k in x)],j]=1
    P=s.kronecker_product(inv,beta);T=s.kronecker_product(s.eye(V),T0)
    PT=P*T
    assert P*H*P.H!=H and T*H.conjugate()*T.H!=H
    assert PT*H.conjugate()*PT.H==H and PT*PT.conjugate()==-s.eye(4*V)
    literal_box(1,s.Rational(3,5),s.Rational(4,5))
    # Noncommutation at an open boundary is real, rather than suppressed
    # by the positive-square proof.
    W=s.kronecker_product(D.T*D/2,s.eye(4))
    K=H-s.kronecker_product(s.eye(V),beta)-s.kronecker_product(D.T*D/2,gamma)
    assert K*W-W*K!=s.zeros(4*V)
    return {'actual_open_cube_homotopy_square_identity':True,
        'all_12_actual_current_charge_conjugation_identities':True,
        'derivative_and_Wilson_Laplacian_commute_at_boundary':False,
        'ordinary_parity_preserved':False,'ordinary_time_reversal_preserved':False,
        'free_open_box_combined_PT_preserved':True,'free_PT_square':'-I',
        'combined_PT_extended_to_fixed_charged_reference_or_native_complement':False}


def bulk_controls():
    A,c,b=s.symbols('A c b',real=True);E2=A*A+2+b-2*A*c
    assert s.expand((A*c-1)*E2+(A*c-1)**2+1+b-A*(c*E2-A*(1-c*c)))==0
    w=s.symbols('w0:3',nonnegative=True);W=sum(w);cross=sum(w[i]*w[j] for i in range(3) for j in range(i+1,3))
    new=1+2*W+2*cross;old=1+4*W+2*cross
    assert s.expand(old-new-2*W)==0
    assert new.subs(dict.fromkeys(w,2))==37
    r=s.Symbol('r',real=True);delta=3*(1-r)
    assert s.expand(r*r+delta*r-(3*r-2*r*r))==0
    assert s.expand(3*r-2*r*r-1+(r-1)*(2*r-1))==0
    alpha,beta,gamma=source_matrices();C0=s.I*beta*alpha[1]
    sy=s.Matrix([[0,-s.I],[s.I,0]]);T0=s.diag(s.I*sy,s.I*sy);PT=beta*T0
    assert C0*gamma.conjugate()*C0.H==-gamma
    assert PT*gamma.conjugate()*PT.H==gamma
    return {'fresh_IBP_identity':True,'bare_hopping_magnitude_lower':'2/343',
        'bulk_free_energy_squared_upper':'37','raw_sea_shift_density_lower_over_m':'6/7',
        'photon_vacuum_kinetic_compression':'beta+r alpha.sin(q)+Gamma*(3-r sumcos(q))',
        'compressed_Gamma_rest_term':'3(1-r)','compressed_squared_curvature':'3r-2r^2 >1 for 1/2<r<1',
        'Gamma_mass_term_allowed_by_proper_rotations_C_and_PT':True,
        'compression_identified_with_full_self_energy_or_particle_mass':False}


def scientific():
    parents,_,_=admit();g=F(1,64);u=F(1,2**16);tau=F(1,2**32)
    dressed=(F(2,343)-400*F(99,100)**2049)*(1-F(13,8192))
    assert dressed>F(1,173)
    contrast=u/2*(g/F(173)-F(64,3)*u*u-1024*tau)
    old=parents['native_wilson_maxwell_interaction_witness']
    example=json.loads(json.dumps(old['finite_sufficient_example']))
    example['preparation']['conditional_active_state']='exact perpendicular free sea in the fixed original Psi/rho/Y frame, with photon vacuum'
    example['preparation']['additional_physical_spectral_transport_in_routing_slot']=True
    f=F(example['preparation']['complete_assembly_failure_upper'])
    d=F(old['bounds_exact']['complete_finite_cap_record_error_upper'])
    separation=(1-f)*contrast/2-2*d
    # A source-only quantitative free-volume obstruction. No enormous
    # graph, matrix or native state is enumerated.
    J=2048;N=2**24;b=F(6*(2*J+2),N)
    density=F(6,7)-8000*F(99,100)**(J+1)-80*b
    assert density>F(1,2) and density*N**3>F(21,2)+1
    bridge=parents['native_orthogonal_wilson_bridge']['full_physical_ground']
    assert F(bridge['native_physical_complement']['actual_neutral_triplet_gap_interval'][1])<F(bridge['source_selected_mass_interval'][0])
    r_lower=1-F(13,8192)
    x_lower=g*g/(6*bounds.root(12))
    s_lower=x_lower/(1+x_lower)
    curvature_lower=s_lower*(1-2*F(13,8192))
    assert r_lower>F(1,2) and curvature_lower>0
    _,beta,gamma=source_matrices();atomic=(s.eye(4)-beta)/2
    assert atomic*atomic==atomic and s.trace(atomic)==2
    assert s.trace(atomic*beta)==-2 and s.trace(atomic*gamma)==0
    vals={'finite_dressed_hopping_magnitude_lower':dressed,
        'matched_setting_probability_contrast_lower':contrast,
        'complete_finite_law_TV_separation_lower':separation,
        'limiting_complete_law_TV_separation_lower':contrast/2,
        'explicit_free_large_box_raw_sea_shift_density_lower':density,
        'compressed_squared_curvature_minus_one_lower':curvature_lower}
    exact={k:str(bounds.rounded(v,False)) for k,v in vals.items()}
    return {'schema':'oph.research.orthogonal_wilson_spatial_join.v1',
        'finite_controls':finite_controls(),'bulk_controls':bulk_controls(),
        'bounds_exact':exact,'diagnostic_binary64':{k:float(F(v)) for k,v in exact.items()},
        'declared_family':{'m_and_charge_frame':'same native spectral m, original Psi/rho/Y',
            'a':'1/m','g':'1/64','positive_completion_c':'21/2',
            'old_beta_sea_raw_scalar_retained':True,'new_perp_vacuum_subtraction':False,
            'source_state':'new perpendicular free sea after finite endpoint spectral transport in fixed charge frame',
            'all_native_complements_references_and_failure_flags_retained':True},
        'homotopy':{'positive_square':'h_theta^2/m^2=I+2 cos(theta)W+A_theta^2',
            'one_body_gap_lower_over_m':'1','projector_derivative_norm_upper':'3',
            'Kato_transport_is_mathematical_rank_and_gap_argument_only':True,
            'chosen_preparation_service':'lexicographic algebraic bases of both endpoint projector ranges; U_alg maps matching bases',
            'physical_homotopy_clock_or_continuous_carrier_added':False,
            'free_sea_rotation_in_same_routing_slot':True},
        'same_fixed_spacing_complete_local_record_limit':True,
        'finite_witness':example,
        'free_complement_countercontrol':{'vertices_per_axis':N,'polynomial_order':J,
            'fixed_c21_over2_selects_code_ground_uniformly':False,
            'native_complement_lower_at_g_zero_for_this_sufficient_box':True,
            'nonzero_g_ground_assignment_inferred':False},
        'separate_source_only_completion_option':{
            'c_N':'18V+1','slope':18,
            'variational_state':'exact beta-negative two-particle onsite occupations, photon vacuum, zero rho and reference labels',
            'variational_code_trial_energy_over_raw_baseline':'m*(S_beta,n-2V)',
            'variational_code_trial_energy_over_raw_baseline_upper':'18mV',
            'complement_separation_from_code_ground_lower':'m',
            'all_caps_K_at_least_3':True,'any_coupling_g':True,'retained_raw_offset':True,
            'admitted_code_record_laws_unchanged':True,
            'raw_norm_addition_upper_for_this_separate_option':'34Vm+Vm c(N)',
            'actual_ground_vector_or_internal_gap_or_volume_limit_proved':False},
        'clock_and_resource_changes':{'complete_elapsed_time':'5+tau/m',
            'raw_norm_addition_upper':'(89/2) V m',
            'raw_norm_addition_applies_to':'fixed c=21/2 family',
            'all_source_reference_pointer_capacity_unchanged_from_finite_witness':True,
            'spectral_transport_coefficient_workspace_included':False},
        'boundaries':{'common_free_cone_is_interacting_particle_dispersion':False,
            'fixed_cube_ground_theorem_inherited_in_volume':False,
            'relativistic_or_UV_limit_or_SM_mass':False}}


def certificate():
    value=scientific();_,public,research=admit()
    value['input_sha256']=public;value['research_input_sha256']=research
    value['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    return value


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();value=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:assert value==json.loads(path.read_text())
    print('PASS: fresh Gamma family, same charge/Gauss frame, controlled interacting local records; ground limits explicit')
