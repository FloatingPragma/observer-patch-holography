"""Exact charged-candidate discriminators on the primary Wilson source."""
from functools import lru_cache
from fractions import Fraction as F
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_channel_labels import channels
from model.particle_primary_regional_fermion import model as regional
from model.particle_primary_odd_sum_rules import model as parity
from model.particle_primary_selected_odd_response import model as selected_odd
from model.particle_spatial_modular_response import regional as support
from model.native_quantum_gauss import algebra as native
from model.particle_local_wilson_density import density

SOURCE_INPUTS=(
    'model/particle_candidate_charged/__init__.py',
    'model/particle_candidate_charged/__main__.py',
    'model/particle_candidate_charged/model.py',
    'model/particle_candidate_charged/test_model.py',
    'model/particle_candidate_charged/PROOF.md',
    'model/particle_candidate_charged/README.md',
    'Lean/QFT/GaugeIrreducibleBorn.lean',
    'Lean/Screen/ElectroweakBreakingComposition.lean',
    'model/standard_model_source_map/KINETIC_ATTACHMENT.md',
    *regional.SOURCE_INPUTS,*parity.SOURCE_INPUTS,*channels.SOURCE_INPUTS,*selected_odd.SOURCE_INPUTS,
)

NAMES=('u_conjugate','lepton_doublet','quark_doublet','d_conjugate','e_conjugate')
CHARGES=(-4,-3,1,2,6)
DIMS=(3,2,6,3,1)


def clean(A):
    return s.SparseMatrix(A).applyfunc(s.simplify)


def blocks():
    result=tuple(sorted((int(-2*s.Rational(q)),tuple(idx)) for q,idx in native.native_blocks()))
    if tuple((q,len(idx)) for q,idx in result)!=tuple(zip(CHARGES,DIMS)):
        raise ArithmeticError('the actual native representation changed')
    return result


def commutant_system(matrices):
    """Row-major vectorization of X T - T X, with exact algebraic entries."""
    n=matrices[0].rows
    return s.Matrix.vstack(*(s.kronecker_product(s.eye(n),T.T)-s.kronecker_product(T,s.eye(n)) for T in matrices))


@lru_cache(maxsize=1)
def coefficient_control():
    currents=native.native_currents();rows=[]
    center=clean(-2*currents[-1]/s.I)
    for name,(q,idx) in zip(NAMES,blocks()):
        matrices=tuple(T.extract(idx,idx) for T in currents)
        system=commutant_system(matrices);kernel=system.nullspace();d=len(idx)
        if len(kernel)!=1 or clean(kernel[0].reshape(d,d))!=s.eye(d):
            raise ArithmeticError('a native block has a non-scalar coefficient commutant')
        P=s.eye(15)
        for other in CHARGES:
            if other!=q:P=P*(center-other*s.eye(15))/(q-other)
        if clean(P)!=s.diag(*(int(i in idx) for i in range(15))):
            raise ArithmeticError('central spectral projectors changed')
        rows.append({'native_name':name,'Q6':q,'internal_dimension':d,
            'exact_constraint_rank':d*d-1,'commutant_nullity':1,
            'native_mask_rows':list(idx),'spinor_fiber_dimension':4*d})
    g=density.gamma_data();beta=g[3]
    # beta and both orientations of each literal hopping coefficient generate gamma_i.
    coefficients=(beta,)+tuple((-beta+beta*x)/2 for x in g[:3])+tuple((-beta+beta*x).H/2 for x in g[:3])
    spin_system=commutant_system(coefficients);spin_kernel=spin_system.nullspace()
    if len(spin_kernel)!=1 or spin_kernel[0].reshape(4,4)!=s.eye(4):
        raise ArithmeticError('the literal Wilson spin coefficients have a larger commutant')
    return {'native_blocks':rows,'native_gauge_generators':len(currents),
        'spin_coefficient_constraint_rank':15,'spin_coefficient_commutant_nullity':1,
        'gauge_only_local_commutant_dimension':80,
        'gauge_and_each_Wilson_coefficient_commutant_dimension':5,
        'generated_local_complex_algebra_dimension':sum((4*d)**2 for d in DIMS),
        'full_H_commutant_computed':False,'onsite_three_flavor_tensor_factor_present':False,
        'composite_or_spatial_generations_excluded':False,
        'alternate_declared_three_flavor_extension_is_same_source':False}


@lru_cache(maxsize=1)
def mass_control():
    g=density.gamma_data();beta=g[3];alpha=tuple(s.I*beta*x for x in g[:3])
    gamma5=clean(g[0]*g[1]*g[2]*beta)
    anti=s.Matrix.vstack(*(s.kronecker_product(s.eye(4),x.T)+s.kronecker_product(x,s.eye(4)) for x in alpha))
    even=s.Matrix.vstack(anti,commutant_system((beta,)))
    k=anti.nullspace();ke=even.nullspace()
    if len(k)!=2 or len(ke)!=1:raise ArithmeticError('the Clifford mass coefficient space differs')
    for M in (beta,s.I*beta*gamma5):
        if M.H!=M or any(clean(M*a+a*M)!=s.zeros(4) for a in alpha):
            raise ArithmeticError('a proposed Hermitian Clifford mass coefficient is invalid')
    if s.Matrix.hstack(ke[0],s.Matrix(beta).reshape(16,1)).rank()!=1:
        raise ArithmeticError('parity-even mass space is not beta')
    if clean(beta*gamma5*beta)!=-gamma5:raise ArithmeticError('native chirality parity changed')
    return {'kinetic_matrices':'alpha_i=i beta gamma_i',
        'unrestricted_exact_constraint_rank':14,'unrestricted_mass_space_dimension':2,
        'Hermitian_mass_basis':['beta','i beta gamma5'],
        'parity_even_exact_constraint_rank':15,'parity_even_mass_space_dimension':1,
        'parity_even_gauge_preserving_mass_space':'beta tensor sum_r m_r P_r',
        'real_mass_coefficients_allowed_by_local_symmetries':5,
        'actual_onsite_coefficient':'3 beta/a',
        'onsite_coefficient_is_a_rest_mass':False,
        'actual_native_chirality_commutes_with_H':False,
        'native_chirality_is_physical_Lorentz_chirality':False,
        'bare_cross_irrep_quadratic_mass_allowed':False}


def ward_control():
    e=tuple(tuple(int(i==j) for i in range(5)) for j in range(5))
    pairs=(('charged_lepton',1,4,-1),('down_quark',2,3,-1),('up_quark',2,0,1))
    rows=[]
    for name,a,b,sign in pairs:
        delta=tuple(sign*(e[a][i]+e[b][i]) for i in range(5))
        rows.append({'channel':name,'bare_factors':[NAMES[a],NAMES[b]],
            'required_H_dagger_native_increment':list(delta),
            'required_H_dagger_Q6':sum(q*n for q,n in zip(CHARGES,delta)),
            'bare_pair_native_increment':list(tuple(e[a][i]+e[b][i] for i in range(5))),
            'bare_pair_scalar_anomalous_kernel_zero':True})
    if len({tuple(x['required_H_dagger_native_increment']) for x in rows})!=3 or any(x['required_H_dagger_Q6']!=-3 for x in rows):
        raise ArithmeticError('the common-Higgs native Ward comparison changed')
    return {'native_order':list(NAMES),'conserved_native_occupations':5,
        'normal_kernel_selection':'Delta_right=Delta_left',
        'anomalous_kernel_selection':'Delta_left+Delta_right=0',
        'kernel_domain':'scalar source-marginal correlations with neutral H functions',
        'charge_coherent_reference_or_compensator_excluded_from_domain':True,
        'Wilson_link_alone_cancels_native_occupation':False,
        'bare_Yukawa_channels':rows,'one_sharp_native_number_H_dagger_for_all_three':False,
        'same_hypercharge_implies_same_native_sector':False,
        'dressed_fields_or_selected_broken_states_excluded':False}


def second_quantize(T,creators):
    out=s.zeros(creators[0].rows)
    for i,j in product(range(T.rows),repeat=2):
        if T[i,j]:out+=T[i,j]*creators[i]*creators[j].H
    return clean(out)


@lru_cache(maxsize=1)
def dressing_control():
    """Full six-mode CAR algebra, including all partially occupied inputs."""
    C=tuple(regional.creator(j,6) for j in range(6));a=tuple(x.H for x in C)
    phi=tuple(clean((a[2*j]*a[5]-a[2*j+1]*a[4])/s.sqrt(2)) for j in range(2))
    D=tuple(clean(sum((C[2*j+t]*phi[j] for j in range(2)),s.zeros(64))) for t in range(2))
    N=tuple(C[j]*a[j] for j in range(6));NL=sum(N[:4],s.zeros(64));NE=sum(N[4:],s.zeros(64))
    eps=s.Matrix([[0,1],[-1,0]]);anchor=s.eye(64)[:,63]
    Lidx=blocks()[1][1];eidx=blocks()[4][1]
    for T in native.native_currents():
        TL=T.extract(Lidx,Lidx);te=T[eidx[0],eidx[0]]
        J=second_quantize(s.diag(s.kronecker_product(TL,s.eye(2)),te*s.eye(2)),C)
        for j in range(2):
            expected=-sum((TL[j,k]*phi[k] for k in range(2)),s.zeros(64))-te*phi[j]
            if clean(J*phi[j]-phi[j]*J)!=clean(expected):raise ArithmeticError('composite scalar gauge covariance failed')
            if clean(J*D[j]-D[j]*J)!=clean(-te*D[j]):raise ArithmeticError('dressed field gauge covariance failed')
    sigma=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    for K in (s.I*x/2 for x in sigma):
        J=second_quantize(s.diag(s.kronecker_product(s.eye(2),K),K),C)
        for j in range(2):
            if clean(J*phi[j]-phi[j]*J)!=s.zeros(64):raise ArithmeticError('the composite scalar carries source spin')
            expected=sum((K[t,j]*D[t] for t in range(2)),s.zeros(64))
            if clean(J*D[j]-D[j]*J)!=clean(expected):raise ArithmeticError('the dressed odd field is not a source spinor')
    for g in (row['endpoint'] for row in regional.spin_control()['rows']):
        if clean(g.T*eps*g)!=eps:raise ArithmeticError('the actual binary cube lift does not preserve epsilon')
    for j in range(2):
        if clean(NL*phi[j]-phi[j]*NL)!=-phi[j] or clean(NE*phi[j]-phi[j]*NE)!=-phi[j]:
            raise ArithmeticError('the scalar native increments differ')
        if clean(NL*D[j]-D[j]*NL)!=s.zeros(64) or clean(NE*D[j]-D[j]*NE)!=-D[j]:
            raise ArithmeticError('the dressed odd native increments differ')
        target=s.sqrt(2)*sum((eps[j,u]*a[4+u]*anchor for u in range(2)),s.zeros(64,1))
        if clean(D[j]*anchor)!=clean(target):raise ArithmeticError('the full-occupation dressing identity failed')
    effect=clean(sum((x.H*x for x in D),s.zeros(64)))
    ev=effect.eigenvals();phi_ev=tuple(clean(x.H*x).eigenvals() for x in phi);D_ev=tuple(clean(x.H*x).eigenvals() for x in D)
    if max(ev)!=s.Rational(9,2) or any(max(x)!=1 for x in phi_ev) or any(max(x)!=3 for x in D_ev):
        raise ArithmeticError('the exact full-CAR norm constants changed')
    P=s.diag(*((-1)**mask.bit_count() for mask in range(64)))
    if any(P*x!=x*P for x in phi) or any(P*x!=-x*P for x in D):raise ArithmeticError('graded parity changed')
    return {'creators':C,'phi':phi,'D':D,'effect':effect,'anchor':anchor,'NL':NL,'NE':NE,
        'effect_spectrum':{str(k):int(v) for k,v in sorted(ev.items())},
        'summary':{'full_local_CAR_modes':6,'full_local_Fock_dimension':64,
            'actual_native_generators_checked':12,'actual_binary_cube_lifts_checked':48,
            'scalar_native_increment':[0,-1,0,0,-1],'scalar_Q6':-3,
            'scalar_internal_representation':'(1, anti-2) with Y=-1/2',
            'scalar_source_spin_irrep':'A1','scalar_CAR_parity':0,'scalar_norm':1,
            'scalar_original_one_point':0,'scalar_selected_Higgs_vacuum_proved':False,
            'zero_charged_one_point_excludes_gauge_invariant_Higgs_phase':False,
            'dressed_native_increment':[0,0,0,0,-1],'dressed_Q6':-6,
            'dressed_internal_representation':'(1,1) with Y=-1',
            'dressed_source_spin_irrep':'G1','dressed_CAR_parity':1,
            'each_dressed_squared_norm':3,'summed_dressed_effect_norm':F(9,2),
            'anchor_overlap_with_bare_e_conjugate_hole_squared':2,
            'D_equals_constant_times_bare_hole_on_full_Fock':False,
            'physical_Higgs_or_charged_lepton_identified':False}}


@lru_cache(maxsize=1)
def compensated_control():
    """Four source-neutral rows on the complete eight-mode CAR/link control."""
    d=dressing_control();C=tuple(regional.creator(j,8) for j in range(8))
    D=tuple(s.SparseMatrix(s.kronecker_product(s.eye(4),x)) for x in d['D'])
    rows=tuple(clean(C[p]*D[t]/3) for p,t in product((6,7),range(2)))
    Np=C[6]*C[6].H+C[7]*C[7].H
    ED=s.SparseMatrix(s.kronecker_product(s.eye(4),d['effect']));E=clean(sum((B.H*B for B in rows),s.SparseMatrix(256,256,{})))
    if E!=clean((2*s.eye(256)-Np)*ED/9):raise ArithmeticError('compensated full-CAR effect differs')
    # The two factors act on disjoint even algebras, hence their exact joint spectrum is a product.
    spectrum={s.Rational(n)*s.Rational(ev)/9 for n in (0,1,2) for ev in d['effect_spectrum']}
    if min(spectrum)!=0 or max(spectrum)!=1:raise ArithmeticError('the complete compensated menu is not a contraction')
    Qy=s.SparseMatrix(s.kronecker_product(s.eye(4),-3*d['NL']+6*d['NE']))
    for B in rows:
        if clean(6*Np*B-B*6*Np)!=6*B or clean(Qy*B-B*Qy)!=-6*B:
            raise ArithmeticError('endpoint charges fail to match the Wilson character')
    anchor=s.eye(256)[:,63];targets=s.Matrix.hstack(*(B*anchor for B in rows))
    if clean(targets.H*targets)!=s.Rational(2,9)*s.eye(4):raise ArithmeticError('physical-anchor row weights changed')
    # chi P0 has link shift +1 and initial constant projection. These commutators close both endpoints.
    link=s.zeros(3);link[2,1]=1;electric=s.diag(-1,0,1);P0=s.diag(0,1,0)
    if electric*link-link*electric!=link or link.H*link!=P0:
        raise ArithmeticError('the admitted charge-six character shift differs')
    return {'rows':rows,'effect':E,'anchor':anchor,'Qy':Qy,'Np':Np,
        'summary':{'inside_CAR_modes':6,'outside_CAR_modes':2,'complete_local_Fock_dimension':256,
            'link_control_dimension':3,'gauge_complete_rows':4,'complete_probe_outcomes':5,
            'row_definition':'K_ps=chi_xy P0_xy c_x,p^dagger D_y,s /3',
            'failure':'sqrt(I-sum_ps K_ps^dagger K_ps)',
            'four_row_effect_norm':1,'each_anchor_success_weight':F(2,9),
            'total_anchor_success_weight':F(8,9),
            'inside_Q6_increment':-6,'outside_Q6_increment':6,'total_native_increment':[0,0,0,0,0],
            'complete_pair_CAR_parity':0,'individual_regional_factor_CAR_parity':1,
            'inside_factor_alone_preserves_physical_Gauss_carrier':False,
            'opposite_charge_and_reference_retained':True,'source_action_extended':False,
            'source_spin_and_gauge_labels_are_probe_labels':True,
            'physical_single_particle_energy_or_mass_computed':False}}


@lru_cache(maxsize=1)
def coherent_overlap_control():
    """Positive original bare/dressed comparison, without charged intermediate projection."""
    d=dressing_control();comp=compensated_control();C=d['creators']
    h=(C[5].H,-C[4].H)
    scalar_effect=clean(sum((x.H*x for x in d['phi']),s.zeros(64)))
    mixed=clean(sum((h[j].H*d['D'][j] for j in range(2)),s.zeros(64)))
    if mixed!=clean(s.sqrt(2)*scalar_effect):raise ArithmeticError('the scalar mixed-overlap identity failed')
    C8=tuple(regional.creator(j,8) for j in range(8))
    h8=tuple(s.SparseMatrix(s.kronecker_product(s.eye(4),x)) for x in h)
    bare=tuple(clean(C8[p]*h8[t]/2) for p,t in product((6,7),range(2)))
    bare_effect=clean(sum((A.H*A for A in bare),s.zeros(256)))
    NE=s.SparseMatrix(s.kronecker_product(s.eye(4),d['NE']))
    if bare_effect!=clean((2*s.eye(256)-comp['Np'])*NE/4):raise ArithmeticError('the bare comparison menu is not complete')
    mixed_pair=clean(sum((A.H*B for A,B in zip(bare,comp['rows'])),s.zeros(256)))
    phi8=s.SparseMatrix(s.kronecker_product(s.eye(4),scalar_effect))
    positive=clean(s.sqrt(2)*(2*s.eye(256)-comp['Np'])*phi8/6)
    if mixed_pair!=positive or clean(comp['anchor'].H*mixed_pair*comp['anchor'])[0]!=2*s.sqrt(2)/3:
        raise ArithmeticError('the compensated original overlap differs')
    plus=tuple(clean((A+B)/2) for A,B in zip(bare,comp['rows']))
    minus=tuple(clean((A-B)/2) for A,B in zip(bare,comp['rows']))
    total=clean(sum((L.H*L for L in (*plus,*minus)),s.zeros(256)))
    difference=clean(sum((L.H*L for L in plus),s.zeros(256))-sum((L.H*L for L in minus),s.zeros(256)))
    if total!=clean((bare_effect+comp['effect'])/2) or difference!=positive:
        raise ArithmeticError('the coherent comparison complete effect or contrast failed')
    onsite=clean(3*(comp['Np']-s.kronecker_product(s.eye(4),d['NL']+d['NE'])))
    if any(clean(onsite*B-B*onsite)!=6*B for B in comp['rows']):
        raise ArithmeticError('the dressed onsite transfer is not six')
    return {'bare_rows':bare,'plus_rows':plus,'minus_rows':minus,'effect':total,
        'positive_difference':positive,'local_scalar_effect':scalar_effect,
        'summary':{'local_gauge_invariant_mixed_identity':'sum_s h_s^dagger D_s=sqrt(2) sum_a Phi_a^dagger Phi_a',
            'charged_intermediate_Gauss_projection_inserted':False,
            'positive_compensated_mixed_identity':'sum_ps A_ps^dagger B_ps=(sqrt(2)/6) P0 (2-Np) sum_a Phi_a^dagger Phi_a',
            'normalized_bare_rows':'A_ps=chi P0 c_p^dagger h_s/2',
            'normalized_dressed_rows':'B_ps=chi P0 c_p^dagger D_s/3',
            'complete_success_rows':'L_plusminus_ps=(A_ps +/- B_ps)/2',
            'complete_probe_outcomes':9,
            'failure':'sqrt(I-sum_plusminus_ps L^dagger L)',
            'positive_original_contrast_multiplier_of_density_floor':'2*sqrt(2)/3',
            'sum_plus_probability_minus_sum_minus_probability_strictly_positive':True,
            'nonzero_common_compensated_energy_transition_support_proved':True,
            'common_transition_energy_location_computed':False,
            'common_transition_frequency_proved_nonzero':False,
            'exact_onsite_and_electric_transfer_center':'87/(4a)',
            'hopping_commutator_included_in_that_center':False,
            'center_is_a_rest_mass':False,
            'coherent_source_reference_records_retained':True}}


@lru_cache(maxsize=1)
def weak_hole_residual_control():
    """Exact positive block inequality on every local occupation sector."""
    d=dressing_control();h=(d['creators'][5].H,-d['creators'][4].H)
    residual=tuple(clean(d['D'][j]-s.sqrt(2)*h[j]) for j in range(2))
    G=clean(sum((E.H*E for E in residual),s.zeros(64)))
    rows=[]
    for holes in range(5):
        indices=[bits for bits in range(64) if 4-(bits&15).bit_count()==holes]
        eigen=G.extract(indices,indices).eigenvals()
        if max(eigen)>s.Rational(3,2)*holes:raise ArithmeticError('the complete weak-hole residual bound failed')
        rows.append({'weak_holes':holes,'full_block_dimension':len(indices),
            'residual_effect_spectrum':{str(k):int(v) for k,v in sorted(eigen.items())},
            'largest_eigenvalue':str(max(eigen))})
    hole=4*s.eye(64)-d['NL']
    if clean(G*hole-hole*G)!=s.zeros(64):raise ArithmeticError('residual is not block diagonal in weak holes')
    if [x['largest_eigenvalue'] for x in rows]!=['0','1','3','3','4']:
        raise ArithmeticError('the full occupation residual spectrum changed')
    return {'residuals':residual,'effect':G,'weak_holes':hole,'blocks':rows,
        'optimal_hole_coefficient':F(3,2)}


def dressed_energy_window(member,regional_bound=None):
    """Transfer a recorded event by a norm bound on complete original histories."""
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('the directed-cube dressed window is member zero only')
    rbound=regional.bind_member(member) if regional_bound is None else regional_bound
    regional.verify_bound(member,rbound)
    obound=selected_odd.bind_member(member);selected_odd.verify_bound(member,obound)
    weak_hole_residual_control()
    raw=rbound['result']['raw_spin_apparatus'];p=F(raw['per_fixed_spin_entry_recorded_weight_lower'])
    if p!=4*F(raw['each_of_four_joint_spin_rows_recorded_weight_lower']):
        raise ArithmeticError('the bare event is not the summed four-spin complete menu')
    defect=F(obound['result']['weak_plus_e_defect_count_upper']);epsilon=defect/48
    if p<=F(43,160) or epsilon>=F(129,400)**2:
        raise ArithmeticError('the original selected response margins do not transfer')
    if F(61,125)**2>=F(8,9)*F(43,160):raise ArithmeticError('amplitude lower is not outward')
    lower=(F(61,125)-F(129,400))**2
    if lower<=F(1,40):raise ArithmeticError('the dressed original finite event is not positive')
    edge_rows=regional.directed_menu()
    edges={(tuple(row['target']),tuple(row['compensator'])) for row in edge_rows}
    if len(edges)!=24 or any(sum(y==origin for x,y in edges)!=3 for origin in product(range(2),repeat=3)):
        raise ArithmeticError('the complete member-zero directed edge allocation changed')
    return contract.serial({'same_original_source_state_reference_clock':True,
        'regional_bare_result_sha256':rbound['result_sha256'],
        'selected_defect_result_sha256':obound['result_sha256'],
        'complete_success_outcomes':96,'complete_probe_outcomes':97,
        'directed_edges':24,'spin_rows_per_edge':4,
        'bare_rows':'A_xyps=chi P0 c_p^dagger h_s/(2 sqrt(24))',
        'dressed_rows':'B_xyps=chi P0 c_p^dagger D_s/(3 sqrt(24))',
        'comparison_multiplier':'2 sqrt(2)/3','comparison_multiplier_squared':F(8,9),
        'local_residual_effect_bound':'sum_s (D_s-sqrt(2)h_s)^dagger(D_s-sqrt(2)h_s) <= (3/2)(4-N_L)',
        'exact_local_weak_hole_blocks':weak_hole_residual_control()['blocks'],
        'original_weak_plus_e_defect_count_upper':defect,
        'original_negative_weak_hole_count_upper':defect/2,
        'complete_history_residual_squared_upper':epsilon,
        'original_bare_total_menu_recorded_weight_lower':p,
        'bare_spin_rows_summed':4,
        'original_bare_rational_weight_strict_lower':F(43,160),
        'scaled_bare_amplitude_strict_lower':F(61,125),
        'residual_amplitude_strict_upper':F(129,400),
        'original_dressed_recorded_weight_strict_lower':lower,
        'original_dressed_recorded_simple_strict_lower':F(1,40),
        'reported_transfer_interval':raw['reported_transfer_interval'],
        'complete_clock_plan':rbound['result']['complete_clock_plan'],
        'additional_clock_loss_subtracted':0,
        'full_purified_history_triangle_used':True,
        'different_failure_outcomes_retained':True,
        'initial_clock_conditions_original_state':False,
        'energy_window_is_a_single_endpoint_or_rest_mass':False,
        'physical_charged_lepton_identified':False})


def bind_member(member):
    contract.verify_member(member);channels.native_number_control()
    chirality=parity.bind_chirality_member(member);parity.verify_chirality_bound(member,chirality)
    coefficient=coefficient_control();mass=mass_control();ward=ward_control();dress=dressing_control();comp=compensated_control();overlap=coherent_overlap_control()
    floor=support.common_contract_payload(member,1)['physical_density_floor_factored']
    primary=None;window=None
    if member['source']['level']==0:
        bound=regional.bind_member(member);regional.verify_bound(member,bound)
        primary={'result_sha256':bound['result_sha256'],
            'raw_spin_per_joint_row_original_energy_weight_simple_lower':'1/16',
            'reported_pair_transfer_interval':bound['result']['raw_spin_apparatus']['reported_transfer_interval'],
            'pair_interval_is_single_lepton_mass':False}
        window=dressed_energy_window(member,bound)
    records=[]
    for name,q,d in zip(NAMES,CHARGES,DIMS):
        records.append({'native_name':name,'Q6':q,'internal_dimension':d,
            'supported':'exact native gauge representation and two source beta spinor copies',
            'source_spin':'G1 + G1 before a beta projection',
            'generation_discriminator_in_local_coefficient_commutant':None,
            'single_handed_source_assignment':'excluded for the stated equivariant native-parity/current dictionary on original complete records',
            'physical_Lorentz_spin':None,'physical_generation':None,'mass_GeV':None,
            'identified_physical_particle':None})
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'native_coefficient_control':coefficient,'Clifford_mass_control':mass,
        'native_number_Ward_control':ward,'scalar_and_dressed_odd_control':dress['summary'],
        'compensated_dressed_instrument':comp['summary'],
        'original_coherent_bare_dressed_join':overlap['summary'],
        'original_physical_density_floor':floor,
        'dressed_success_original_weight_lower_factored':{'prefactor':F(8,9),'base':2,'negative_exponent':floor['negative_exponent']},
        'each_dressed_row_original_weight_lower_factored':{'prefactor':F(2,9),'base':2,'negative_exponent':floor['negative_exponent']},
        'all_complete_source_clock_records_retained':True,
        'positive_success_bound_is_uniform_in_capacity':False,
        'dressed_positive_narrow_energy_band_derived':False,
        'dressed_positive_reported_energy_window_derived':window is not None,
        'original_reference_output_Loewner_inequality_claimed':False,
        'source_state_reference_clock_changed':False,
        'native_parity_record_result_sha256':chirality['result_sha256'],
        'primary_member_zero_regional_energy_attachment':primary,
        'original_compensated_dressed_energy_window':window,
        'candidate_discriminators':records,'physical_species_confirmed':[],
        'physical_mass_table_GeV':[],
        'next_dependent_calculation':'resolve an energy/momentum branch or charged threshold, separate endpoint and flux response, and prove capacity control on the same native currents and parity dictionary'})
    return {'schema':'oph.particle_candidate_charged_bound.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('charged candidate source or bound differs')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    result={'schema':'oph.particle_candidate_charged_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('charged candidate inputs changed during replay')
    return result
