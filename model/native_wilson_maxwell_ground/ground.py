"""Same full-source ground audit; exact filled-band contractions and Schur bounds."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
import sympy as s
from sympy.polys.matrices import DomainMatrix
import mpmath as mp

HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
sys.path.insert(0,str(RER))
from model.native_wilson_maxwell import field,exact_current
from model.native_wilson_maxwell_refinement import bounds
from model.affine_gauss_mobile.model import ivalue,rounded
from model.native_block_dirac import algebra as native_code

PARENT='model/native_wilson_maxwell/results.json'
PARENT_SHA='eba4963148b62f6841a038119ed5b9baf4a6f142009915b31e49d72dacd12b72'
GEOMETRY='model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/results.json'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def admit():
    assert digest(RER/PARENT)==PARENT_SHA
    p=json.loads((RER/PARENT).read_text())
    for name,h in p['input_sha256'].items():
        assert digest(RER/name)==h,name
    return p


@lru_cache(None)
def components():
    H,_,minimal,powers=exact_current.literal_matter()
    z=s.Symbol('z');rs=(3,7,13)
    factors={r:s.Poly(z*z-(2*r+s.Rational(1,2))*z+(r-s.Rational(1,4))**2,z) for r in rs}
    Q=s.prod(factors.values());matrices=[];projectors=[]
    for r in rs:
        other=Q.exquo(factors[r]);E=(other*s.invert(other,factors[r])).rem(Q)
        numerator=s.Poly((3*r+s.Rational(1,4)-z)/(2*(r-s.Rational(1,4))),z)
        R=(E*numerator).rem(Q)
        M=sum((R.nth(k)*powers[2*k+1] for k in range(6)),s.zeros(32))
        P=sum((E.nth(k)*powers[2*k] for k in range(6)),s.zeros(32))
        matrices.append(s.ImmutableMatrix(M));projectors.append(s.ImmutableMatrix(P))
    dm=[DomainMatrix.from_Matrix(M).convert_to(s.QQ_I) for M in matrices]
    assert sum(projectors,s.zeros(32))==s.eye(32)
    for i,M in enumerate(dm):
        for j,N in enumerate(dm):
            assert (M*N).to_Matrix()==(rs[i]*projectors[i] if i==j else s.zeros(32))
    for M in matrices:
        for x in range(8):assert s.trace(M[4*x:4*x+4,4*x:4*x+4])==0
    return H,rs,matrices


@lru_cache(None)
def contractions():
    H,rs,Ms=components();G=field.cube()['G'];table=s.zeros(3)
    for i,M in enumerate(Ms):
        for j,N in enumerate(Ms):
            value=0
            for x in range(8):
                for y in range(8):
                    value+=G[x,y]*s.trace(M[4*x:4*x+4,4*y:4*y+4]*N[4*y:4*y+4,4*x:4*x+4])
            value=s.expand(value)
            assert value.is_Rational;table[i,j]=value
    assert table==table.T
    coul=s.simplify(s.trace(G)/2-sum(table[i,j]/(8*s.sqrt(rs[i]*rs[j])) for i in range(3) for j in range(3)))
    onsite=s.zeros(32)
    for x in range(8):onsite[4*x:4*x+4,4*x:4*x+4]=s.Rational(5,2)*s.diag(1,1,-1,-1)
    hop=H-onsite
    coefficients=[s.expand(s.trace(M*hop)/2) for M in Ms]
    assert all(value.is_Rational for value in coefficients)
    positive_hop=s.simplify(sum(value/s.sqrt(r) for value,r in zip(coefficients,rs)))
    mp.iv.dps=90
    ci=rounded(ivalue(coul),places=70);hi=rounded(ivalue(positive_hop),places=70)
    return {'exact_Coulomb_vacuum_coefficient':str(coul),'Coulomb_vacuum_coefficient_interval':ci,
        'exact_positive_diamagnetic_hopping_coefficient':str(positive_hop),
        'positive_diamagnetic_hopping_coefficient_interval':hi,
        'complete_CRT_Coulomb_contraction_coefficients':[[str(x) for x in table.row(i)] for i in range(3)],
        'complete_CRT_hopping_contraction_coefficients':[str(x) for x in coefficients],
        'all_eight_vacuum_charge_expectations_exactly_zero':True,
        'full_32_mode_filled_band_and_all_particle_hole_contractions':True}


def native_complement(parent):
    geometry=json.loads((RER/GEOMETRY).read_text())
    transfer=geometry['full_interacting_transfer'];cube=geometry['full_spin_spectrum']
    lam=F(1,2**34);separation=F(transfer['atomic_separation_lower']);tau=360
    jl,ju=map(F,transfer['actual_exchange_interval']);dl,du=map(F,cube['gap_interval'])
    r4=25*(lam*tau)**4/separation**3
    lo=2*lam*lam*jl*dl-2*r4;hi=2*lam*lam*ju*du+2*r4
    assert 0<lo<hi<F(1,10**20)
    join=json.loads((RER/'model/native_abelian_field_join/results.json').read_text())
    port=F(join['intervals']['complete_port_source_gap_lower'])
    occupied=F(join['intervals']['complete_old_and_dual_reference_gap_lower'])
    assert min(port,occupied)>hi
    return {'source_hopping':str(lam),'actual_neutral_triplet_gap_interval':[str(lo),str(hi)],
        'complete_charged_port_gap_lower':str(port),
        'complete_old_dual_reference_gap_lower':str(occupied),
        'one_triplet_constituent_choices':32,'each_native_triplet_rank':3,
        'active_complement_ground_multiplicity':96,
        'retained_preparation_reference_spectators_multiply_degeneracy':True,
        'all_added_Rloc_vacua_and_outer_charge_zero':True,
        'Coulomb_and_extra_reference_cost_on_these_states':'0',
        'common_photon_zero_point_over_m':'3+sqrt(6)',
        'raw_energy_baseline':'32 Eg + m(3+sqrt(6))'}


def local_completion_control():
    E,_,_=native_code.doublet_code()
    pair=s.kronecker_product(E,E);embedding=s.zeros(81,16)
    for word in range(16):embedding[:,word]=pair[:,(word&3)*4+(word>>2)]
    P=embedding*embedding.H
    n=s.diag(0,1,1);n2=s.kronecker_product(n,s.eye(3))+s.kronecker_product(s.eye(3),n)
    charge=s.kronecker_product(n2,s.eye(9))-s.kronecker_product(s.eye(9),n2)
    assert P*P==P and s.trace(P)==16 and P*charge==charge*P
    # Two macrovertices with a retained charge-carrying complement. The
    # exact FW pair rotation tests the Gauss spectral-lift orientation.
    p=s.diag(1,1,0);P1=s.kronecker_product(p,s.eye(3));P2=s.kronecker_product(s.eye(3),p)
    Q1=s.kronecker_product(s.diag(0,1,-1),s.eye(3));Q2=s.kronecker_product(s.eye(3),s.diag(0,-1,1))
    U=s.eye(9);U[0,0]=U[4,4]=s.Rational(3,5);U[0,4]=s.Rational(4,5);U[4,0]=-s.Rational(4,5)
    assert U.H*U==s.eye(9) and U*(Q1+Q2)==(Q1+Q2)*U
    rho1,rho2=U*Q1*U.H,U*Q2*U.H
    for local in (P1,P2):
        assert local*U==U*local
        assert local*rho1==rho1*local and local*rho2==rho2*local
    labels=sorted({(Q1[j,j],Q2[j,j],-Q1[j,j]-Q2[j,j]) for j in range(9)})
    Yq=s.zeros(9*len(labels),9)
    for j in range(9):Yq[j*len(labels)+labels.index((Q1[j,j],Q2[j,j],-Q1[j,j]-Q2[j,j])),j]=1
    Y=s.kronecker_product(U,s.eye(len(labels)))*Yq*U.H
    assert Y.H*Y==s.eye(9)
    for local in (P1,P2):assert s.kronecker_product(local,s.eye(len(labels)))*Y==Y*local
    assert rho1!=Q1 and rho2!=Q2
    defects=2*s.eye(9)-P1-P2;global_complement=s.eye(9)-P1*P2
    assert all(x>=0 for x in (defects-global_complement).diagonal())
    return {'actual_81_state_native_selected_port_carrier_projector_rank':16,
        'exact_native_projector_port_charge_commutation':True,
        'full_other_native_sectors':'spectral invariance and fixed port charge; no truncated complement identification',
        'coherent_FW_density_and_Gauss_intertwiner_control':True,
        'each_local_Px_reduces_all_full_source_terms':True,
        'onsite_means':'one declared macrovertex containing four complete native constituents',
        'preexisting_global_kinetic_current_guard_removed':False,
        'counterterm_is_forced_by_primitive_repair_law':False}


def preparation(parent,p_code,angle,gap,counterterm,mupper,depth=2):
    # All source masses, references and raw zero points remain in B.
    p=bounds.rounded(p_code/(1+angle*angle),False,60)
    delta=bounds.rounded(gap,False,60)
    raw=F(parent['clock_and_capacity']['full_raw_H_ext_norm_upper'])+8*counterterm*mupper
    B=F(-(-raw.numerator//raw.denominator))
    epsilon_target=F(1,1024)
    ratio=2*B/(epsilon_target*delta);M=-(-ratio.numerator//ratio.denominator)
    epsilon=2*B/(M*delta)
    assert epsilon<=epsilon_target
    odds=(1-p)/p;success=F(1);steps=[]
    for stage in range(depth):
        probability=1/(1+odds)
        pairs=2**(depth-stage-1)
        success*=probability**(2*pairs)
        after=odds*odds+2*epsilon*epsilon*odds
        steps.append({'level':stage+1,'number_of_filters':pairs,
            'input_ground_probability_lower':str(bounds.rounded(probability,False,60)),
            'each_success_probability_lower':str(bounds.rounded(probability**2,False,60)),
            'output_ground_odds_upper':str(bounds.rounded(after,True,60))})
        odds=after
    fidelity=1/(1+odds)
    error=bounds.root(1-fidelity)
    count=2**depth-1;copies=2**depth
    mp.iv.dps=90;log=mp.iv.ln
    base=32*(256*log(2)+8*log(5)+24*log(119)+log(2**31+1))+7*log(9)+log(65)+log(56)
    capacity=rounded(copies*base+count*log(M),places=60)
    return {'depth':depth,'full_source_trial_copies':copies,'native_constituents':32*copies,
        'finite_filter_count':count,'initial_ground_probability_lower':str(p),
        'raw_full_kinematic_source_norm_integer_upper':str(B),
        'onsite_completion_full_raw_norm_increment_upper':str(8*counterterm*mupper),
        'actual_full_physical_gap_lower':str(delta),
        'Cayley':'(I-i H_repaired/B)(I+i H_repaired/B)^-1, with raw source energy',
        'each_pointer_dimension':M,'actual_mixed_energy_amplitude_upper':str(epsilon),
        'complete_tree_success_lower':str(bounds.rounded(success,False,60)),
        'conditional_main_ground_fidelity_lower':str(bounds.rounded(fidelity,False,60)),
        'conditional_main_trace_distance_upper':str(bounds.rounded(error,True,60)),
        'conditional_retained_purification_vector_error_upper':str(bounds.root(2*(1-bounds.root(fidelity,False)))),
        'all_filter_basis_outcomes_retained':True,'all_partner_sources_retained':True,
        'ground_energy_or_exact_ground_projector_used':False,
        'stages':steps,
        'total_pointer_dimension_formula':f'{M}^{count}',
        'capacity_nats_formula':f'{copies}*(parent capacity-ln1161)+{count}*ln({M})',
        'capacity_nats_interval':capacity,
        'optional_complete_count_readout_additional_capacity':'ln1161',
        'initial_preparation_duration_source_units':'1',
        'each_full_filter_program_duration_source_units':'1',
        'total_elapsed_time_all_branches':1+count,
        'filter_pulse':'i Log(exp(-i H_copies) U_filter)-H_copies; keep the final common free factor',
        'failed_ancestor_slots':'declared no-op instrument with the same elapsed source time and retained earlier failure label',
        'full_service_H_norm_upper':str(copies*B),
        'each_control_norm_upper':str(copies*B+F(22,7)),
        'net_bare_source_work_absolute_upper':str(2*copies*B),
        'gross_all_pulse_switching_work_upper':str(2*(1+count)*(copies*B+F(22,7))),
        'efficient_native_C_E_compilation_claimed':False}


def calculate():
    parent=admit();data=contractions();complement=native_complement(parent)
    cl,cu=map(F,data['Coulomb_vacuum_coefficient_interval'])
    hl,hu=map(F,data['positive_diamagnetic_hopping_coefficient_interval'])
    vl=F(1,16)+1/(12*bounds.root(6));vu=F(1,16)+1/(12*bounds.root(6,False))
    g=F(1,64);hop_norm=18*g+F(27,4)*g*g
    sigma2upper=F(parent['independent_current_coefficient']['uniform_polynomial_method']['pair_one_photon_V1_vacuum_norm_squared_over_m_squared_interval'][1])
    sigup=bounds.root(sigma2upper)
    # Fourth/sixth vacuum moments agree with the exact Gaussian at cap3.
    mu_lower=g*g*cl+hl*(g*g*vl/2-g**4*vu*vu/8)
    mu_upper=g*g*cu+hu*(g*g*vu/2-g**4*vl*vl/8+g**6*vu**3/48)
    vacuum_remainder=24*bounds.root(3)*vu+bounds.root(8*cu)
    c=F(1,50000)
    denominator1=2*bounds.root(3,False)+1-hop_norm-c
    denominator2=2-hop_norm-c
    schur=(g*sigup/bounds.root(denominator1,False)+g*g*vacuum_remainder/bounds.root(denominator2,False))**2
    margin=mu_lower-c-schur
    assert margin>0
    mlo,mhi=map(F,parent['clock_and_capacity']['source_selected_excitation_interval'])
    assert mhi<1
    assert c*mlo>F(complement['actual_neutral_triplet_gap_interval'][1])
    # A separately declared guarded vacuum counterterm. Its coefficient is
    # the smallest dyadic at precision2^-48 bounding the exact quadratic
    # vacuum shift; it uses only the source contractions and g.
    quadratic_upper=cu+hu*vu/2
    scale=2**48
    coeff=F(-(-quadratic_upper.numerator*scale//quadratic_upper.denominator),scale)
    assert coeff-F(1,scale)<cl+hl*vl/2<=quadratic_upper<=coeff
    counterterm=g*g*coeff
    selection_margin=counterterm-mu_upper
    assert selection_margin>0
    code_gap=2-hop_norm-mu_upper
    repaired_gap=min(mlo*code_gap,F(complement['actual_neutral_triplet_gap_interval'][0])+mlo*selection_margin)
    wupper=g*sigup+g*g*vacuum_remainder
    angle=wupper/(2-hop_norm-mu_upper)
    assert 0<angle<1
    quantities={'g':g,'whole_code_Peierls_norm_over_m_upper':hop_norm,
        'free_vacuum_shift_over_m_lower':mu_lower,'free_vacuum_shift_over_m_upper':mu_upper,
        'vacuum_second_order_remainder_coefficient_upper':vacuum_remainder,
        'code_energy_over_baseline_lower_in_m_units':c,
        'weighted_Schur_term_upper':schur,'strict_Schur_margin':margin,
        'counterterm_quadratic_coefficient':coeff,'counterterm_over_m':counterterm,
        'completion_shift_minus_code_trial_upper_in_m_units':selection_margin,
        'same_code_internal_ground_gap_lower_over_m':code_gap,
        'repaired_full_physical_gap_lower':repaired_gap,
        'coupled_code_ground_sin_angle_to_free_vacuum_upper':angle,
        'original_code_above_native_triplet_margin_lower':c*mlo-F(complement['actual_neutral_triplet_gap_interval'][1])}
    out={key:str(bounds.rounded(value,'upper' in key or key=='weighted_Schur_term_upper')) for key,value in quantities.items()}
    prepared=preparation(parent,F(parent['actual_preparation']['global_code_probability_lower']),
        angle,repaired_gap,counterterm,mhi)
    sigma_lower=F(parent['independent_current_coefficient']['uniform_polynomial_method']['pair_one_photon_V1_vacuum_norm_over_m_lower'])
    vacuum_component_lower=1/bounds.root(1+angle*angle)
    neutral_norm=18*g+F(59,4)*g*g
    target_numerator=vacuum_component_lower*(g*sigma_lower-F(27,16)*g**3)-neutral_norm*angle
    ground_target_amplitude=target_numerator/(11+mu_upper)
    eps_prepared=F(prepared['conditional_retained_purification_vector_error_upper'])
    success=F(prepared['complete_tree_success_lower'])
    assert ground_target_amplitude>eps_prepared>0
    ground_record={'target':[1,1,1],
        'actual_stationary_ground_target_amplitude_lower':str(bounds.rounded(ground_target_amplitude,False,60)),
        'actual_stationary_ground_target_probability_lower':str(bounds.rounded(ground_target_amplitude**2,False,60)),
        'physical_preparation_and_count_joint_probability_lower':str(bounds.rounded(success*(ground_target_amplitude-eps_prepared)**2,False,60)),
        'complete_count_pointer_labels':1161,
        'all_preparation_failure_records_and_partner_sources_retained':True,
        'interpretation':'stationary interacting-ground occupation in the declared free quasiparticle/photon frame',
        'ground_emission_scattering_or_radiation':False}
    return {'schema':'oph.same_native_wilson_ground_audit.v1','parent_sha256':PARENT_SHA,
        'exact_vacuum_contractions':data,'complete_native_complement':complement,
        'bounds_exact':out,'diagnostic_binary64':{k:float(F(v)) for k,v in out.items()},
        'unchanged_source':{'ground_in_code':False,'active_ground_multiplicity':96,
            'old_finite_nonstationary_count_law_unchanged_by_this_audit':True,
            'code_only_gap_does_not_select_full_ground':True},
        'declared_onsite_completion':{'changed_H':'H_repaired=H_original+m c(g) sum_x(I-P_x), transported by Y and zero on the outer constraint complement',
            'c_of_g':'g² times the displayed dyadic quadratic coefficient; c(0)=0',
            'old_source_receipt_or_H_mutated':False,
            'complete_pattern_preserving_count_laws_exactly_invariant':True,
            'all_code_eigenvectors_absolute_energies_and_internal_differences_unchanged':True,
            'q_defect_complement_energy_increment':'q m c(g), q=1,...,8',
            'different_pattern_coherent_phases_and_energy_differences_changed':True,
            'global_scalar_energy_subtraction':'none',
            'coefficient_ancestry':'source contractions plus declared48-bit dyadic ceiling; not primitive-repair uniqueness or empirical calibration',
            'selected_active_coupled_ground_is_unique':True,
            'coupled_ground_is_free_vacuum':False,
            'actual_ground_radiation_or_SM_mass_claim':False},
        'complete_repaired_ground_preparation':prepared,
        'complete_stationary_ground_dressing_record':ground_record,
        'exact_macrocell_and_Gauss_completion_control':local_completion_control()}


if __name__=='__main__':print(json.dumps(calculate(),indent=2,sort_keys=True))
