"""Original neutral two-point moments bounded by source occupied pairs."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
import sympy as s
from scipy.sparse import csr_matrix

from model.particle_common_contract import contract, provenance
from model.particle_primary_selected_spin_response import model as spin
from model.particle_primary_coherent_odd_response import model as coherent
from model.particle_primary_odd_sum_rules import model as odd_rules
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_original_shell_moments import model as clock_means

SOURCE_INPUTS=tuple(dict.fromkeys((
    'model/particle_candidate_neutral_dispersion/__init__.py',
    'model/particle_candidate_neutral_dispersion/__main__.py',
    'model/particle_candidate_neutral_dispersion/model.py',
    'model/particle_candidate_neutral_dispersion/test_model.py',
    'model/particle_candidate_neutral_dispersion/PROOF.md',
    'model/particle_candidate_neutral_dispersion/README.md',
    *spin.SOURCE_INPUTS,*coherent.SOURCE_INPUTS,*odd_rules.SOURCE_INPUTS,
    *spectroscopy.SOURCE_INPUTS,*clock_means.SOURCE_INPUTS,
)))


def derivative(polynomial, index):
    out={}
    for word,coefficient in polynomial.items():
        if index in word:
            position=word.index(index)
            pair=word[:position]+word[position+1:]
            out[pair]=s.expand(out.get(pair,0)+(-1)**position*coefficient)
    return {word:value for word,value in out.items() if value}


def spectrum(matrix):
    return {str(value):int(multiplicity) for value,multiplicity in matrix.eigenvals().items()}


@lru_cache(maxsize=1)
def occupied_pair_control():
    """Full local CAR identities and exact occupation inequalities."""
    local=spin.local_controls();number=local['positive_number'];rows=[]
    annihilators=[spin.odd.creation_word((i,),6).H for i in range(6)]
    for label,family in local['families'].items():
        d=family['dimension'];J=s.zeros(64);derivatives=0
        for polynomial,A in zip(family['polynomials'],family['operators']):
            for i in range(6):
                B=spin.clean(sum((value*spin.odd.creation_word(word,6)
                    for word,value in derivative(polynomial,i).items()),s.zeros(64)))
                if B!=spin.clean(annihilators[i]*A+A*annihilators[i]):
                    raise ArithmeticError('actual CAR derivative sign differs')
                J+=B*B.H;derivatives+=1
        J=spin.clean(J);deficit=spin.clean(s.Rational(d,2)*number-J)
        if J.H!=J or J*number!=number*J:
            raise ArithmeticError('occupied pair effect is not number preserving')
        if any(value<0 for value in J.eigenvals()) or any(value<0 for value in deficit.eigenvals()):
            raise ArithmeticError('occupied pair number inequality failed')
        if J[63,63]!=3*d or deficit[63,63]!=0:
            raise ArithmeticError('sharp occupied pair constant changed')
        sectors=[]
        for n in range(7):
            indices=[i for i in range(64) if i.bit_count()==n]
            sectors.append({'occupation':n,'dimension':len(indices),
                'occupied_pair_spectrum':spectrum(J.extract(indices,indices))})
        rows.append({'probe_irrep':label,'dimension':d,'derivative_identities_checked':derivatives,
            'occupied_pair_number_coefficient':F(d,2),'occupation_sectors':sectors,
            'deficit_spectrum':spectrum(deficit),'coefficient_sharp_at_full_occupation':True})
    return {'local_CAR_dimension':64,'native_modes_per_site':60,
        'active_source_modes':local['active_source_modes'],
        'all_other_native_modes_retained_as_spectators':True,
        'effect_order':'sum B B*, occupied pairs','families':rows}


def _integer_operator(polynomial, offset=0):
    result=csr_matrix((4096,4096),dtype=np.int64)
    for word,coefficient in polynomial.items():
        result+=int(coefficient)*spin._creation_sparse(tuple(i+offset for i in word),12)
    return result


@lru_cache(maxsize=1)
def neighbor_car_control():
    """Exact full twelve-mode identity for each actual active link entry."""
    # Two spin-preserving weak blocks and two singlet blocks.
    entries=tuple((i,k) for offset in (0,3)
        for block in ((offset,offset+1),(offset+2,)) for i in block for k in block)
    checks=0;reverse_checks=0
    for family in spin.local_controls()['families'].values():
        for polynomial in family['polynomials']:
            square,integer=spin._integer_polynomial(polynomial)
            if square<=0:raise ArithmeticError('normalized polynomial weight changed')
            Fx=_integer_operator(integer);Fy=_integer_operator(integer,6)
            for i,k in entries:
                hop=spin._creation_sparse((i,),12)@spin._creation_sparse((k+6,),12).T
                commutator=hop@Fy-Fy@hop
                actual=Fx.T@commutator+commutator@Fx.T
                Bx=_integer_operator(derivative(integer,i));By=_integer_operator(derivative(integer,k),6)
                difference=actual-Bx.T@By;difference.eliminate_zeros()
                if difference.nnz:raise ArithmeticError('neighbor graded CAR coefficient differs')
                reverse=hop.T@Fy-Fy@hop.T;reverse.eliminate_zeros()
                if reverse.nnz:raise ArithmeticError('reverse hop unexpectedly creates at source site')
                checks+=1;reverse_checks+=1
    actual_wilson=odd_rules.neighbor_graded_control()
    return {'full_CAR_dimension':4096,'active_modes_per_site':6,
        'native_modes_per_site':60,'allowed_active_link_entries':entries,
        'exact_graded_coefficient_identities':checks,'reverse_commutators_exactly_zero':reverse_checks,
        'actual_source_axes_checked':actual_wilson['actual_source_axes_checked'],
        'actual_positive_beta_Wilson_coefficient':F(-1,2),
        'unrestricted_occupation_sectors':True,
        'positive_witness_norm_used_as_signed_expectation':False,
        'link_and_Gauss_cutoffs_applied_after_parent_CAR_identity':True}


def centroid_bounds(dimension,delta,residual_cap):
    d,delta,C=map(F,(dimension,delta,residual_cap))
    if d not in (2,4) or not 0<=delta<5 or C not in (19,36):
        raise ValueError('member-zero neutral multiplet constants required')
    if (d,C) not in ((F(2),F(19)),(F(4),F(36))):
        raise ValueError('residual cap belongs to a different multiplet')
    w=1-delta/32;z=d*w;m1=d*delta/64
    adjacent=2*m1/z;spread=6*m1/z;spatial=3*m1/z
    if not (w>F(27,32) and C/z<F(27,8)**2 and adjacent<F(5,27)
            and spread<F(5,9) and spatial<F(5,18)):
        raise ArithmeticError('original centroid bounds failed')
    radius=F(27,8)+F(5,18)
    if radius!=F(263,72) or radius>=F(11,3):
        raise ArithmeticError('coarse centroid interval failed')
    return {'dimension':int(d),'source_defect_expectation_upper':delta,
        'signed_CAR_weight_lower':z,'signed_CAR_weight_upper':d,
        'neighbor_first_moment_absolute_upper':m1,
        'onsite_centered_residual_squared_upper':C,
        'onsite_centroid_squared_distance_from_nine_upper':C/z,
        'adjacent_centroid_difference_absolute_upper':adjacent,
        'centroid_span_upper':spread,'simple_centroid_span_strict_upper':F(5,9),
        'centroid_strict_interval':[F(16,3),F(38,3)],
        'within_multiplet_centroid_ratio_strict_upper':F(53,48),
        'neighbor_first_moment_value_or_sign_computed':False,
        'centroids_are_isolated_particle_energies':False}


def complete_mean_estimator(member,epsilon=F(1,1024)):
    epsilon=F(epsilon)
    if not 0<epsilon<1:raise ValueError('positive absolute tolerance below one required')
    source=member['preparation_plan']['source']
    R=max(abs(F(source['spectral_lower'])),abs(F(source['spectral_upper'])))
    eta=epsilon/128;allowance=epsilon/(14336*R)
    plan=spectroscopy.clock_plan(member,eta,allowance,allowance)
    tail=F(plan['clock_plan']['complete_clock_tail_upper'])
    error=F(plan['two_complete_clock_error_upper'])
    row=clock_means.finite_record_error(R,eta,tail,error)
    if 16*row>epsilon/2:raise ArithmeticError('neighbor estimator deterministic bias exceeds target')
    return {'source_clock_plan':plan,'requested_absolute_tolerance':epsilon,
        'one_Walsh_row_mean_error_upper':row,'neighbor_estimator_bias_upper':16*row,
        'simple_neighbor_estimator_bias_upper':epsilon/2,
        'estimator':'8*sigma*omega_hat*(1[q=000]-1[q=111]) within the chosen multiplet; zero on other records',
        'ideal_identity':'m1=8*(a0-a3), a_s=(m0+(3-2s)*m1)/48',
        'complete_success_records':96,'retained_address_failure_records':8,
        'all_failures_both_energy_records_reference_and_memory_retained':True,
        'full_original_histogram_executed':False,
        'finite_sample_statistical_error_included':False,
        'source_clock_generator_or_energy_zero_changed':False}


def bind_member(member):
    contract.verify_member(member)
    if member['source']['level']!=0:
        raise ValueError('occupied-pair spatial moment bound is instantiated only at member zero')
    inherited=coherent.bind_member(member);source=inherited['result']
    delta=32*(1-F(source['original_total_signed_success_lower']))
    occupied_pair_control();neighbor_car_control()
    channels=[{'probe_irrep':label,**centroid_bounds(d,delta,C)}
        for label,d,C in [('G1',2,19),('H',4,36)]]
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'original_coherent_response_receipt_sha256':contract.digest(inherited),
        'original_polynomial_density_and_full_physical_source':True,
        'source_CAR_parity_even_and_odd_one_point_mean_zero':True,
        'native_increment':[0,2,0,0,1],'spatial_first_moment_shells':['m0','m1','0','0'],
        'Walsh_first_moments':['m0+3*m1','m0+m1','m0-m1','m0-3*m1'],
        'Walsh_multiplicities':[1,3,3,1],'complete_menu_moment_divisor':48,
        'channels':channels,'complete_neighbor_mean_estimator':complete_mean_estimator(member),
        'finite_Walsh_labels_are_conserved_operational_momentum':False,
        'nonconstant_original_odd_two_point_first_moments_proved':False,
        'positive_threshold_or_isolated_pole_proved':False,
        'relativistic_dispersion_or_refinement_limit_proved':False,
        'physical_species_map':[],'physical_mass_table_GeV':[]})
    return {'schema':'oph.particle_candidate_neutral_dispersion.v1',
        'result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('neutral moment result or fixed source binding changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.particle_candidate_neutral_dispersion_certificate.v1',
        'exact_occupied_pair_control':occupied_pair_control(),
        'exact_neighbor_CAR_control':neighbor_car_control(),
        'members':[bind_member(contract.member(0))],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:
        raise ValueError('source changed during neutral moment certificate')
    return contract.serial(value)
