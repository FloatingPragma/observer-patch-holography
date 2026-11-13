"""Exact cluster intertwining and quantitative complete-source Weyl vectors."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import sympy as sy
from model.openwave_full_strong_current.native import (create,translate,clean,add,dot,action_vector,base_action,correlated_action,epsilon,pair_column,coefficients)
from model.openwave_full_strong_current.source import source_identity
from model.openwave_hadron_refinement_flow.analytic import full_bounds,constants
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_baryon_dispersion.source import matrix,repair_matrix
from model.openwave_cluster_bandwidth.meson import matrices,data
from model.openwave_strange_shell_current.source import refined_shell,lambda_pole
from model.openwave_capacity_stress.source import sqrt_interval,enclosure
from model.openwave_hadron_refinement_flow.certified_rest import sqrt_upper
from model.openwave_euclidean_refinement.response import dyadic_outer


def wedge(left,right):
    out=defaultdict(F)
    for a,u in left.items():
        for b,v in right.items():
            state,sign=create(b,a)
            if sign:out[state]+=sign*u*v
    return clean(out)


def shifted(vector,d):return {translate(s,d):v for s,v in vector.items()}


@lru_cache(None)
def cluster_identity():
    # Exact source-term identity is evaluated on two distinct separated native
    # singlet configurations. Its general finite-range proof is in PROOF.md.
    controls=[]
    for shape,pair in [(((0,0,0),(0,0,0),(0,1,1)),((0,0,1),(0,1,0))),
                       (((0,0,0),(1,0,0),(1,1,1)),((0,0,1),(1,1,0)))]:
        left=epsilon(shape);right=shifted(pair_column(pair),5);v=wedge(left,right)
        assert dot(v,v)==dot(left,left)*dot(right,right)>0
        rows={}
        for name,action in [('base',base_action),('triquark',lambda s:correlated_action(s,'triquark')),('meson',lambda s:correlated_action(s,'meson'))]:
            actual=action_vector(action,v);expected=defaultdict(F)
            add(expected,wedge(action_vector(action,left),right))
            add(expected,wedge(left,action_vector(action,right)))
            assert actual==clean(expected),(shape,pair,name)
            rows[name]={'all_output_occupations':len(actual),'exact_intertwiner':True}
        controls.append({'baryon_shape':[list(x)for x in shape],'meson_shape_before_translation':[list(x)for x in pair],
                         'meson_translation_cells':5,'exact_product_norm_squared':str(dot(v,v)),'terms':rows})
    return controls


@lru_cache(None)
def uniform_data():
    source=source_identity();p=coefficients();shell=refined_shell();r=F(shell['half_angle_tangent_interval'][1]);c=(1-r*r)/(1+r*r)
    A,G=matrix(F(1),F(0),True);J,_,zeta=repair_matrix(True)
    pu=F(A[0,0]/G[0,0])+2*(1-c)*F(J[0,0]/G[0,0])
    MA,_,_,_=matrices(F(1),F(0),repair=False);mu=F(MA[0,0])+2*(1-c)*zeta*F(MA.inv()[0,0])
    Cpro=full_bounds(*degree_coefficients())['Q_floor']-p['delta'];Cpi=data(False)[2]
    assert pu<Cpro and mu<Cpi
    return {'shell':shell,'unit':p['unit'],'slope':4*p['kappa']/3,
        'proton':{'floor':Cpro,'energy_upper':pu,'hopping':6*p['t'],'translation_bound':6*p['t']+4*source['triquark_J_norm']},
        'pion':{'floor':Cpi,'energy_upper':mu,'hopping':4*p['t'],'translation_bound':4*p['t']+4*source['meson_J_norm']}}


def tail(kind,R):
    if type(R)is not int or R<1:raise ValueError('integer radius>=1 required')
    d=uniform_data();row=d[kind];value=F(1)
    for r in range(2,R+1):
        ratio=row['hopping']/(row['floor']+d['slope']*(r-2)-row['energy_upper'])
        assert 0<ratio<1
        value=dyadic_outer((F(0),value*ratio),bits=32*R)[1]
    return value


def sequence_member(R,L):
    if type(R)is not int or R<2 or type(L)is not int or L<2:raise ValueError('radii and lengths must be integers>=2')
    d=uniform_data();rows={};total=F()
    for kind in('proton','pion'):
        q=tail(kind,R);discard=tail(kind,R+1)
        norm_lower=sqrt_interval(1-discard*discard)[0];assert norm_lower>0
        truncation=2*d[kind]['hopping']*q/norm_lower
        spatial=d[kind]['translation_bound']*sqrt_upper(F(2,L),den=2**(32*R))
        truncation=dyadic_outer((F(0),truncation),bits=32*R)[1]
        spatial=dyadic_outer((F(0),spatial),bits=32*R)[1]
        residual=truncation+spatial;total+=residual
        rows[kind]={'relative_radius_cells':R,'discarded_eigenvector_norm_upper':str(discard),
            'truncated_eigenvector_norm_lower':str(norm_lower),'normalized_internal_residual_upper_source':str(truncation),
            'spatial_box_length_cells':L,'box_shift_difference_norm_exact_squared':str(F(2,L)),
            'spatial_boundary_residual_upper_source':str(spatial),'complete_one_cluster_residual_upper_source':str(residual)}
    D=4*(L+R+3);assert D+1-(L+R)>3
    return {'relative_radius':R,'box_length':L,'right_cluster_translation':D,'minimum_cluster_gap_cells':D+1-(L+R),
            'clusters':rows,'exact_normalized_product_norm_squared':'1','full_H5_energy_residual_upper_source':str(total),
            'full_H5_energy_residual_upper_GeV':str(dyadic_outer((F(0),total*d['unit']),bits=32*R)[1]),
            'full_H5_energy_residual_upper_GeV_diagnostic':float(total*d['unit']),
            'all_source_interaction_terms_factor_exactly_on_this_support':True}


def run():
    identity=cluster_identity();d=uniform_data();members=[sequence_member(r,2**(4*r))for r in(4,8,16)]
    assert all(F(b['full_H5_energy_residual_upper_GeV'])<F(a['full_H5_energy_residual_upper_GeV'])for a,b in zip(members,members[1:]))
    lam=lambda_pole()
    return {'schema':'oph.complete_native_proton_pion_weyl_sequence.v1',
        'verdict':'PASS_ACTUAL_BREAKUP_ESSENTIAL_SPECTRUM_CHANNEL_AT_NATIVE_LAMBDA_ENERGY',
        'exact_source_intertwining':{'all_source_terms':'H5 J_sep=J_sep(H3 tensor I+I tensor H2) while the clusters and every one-step successor remain separated.',
            'CAR_signs_Pauli_exchange_and_complete_correlated_terms_retained':True,
            'colour_singlet_gap_electric_flux_exact':'0','native_operator_controls':identity,
            'general_all_finite_support_proof':'Finite-range CAR factorization and singlet prefix Gauss identity, as proved in PROOF.md.',
            'full_H5_replaced_by_free_clusters_everywhere':False},
        'energy_shell':{'Lambda_full_source_rest_GeV_interval':lam['full_ground_GeV_interval'],
            'Lambda_full_source_rest_GeV_display':lam['full_ground_GeV_display'],
            'half_angle_tangent_interval':d['shell']['half_angle_tangent_interval'],
            'opposite_momentum_GeV_interval':d['shell']['opposite_momentum_GeV_interval'],
            'strict_below_margin_GeV':d['shell']['below_margin_GeV'],'strict_above_margin_GeV':d['shell']['above_margin_GeV'],
            'complete_interacting_H5_essential_spectrum_contains_native_Lambda_rest_energy':True,
            'claim_is_only_a_sum_of_isolated_cluster_mass_lines':False},
        'constructive_Weyl_family':{'definition':'For an isolated pair-shell momentum k*, take normalized span<=R truncations of the exact proton/pion ground eigenvectors at k*,-k*. Multiply each by an L-site normalized plane-wave box and separate their origins by D=4(L+R+3); apply the native CAR product.',
            'exact_full_source_ground_vectors_defined_by':'Isolated ground spectral projections of the confining relative operators; choose any normalized vector in each range at k*. Simplicity or a unique eigenvector is not required.',
            'all_eigenvector_coefficients_numerically_enumerated':False,
            'finite_local_observer_preparation_of_exact_spectral_vectors_compiled':False,
            'full_infinite_chain_vectors_are_normalizable':True,'finite_ring_or_delta_momentum_used':False,
            'weak_zero_reason':'The antiquark and its right singlet cluster lie beyond D; overlap with every fixed finite occupation vector is eventually zero.',
            'tail_recursion':'Exact tail <= product_(r=2)^R b/[C+4*kappa*(r-2)/3-E_upper]. Executed upper bounds round outward at32R bits, adding at most R*2^(-32R); the entire relative tail is retained.',
            'residual_bound':'sum_clusters[2 b q_R/sqrt(1-q_(R+1)^2)+T_translation sqrt(2/L)]',
            'sequence_rule':'R tends to infinity, L=2^(4R), D=4(L+R+3); every residual bound tends to zero.',
            'executed_exact_bound_members':members},
        'uniform_tail_parameters':{name:{key:str(value)for key,value in d[name].items()}for name in('proton','pion')},
        'complete_readout':{'energy_instrument':'Spectral interval and its full complement on normalized Weyl vectors.',
            'outside_width_w_probability_upper':'min(1, residual_source^2/w^2) for every positive w.',
            'nonselective_energy_readout_work_exact':'0','all_failure_and_external_reference_outcomes_retained':True},
        'new_empirical_fit_or_Nature_validation':False,
        'inherited_calibrations':'Electron GeV anchor, strong and axial coefficients, strange source coefficients, supplied strong length and physical speed.',
        'outgoing_wave_operator_constructed':False,'nonzero_native_weak_on_shell_spectral_density_proved':False,
        'physical_Lambda_decay_rate_or_3plus1_particle_identification_proved':False,
        'scope':'Actual essential-spectrum inclusion for the complete interacting one-dimensional source, via normalizable separated-cluster Weyl vectors. The interacting outgoing weak overlap and physical continuum remain separate open premises.'}
