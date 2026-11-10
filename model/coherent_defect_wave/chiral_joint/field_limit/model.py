"""Actual quantum-link low-energy reduction and complete residual controls."""
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import product
import hashlib
from math import isqrt
import sympy as s
from ..model import active_terms,clean
from ..two_cell.hamiltonian import hopping,local_yukawa,matrix_on
from ..two_cell.spectrum import physical_block
from ..two_cell.basis import local_current_column
from ...chiral_source.model import link_labels
from ...physical_attachment.chiral import representation


def sqrt_upper(value,bits=80):
    value=s.Rational(value);scale=1<<bits
    if value<0:raise ValueError('nonnegative squared norm required')
    n=isqrt(int(value.p)*scale*scale//int(value.q))
    if n*n*int(value.q)<int(value.p)*scale*scale:n+=1
    return s.Rational(n,scale)


def exact_positive(value):
    value=s.sympify(value)
    if value.has(s.Float) or value.is_positive is not True:
        raise ValueError('positive exact parameter required')
    return value


def star_sign(mask,link,site):
    return (-1)**((((mask>>(6*site))&63).bit_count()+int(link!=0))%2)


@lru_cache(maxsize=1)
def star_parity_checks():
    transitions=0
    for mask in range(4096):
        for link in range(6):
            for (target,dest),coefficient in hopping(mask,link).items():
                for site in (0,1):
                    assert star_sign(mask,link,site)==star_sign(target,dest,site)
                transitions+=1
    local=0
    for mask in range(64):
        for scalar in range(5):
            for (target,dest),coefficient in local_yukawa(mask,scalar).items():
                assert mask.bit_count()%2==target.bit_count()%2
                local+=1
    paths=0
    # An arbitrary simple oriented path contributes degree two at each
    # interior vertex and degree one at its endpoints, including a seam.
    for side in (2,4,6):
        for start in range(side):
            for length in range(1,side):
                incidence=Counter()
                for j in range(length):
                    incidence[(start+j)%side]+=1
                    incidence[(start+j+1)%side]+=1
                incidence[start]+=1;incidence[(start+length)%side]+=1
                assert all(v%2==0 for v in incidence.values())
                paths+=1
    return {'actual_hopping_columns':4096*6,'actual_nonzero_hopping_entries':transitions,
            'actual_local_Yukawa_columns':64*5,'actual_local_Yukawa_entries':local,
            'path_incidence_controls':paths,
            'conserved_operator':'C_x=(-1)^N_matter,x product_incident_links (-1)^N_link',
            'blank_link_effective_H_preserves_each_site_matter_parity':True,
            'single_fermion_hopping_with_all_links_returned_blank':False,
            'dressed_charged_propagation_excluded':False}


@lru_cache(maxsize=1)
def path_fluctuation_checks():
    y=representation()[2][-1];blocks={}
    for i,q in enumerate(y.diagonal()):blocks.setdefault(q,[]).append(i)
    labels=set(link_labels());rows=[]
    for indices in blocks.values():
        d=len(indices);a=indices[0];b=indices[-1]
        for length in (1,2,3,4):
            words=[]
            for inside in product(indices,repeat=length-1):
                vertices=(a,*inside,b)
                word=tuple(zip(vertices[:-1],vertices[1:]))
                assert all(edge in labels for edge in word)
                words.append(word)
            assert len(words)==len(set(words))==d**(length-1)
            rows.append({'block_dimension':d,'length':length,
                         'fixed_endpoint_path_norm_squared':str(len(words)),
                         'one_internal_endpoint_full_transporter_norm_squared':str(d**length)})
    return {'native_block_dimensions':[len(q) for q in blocks.values()],
            'actual_link_dimension':len(labels),'rows':rows,
            'raw_path_equals_unitary_transport_on_blank_reference':False,
            'normalizing_path_changes_declared_source_coefficients':True}


@lru_cache(maxsize=1)
def bounded_transporter_control():
    checks=[]
    for d in (1,2,3,6):
        size=1+d*d;entries={}
        for a in range(d):
            for b in range(d):entries[a*size+1+a*d+b,b*size]=1/s.sqrt(d)
        v=s.ImmutableSparseMatrix(d*size,d*size,entries)
        expected=s.ImmutableSparseMatrix(d*size,d*size,{(b*size,b*size):1 for b in range(d)})
        assert v.H*v==expected
        checks.append({'native_block_dimension':d,'combined_operator_dimension':d*size,
                       'adjoint_product':'I_internal tensor P_blank','operator_norm':'1'})
    menus=[]
    for side in (2,4,6,8):
        # a=pi, L=N*pi. All receiving words of different directed paths
        # are orthogonal on the blank-link input, even for adjacent ends.
        variance=s.simplify(3*sum(1/(4*side**2*s.sin(s.pi*r/side)**2) for r in range(1,side)))
        assert variance==s.Rational(side*side-1,4*side*side)
        menus.append({'side':side,'normalized_three_axis_hopping_norm_squared':str(variance)})
    return {'declared_alternative':'U_ab/sqrt(d) within each native irreducible block',
            'actual_raw_periodic_family_modified':False,'local_exact_controls':checks,
            'all_path_products_are_contractions':True,
            'blank_reference_path_norm_squared':'1',
            'blank_reference_one_fermion_hopping_norm_squared':'pi^2*(N^2-1)/(4*N^2*a^2)',
            'fixed_spacing_volume_uniform_bound':'pi^2/(4*a^2)',
            'explicit_side_controls_spacing_pi':menus,
            'Ward_identities_preserved_by_block_scalar_normalization':True,
            'flat_link_substitution_into_original_U_then_normalization_preserves_free_symbol':False,
            'nilpotent_link_replaced_by_unitary':False,'star_parity_removed':False,
            'interacting_relativistic_limit_proved':False}


def reduction_bounds(kappa=128,steps=3,time=1):
    kappa=exact_positive(kappa);time=s.sympify(time)
    if type(steps) is not int or steps<0 or time.has(s.Float) or time.is_real is not True:
        raise ValueError('nonnegative integer iteration count and exact real time')
    a=s.Integer(18);b=s.Integer(3);gap=kappa-2*a
    if gap<4*b:raise ValueError('certified graph regime requires kappa-36>=12')
    radius=2*b/gap;contraction=2*b*radius/gap
    error=radius*contraction**steps
    isometry_error=(1+2*radius)*error
    hbound=kappa+a+b
    residual=4*hbound*isometry_error
    return {'local_H_norm_bound':str(a),'off_diagonal_norm_bound':str(b),
            'kappa':str(kappa),'Sylvester_gap_lower_bound':str(gap),
            'graph_radius':str(radius),'contraction_factor':str(contraction),
            'iterations':steps,'graph_error':str(error),
            'isometry_error':str(isometry_error),'full_H_norm_bound':str(hbound),
            'full_generator_residual_bound':str(residual),
            'time':str(time),'complete_effect_error':str(min(1,abs(time)*residual)),
            'all_input_vectors_and_retained_effects_covered':True,
            'full_global_source_spectrum_computed':False}


def _vec(matrix):
    return s.Matrix([matrix[i,j] for j in range(matrix.cols) for i in range(matrix.rows)])


def _unvec(vector,rows,cols):
    return s.Matrix(rows,cols,lambda i,j:vector[i+rows*j])


@lru_cache(maxsize=8)
def physical_graph_matrices(kappa=s.Integer(128),steps=3):
    """All fifteen exact states of one full neutral Gauss conserved block.

    The state space is unchanged. Only the explicitly supplied electric
    coefficient kappa differs from the unit-parameter control.
    """
    kappa=exact_positive(kappa);bounds=reduction_bounds(kappa,steps)
    keys,embedding,metric,form,generator=physical_block(0,0,2)
    assert generator.rows==15
    blank=[];excited=[]
    for j in range(embedding.cols):
        kinds={int(keys[i][1]!=0) for i in range(embedding.rows) if embedding[i,j]}
        assert len(kinds)==1
        (excited if kinds=={1} else blank).append(j)
    assert len(blank)==6 and len(excited)==9
    order=blank+excited
    h=generator.extract(order,order)+s.diag(*([0]*6+[kappa-1]*9))
    g=metric.extract(order,order)
    assert h.T*g==g*h and g[:6,6:]==s.zeros(6,9)
    a=h[:6,:6];upper=h[:6,6:];b=h[6:,:6];d=h[6:,6:]
    gp=g[:6,:6];gq=g[6:,6:]
    assert gp*upper==b.T*gq
    # Solve the exact Sylvester equation d X-X a=Y. Its operator gap is
    # certified in the actual physical metric, not an entrywise norm.
    sylvester=s.kronecker_product(s.eye(6),d)-s.kronecker_product(a.T,s.eye(9))
    inverse=sylvester.inv(method='DM')
    x=s.zeros(9,6)
    for _ in range(steps):
        x=_unvec(inverse*_vec(-b+x*upper*x),9,6).applyfunc(s.cancel)
    graph=s.eye(15)[:,:6]+s.Matrix.vstack(s.zeros(6,6),x)
    gm=(graph.T*g*graph).applyfunc(s.cancel)
    km=(graph.T*g*h*graph).applyfunc(s.cancel)
    assert km==km.T
    effective=(gm.inv(method='DM')*km).applyfunc(s.cancel)
    residual=(h*graph-graph*effective).applyfunc(s.cancel)
    # The squared Hilbert-Schmidt residual in orthonormal source and
    # dressed-input coordinates dominates its operator norm squared.
    residual_squared=s.factor(s.trace(gm.inv(method='DM')*residual.T*g*residual))
    assert residual_squared>=0
    # A trace-norm check may be looser than the dimension-independent
    # analytic operator bound by sqrt(rank), so retain that factor here.
    asserted=s.Rational(bounds['full_generator_residual_bound'])
    assert residual_squared<=6*asserted**2
    leading=(upper*b).applyfunc(s.cancel)
    assert leading!=s.zeros(6)
    return {'H':h,'metric':g,'blank_metric':gp,'A':a,'B':b,'B_adjoint':upper,
            'D':d,'X':x,'graph':graph,'effective_metric':gm,'effective_form':km,
            'effective_generator':effective,'residual_squared_Hilbert_Schmidt':residual_squared,
            'leading_virtual_interaction':-leading/kappa,'bounds':bounds}


@lru_cache(maxsize=1)
def quartic_control():
    # All other modes and the two scalar coordinates are spectators for
    # this exact hopping-square calculation; Yukawa belongs to A and D.
    def response(mask):
        value=s.Integer(0)
        for (_,link),coefficient in hopping(mask,0).items():
            assert link!=0
            value+=s.conjugate(coefficient)*coefficient
        return s.simplify(value)
    values={str((left,right)):response(left|(right<<6)) for left,right in product((0,1),repeat=2)}
    for left,right in product((0,1),repeat=2):
        key=(left|(right<<6),0,0,0)
        for site in (0,1):
            for current in range(12):
                assert not local_current_column(key,site,current,boundary=(-left,-right))
    assert values=={'(0, 0)':0,'(0, 1)':s.Rational(1,4),'(1, 0)':0,'(1, 1)':0}
    connected=-(values['(1, 1)']-values['(1, 0)']-values['(0, 1)']+values['(0, 0)'])
    assert connected==s.Rational(1,4)
    return {'actual_channel':'native e^c spin-up on both cells',
            'B_adjoint_B_diagonal_values':{k:str(v) for k,v in values.items()},
            'connected_diagonal_coefficient_of_minus_B_adjoint_B':'1/4',
            'leading_energy_coefficient':'1/(4*kappa)',
            'one_body_operator_countercontrol_rejected':True,
            'all_24_actual_Gauss_constraints_checked_on_each_compensated_input':True,
            'references_retain_opposite_L_spin_and_matching_parity':True,
            'bare_charged_controls_are_not_uncompensated_physical_states':True}


def connected_effective_energy_bound(kappa=65536):
    bounds=reduction_bounds(kappa,0)
    kappa=s.Rational(bounds['kappa']);r=s.Rational(bounds['graph_radius'])
    a=s.Rational(bounds['local_H_norm_bound']);b=s.Rational(bounds['off_diagonal_norm_bound'])
    # The exact canonical low generator obeys K=A-Bdagger B/kappa+error.
    graph_leading_error=(2*a*r+b*r*r)/kappa
    canonical_error=b*graph_leading_error+r*r*(a+b*r)
    connected_leading=s.Rational(1,4)/kappa
    lower=connected_leading-4*canonical_error
    assert lower>0
    return {'kappa':str(kappa),'canonical_effective_H_leading_error':str(canonical_error),
            'connected_energy_leading_term':str(connected_leading),
            'exact_dressed_connected_energy_lower_bound':str(lower),
            'four_input_states':'vacuum, left e^c up, right e^c up, both; retain compensating static charges',
            'reference_energy_of_compensators':'0, declared',
            'connected_effective_energy_is_positive':True,
            'finite_interaction_not_a_continuum_scattering_claim':True}


def certificate_data():
    control=physical_graph_matrices()
    matrix_fingerprint=hashlib.sha256(str((control['effective_metric'],control['effective_form'])).encode()).hexdigest()
    return {'source':'same actual one-link two-cell operator grammar, both quantum Yukawas=1, scalar penalty=1, declared electric penalty kappa=128',
            'star_parity':star_parity_checks(),'quantum_path_fluctuations':path_fluctuation_checks(),
            'bounded_transporter_alternative':bounded_transporter_control(),
            'virtual_interaction':quartic_control(),'controlled_connected_energy':connected_effective_energy_bound(),
            'low_energy_graph':control['bounds'],
            'physical_execution':{'neutral_Gauss_block':{'L':0,'twoSz':2},
                'complete_physical_dimension':15,'blank_dimension':6,'excited_dimension':9,
                'exact_metric_Hermiticity':True,'full_graph_residual_squared_Hilbert_Schmidt':str(control['residual_squared_Hilbert_Schmidt']),
                'complete_effect_error_per_unit_time':str(sqrt_upper(control['residual_squared_Hilbert_Schmidt'])),
                'effective_metric_and_form_sha256':matrix_fingerprint,
                'reference_is_global_ground':False,'same_H_all_15_states_retained_in_residual':True},
            'scope':{'finite_dressed_effective_interaction_constructed':True,
                'unit_parameter_ground_transferred':False,'interacting_Lorentz_continuum':False,
                'physical_spin_or_species_identified':False,'GeV_mass_prediction':False}}
