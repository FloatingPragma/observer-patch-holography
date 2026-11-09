"""Exact source-scoped phase and C/E orientation constraints."""
from functools import lru_cache
from itertools import combinations

import sympy as s

from model.common_source import native_family_record_memory as memory
from model.standard_model_source_map import koide_source as coherent


def clean(matrix):
    return s.ImmutableMatrix(matrix.applyfunc(lambda x:s.expand(s.radsimp(x))))


def _zero(matrix):
    return all(x == 0 for x in clean(matrix))


@lru_cache(maxsize=1)
def _native_geometry():
    data,(gc,ge)=memory.geometry()
    frame=data['frame']; metric=(frame.T*frame)[0,0]
    def restrict(g):
        return clean(frame.T*memory.objects.permutation_matrix(g)*frame/metric)
    c,e=restrict(gc),restrict(ge)
    r=clean(e*c); reverse=clean(c*e)
    if not all(_zero(x) for x in (c*c-s.eye(3),e*e-s.eye(3),r**3-s.eye(3),
                                  reverse-r.T,c*r*c-r.T,e*r*e-r.T)):
        raise ArithmeticError('the actual half-turn orientation identities failed')
    if s.simplify(s.trace(r)) != 0 or s.Matrix.hstack(*(x.reshape(9,1) for x in (s.eye(3),r,r*r))).rank() != 3:
        raise ArithmeticError('the actual cyclic commutant has the wrong dimension')
    group={tuple(range(12))}; todo=list(group)
    while todo:
        old=todo.pop()
        for g in (gc,ge):
            composed=tuple(g[old[j]] for j in range(12))
            if composed not in group:
                group.add(composed);todo.append(composed)
    if len(group) != 6:
        raise ArithmeticError('the actual two half-turns do not generate D3')
    axis=clean((s.eye(3)+r+r*r)/3)
    extra=next((g for g in sorted(data['rotations']) if not _zero(restrict(g)*axis-axis*restrict(g))),None)
    if extra is None or not _zero(axis*axis-axis) or s.simplify(s.trace(axis)-1)!=0:
        raise ArithmeticError('full carrier rotations did not move the cyclic axis')
    return {'C':c,'E':e,'R':r,'R_reverse':reverse,'axis':axis,
            'extra_permutation':extra,'extra':restrict(extra),
            'gC':gc,'gE':ge,'group_order':len(group),'carrier_rotation_count':len(data['rotations'])}


def native_geometry():
    # Cached values are immutable matrices or tuples; callers receive a fresh map.
    return dict(_native_geometry())


def reflection_identity():
    data=native_geometry();r=data['R'];reflection=data['C']
    a,br,bi=s.symbols('a br bi',real=True)
    c=a*s.eye(3)+br*(r+r*r)+s.I*bi*(r-r*r)
    if not _zero(reflection*c*reflection-c+2*s.I*bi*(r-r*r)) or _zero(r-r*r):
        raise ArithmeticError('reflection did not isolate the imaginary circulant coefficient')
    x,y,b=s.symbols('x y b',real=True)
    root=s.sqrt(2)*b*s.eye(3)+b*((x+s.I*y)*r+(x-s.I*y)*r*r)
    action=clean(root*root/3)
    expected=2*b*b*y*(s.sqrt(2)-x)*(r-r*r)/3
    if not _zero((action-action.conjugate())/(2*s.I)-expected):
        raise ArithmeticError('the action reflection coefficient differs')
    # The same literal source reflection exchanges z and conjugate(z).
    if not _zero(reflection*root*reflection-root.conjugate()):
        raise ArithmeticError('ordered C/E reversal is not phase conjugation')
    return {'actual_halfturn_group_order':data['group_order'],
            'actual_carrier_rotation_count':data['carrier_rotation_count'],
            'gC':list(data['gC']),'gE':list(data['gE']),
            'native_R':[[str(x) for x in row] for row in r.tolist()],
            'ordered_EC_is_inverse_of_ordered_CE':True,
            'reflection_conjugates_coherent_phase':True,
            'unordered_action_spectrum_preserved_by_orientation_reversal':True,
            'C3_commutant_complex_dimension':3,
            'Hermitian_C3_commutant_real_dimension':3,
            'halfturn_invariant_Hermitian_commutant_real_dimension':2,
            'Fourier_active_labels_only_permute_unordered_action_spectrum':True,
            'reflection_invariant_balanced_positive_phase':'0 modulo 2*pi in the fixed R convention',
            'reflection_invariant_balanced_positive_response_has_repeated_pair':True,
            'full_A5_invariant_triplet_response_scalar':True,
            'extra_axis_moving_native_permutation':list(data['extra_permutation']),
            'full_A5_invariance_assumed_for_actual_reset_memory':False}


def fixed_memory_control():
    source=memory.architecture();data=native_geometry();r=data['R']
    projectors=source['projectors'];pi=source['stationary']
    records=memory.feedback.native_records()
    weights={name:s.zeros(12,1) for name in ('always','C','E','CE')}
    count=0
    for cc,cr,ec,er,_,weight in records:
        c=cc==0 and cr==1;e=ec==0 and er in (4,5)
        for old in range(12):
            destination=memory._destination(old,c,e,2,source['gC'],source['gE'])
            for name in weights:
                if name=='always' or name=='C' and c or name=='E' and e or name=='CE' and c and e:
                    weights[name][destination]+=pi[old]*weight
            count+=1
    coefficients={name:clean(sum((w*p for w,p in zip(values,projectors,strict=True)),s.zeros(3)))
                  for name,values in weights.items()}
    rows=[]
    for name,trace in zip(weights,(s.Integer(1),s.Rational(1,12),s.Rational(1,6),s.Rational(1,72)),strict=True):
        coefficient=coefficients[name]
        if not _zero(coefficient-coefficient.T) or any(s.im(v)!=0 for v in coefficient) or s.simplify(s.trace(coefficient)-trace)!=0:
            raise ArithmeticError('an original source duty coefficient changed')
        rows.append({'duty':name,'trace':str(trace),'rank':coefficient.rank(),
                     'real_symmetric':True,'commutes_with_actual_C3':_zero(coefficient*r-r*coefficient)})
    if [row['rank'] for row in rows] != [3,3,2,2] or rows[0]['commutes_with_actual_C3']:
        raise ArithmeticError('the actual memory rank/symmetry control changed')
    return {'control':'base native C/E law; ungated p=1/12,q=1/6,reset=2',
            'all_fine_records':len(records),'old_memory_inputs':12,'fine_memory_rows':count,
            'duties':rows,'raw_C_probability':'1/12','raw_E_probability':'1/6',
            'exchanging_all_C_E_record_roles_preserves_source_law':False,
            'memory_response_assumed_C3_invariant':False,
            'general_current_gated_positive_p_q_rank_result_owner':'model/native_flavor_orbit/algebra.py:memory_certificate'}


def ordered_feedback_grade():
    x=s.Matrix([[0,1],[1,0]]);z=s.diag(1,-1)
    x6=s.kronecker_product(x,s.eye(3));z6=s.kronecker_product(z,s.eye(3))
    projectors=memory.architecture()['projectors']
    cosine,sine=s.Rational(99,101),s.Rational(20,101)
    pulses=[]
    for p in projectors:
        if not _zero(p-p.T) or not _zero(p*p-p) or any(s.im(v)!=0 for v in p):
            raise ArithmeticError('an actual feedback projector left the real class')
        u=clean(s.eye(6)+(cosine-1)*s.kronecker_product(s.eye(2),p)
                -s.I*sine*s.kronecker_product(x,p))
        if not _zero(u.H*u-s.eye(6)) or not _zero(u*x6-x6*u) or not _zero(z6*u.conjugate()*z6-u):
            raise ArithmeticError('an actual feedback pulse lost its grade')
        pulses.append(u)
    ordered=clean(pulses[3]*pulses[2]);reverse=clean(pulses[2]*pulses[3])
    if _zero(ordered-reverse) or not _zero(z6*ordered.conjugate()*z6-ordered):
        raise ArithmeticError('ordering control or its preserved grade failed')
    return {'actual_real_projectors_checked':12,'commutes_with_Nambu_X':True,
            'Z_conjugation_Z_fixes_every_pulse_and_word':True,
            'noncommuting_ordered_pulse_control_retained':True,
            'differentiable_unitary_limit_generator':'I tensor iA + X tensor S; A real skew, S real symmetric',
            'pure_X_tensor_mass_in_fixed_native_frame_requires_real_family_matrix':True,
            'whole_fine_current_ancilla_archive_grading_proved':False,
            'all_family_conjugates_or_mass_spectra_excluded':False,
            'finite_step_reduced_channel_is_unitary':False}


def all_lag_control():
    """Reconstruct the finite all-lag potential from the actual twelve-state law."""
    source=memory.architecture()
    transition=s.Matrix(source['T']);pi=s.Matrix(source['stationary'])
    stationary=pi*s.ones(1,12)
    inverse=(s.eye(12)-transition+stationary).inv()-stationary
    if not _zero((s.eye(12)-transition)*inverse-s.eye(12)+stationary):
        raise ArithmeticError('the actual memory group inverse differs')
    k=s.diag(*pi)/2-pi*pi.T/2+transition*inverse*s.diag(*pi)
    projectors=source['projectors']
    potential=clean(-s.I/2*sum((k[i,j]*(projectors[i]*projectors[j]-projectors[j]*projectors[i])
                               for i in range(12) for j in range(12) if k[i,j]!=0),s.zeros(3)))
    expected=s.I*s.Rational(603,6666724)*s.Matrix(
        [[0,3+s.sqrt(5),0],[-3-s.sqrt(5),0,-2],[0,2,0]])
    if not _zero(potential-expected) or not _zero(potential-potential.H):
        raise ArithmeticError('the all-lag native ordering potential differs')
    rotation=native_geometry()['R']
    cyclic_coefficient=s.Rational(603,6666724)*(1+s.sqrt(5))
    if not _zero(potential*rotation-rotation*potential) or not _zero(
            potential-s.I*cyclic_coefficient*(rotation-rotation.T)):
        raise ArithmeticError('the actual ordering potential is not its cyclic orientation generator')
    diffusion=k+k.T
    support=tuple(i for i in range(12) if pi[i]!=0)
    principal=[diffusion.extract(indices,indices).det()
               for size in range(1,len(support)+1) for indices in combinations(support,size)]
    if any(x<0 for x in principal) or _zero(diffusion):
        raise ArithmeticError('the actual all-lag diffusion is not nonzero positive semidefinite')
    if any(diffusion[i,j]!=0 for i in range(12) for j in range(12)
           if i not in support or j not in support):
        raise ArithmeticError('a diffusion entry leaves the complete stationary support')
    return {'memory_dimension':12,'stationary_support':list(support),
            'original_record_weights_used':True,'group_inverse_exact':True,
            'potential_family_matrix':[[str(x) for x in row] for row in potential.tolist()],
            'potential_placement':'I2 tensor A',
            'commutes_with_actual_C3':True,
            'cyclic_orientation_coefficient':str(cyclic_coefficient),
            'cyclic_orientation_identity':'A=i*kappa*(R-R^-1); kappa=603*(1+sqrt(5))/6666724',
            'coherent_intertwining_frame_identity':'O^T*A*O=i*kappa*(P-P^T) whenever R*O=O*P',
            'orientation_generator_supplies_scalar_phase_without_control':False,
            'diffusion_positive_semidefinite':True,'diffusion_nonzero':True,
            'principal_minors_checked':len(principal),
            'corrected_reduced_generator':'G0+h*G1; G1 includes ordering and diffusion',
            'fixed_gain_fixed_t_ordering_survives_h_to_zero':False,
            'retained_potential_without_diffusion_is_same_source':False,
            'second_order_complete_archive_comparison_claimed':False,
            'owner':'model/particle_memory_response/PROOF.md'}


def mcpr_symbolic_scope():
    """Solve the stipulated affine phase law; no stored endpoint is consumed."""
    u,alpha_p,chi,zeta=s.symbols('u alpha_p chi zeta',real=True)
    chi_image=u*u/512+u*chi/77
    zeta_image=alpha_p*alpha_p/21+(u/27+u*u/135)*zeta
    chi_value=s.solve(chi-chi_image,chi)[0]
    zeta_value=s.solve(zeta-zeta_image,zeta)[0]
    if s.cancel(chi_value-(u*u/512)/(1-u/77))!=0 or s.cancel(
            zeta_value-(alpha_p*alpha_p/21)/(1-u/27-u*u/135))!=0:
        raise ArithmeticError('the declared MCPR fixed-point identity differs')
    rho=s.symbols('rho',positive=True)
    # Three C3 cosines have zero sum and squared sum 3/2.
    quotient=s.cancel((3+s.Rational(3,2)*rho*rho)/9)
    if s.simplify(quotient.subs(rho,s.sqrt(2)*s.exp(-chi))-(1+s.exp(-2*chi))/3)!=0:
        raise ArithmeticError('the declared amplitude/positive-root quotient differs')
    return {'law_owner':'code/particles/leptons/derive_charged_mcpr_completion_conditional.py',
            'phase_law':'2/9 + ((P-phi)^2/(21*pi))/(1-u/27-u^2/135)',
            'chi_law':'(u^2/512)/(1-u/77)',
            'amplitude_law':'sqrt(2)*exp(-chi)',
            'positive_root_quotient':'(1+exp(-2*chi))/3',
            'positive_root_quotient_additional_premise':'all three response roots r_j are nonnegative',
            'exact_affine_fixed_points_checked':True,
            'sign_premises':['P real','u>0','1-u/77>0','1-u/27-u^2/135>0'],
            'phase_at_least_two_ninths_under_sign_premises':True,
            'phase_equality_condition':'P=phi',
            'positive_root_quotient_strictly_below_two_thirds_for_positive_u':True,
            'native_memory_realization_owner':'model/mcpr_native_memory_bridge/PROOF.md',
            'pixel_closure_owner':'model/pixel_mass_closure/PROOF.md',
            'configured_coordinates_and_architecture_determine_phase':True,
            'architecture_historically_target_informed':True,
            'primitive_repair_selects_phase_offset_proved':False,
            'equals_exact_balanced_Koide_action':False,
            'retained_coordinate_receipts_consumed':False,
            'parent_memory_or_pixel_closure_replayed':False}


def certificate():
    default=coherent.exact_source()
    a,b=s.sqrt(coherent.WEIGHTS[0]),s.sqrt(coherent.WEIGHTS[1])
    if default['z'] != 1 or not _zero(default['C_exclusive']-a*s.eye(3)-b*(coherent.P+coherent.P**2)):
        raise ArithmeticError('the fixed coherent default phase changed')
    return {'schema':'oph.charged_lepton_response.source_phase_selection.v1',
            'fixed_coherent_program_default':{'z':'1','theta':'0',
                'action_ratios_to_smallest':['1','1','34+24*sqrt(2)'],
                'readable_program_is_ambiguous':False},
            'source_halfturns':reflection_identity(),
            'complete_base_memory_duties':fixed_memory_control(),
            'ordered_feedback':ordered_feedback_grade(),
            'actual_all_lag_memory':all_lag_control(),
            'conditional_MCPR_phase_law':mcpr_symbolic_scope(),
            'phase_from_unchanged_C_E_memory_program_derived':False,
            'changed_coherent_control_is_same_complete_program':False,
            'geometric_coherent_line_phase_owner':'model/common_source/native_mass_holonomy.py',
            'geometric_single_face_equal_link_choices':['pi/30','pi/15'],
            'geometry_loop_lift_to_charged_phase_is_automatic':False,
            'empirical_inputs_used':False,
            'universal_phase_selection_impossibility_claimed':False}
