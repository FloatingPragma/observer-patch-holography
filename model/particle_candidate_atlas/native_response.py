"""Complete C-duty finite-step response and its source-scoped obstruction."""
from functools import lru_cache
import json

import sympy as s

from . import memory_diffusion as source


def _exact(value):
    if isinstance(value, (bool, float)):
        raise ValueError('exact rational inputs are required')
    result = s.sympify(value)
    if result.is_Rational is not True:
        raise ValueError('exact rational inputs are required')
    return result


def _clean(matrix):
    return s.ImmutableMatrix(s.Matrix(matrix).applyfunc(
        lambda x: s.factor(s.radsimp(s.expand(x)))))


def _entries(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def c_duty_weights(p, q):
    """Original disabled/port2/port3/port10 probabilities in stationarity."""
    p, q = _exact(p), _exact(q)
    if not 0 < p < 1 or not 0 < q < 1:
        raise ValueError('strictly positive reset and enabled branches required')
    a, b = p*(1-q), p*q
    denominator = 1-a*a+a*b
    values = (1-p, ((1-a)*(b+a*a)+a*b*b)/denominator,
              a*(1-a)/denominator, a*b/denominator)
    if sum(values) != 1 or any(x <= 0 for x in values):
        raise ArithmeticError('the complete C-duty weights failed')
    return values


def purity_floor(p, q, tangent):
    """Lower bound for pulse angle theta=2 atan(tangent), with all records."""
    tangent = _exact(tangent)
    if tangent == 0:
        raise ValueError('the nontrivial-clock theorem excludes zero angle')
    weights = c_duty_weights(p, q)
    return s.Rational(32,25)*weights[0]*weights[2]*tangent*tangent/(1+tangent*tangent)


def finite_control(p=s.Rational(1,12), q=s.Rational(1,6), tangent=s.Rational(1,10)):
    """Exact source projectors, four complete coarse records, no postselection."""
    p,q,tangent = _exact(p),_exact(q),_exact(tangent)
    floor = purity_floor(p,q,tangent)
    row = source.exact_response('C')
    projectors = tuple(row['projectors'][i] for i in (2,3,10))
    if s.simplify(s.trace(projectors[0]*projectors[1])) != s.Rational(1,5):
        raise ArithmeticError('the literal native port overlap changed')
    weights = c_duty_weights(p,q)
    phase = (1-s.I*tangent)**2/(1+tangent*tangent)
    pulses = (s.eye(3),)+tuple(_clean(s.eye(3)+(phase-1)*a) for a in projectors)
    gram = _clean(sum((w*u.H*u for w,u in zip(weights,pulses,strict=True)),s.zeros(3)))
    if gram != s.eye(3):
        raise ArithmeticError('complete record effects do not sum to identity')
    initial = projectors[0]
    histories = tuple(_clean(w*u*initial*u.H)
                      for w,u in zip(weights,pulses,strict=True))
    output = _clean(sum(histories,s.zeros(3)))
    loss = s.simplify(1-s.trace(output*output))
    if s.trace(output) != 1 or s.simplify(loss-floor) < 0 or floor <= 0:
        raise ArithmeticError('the finite source purity witness failed')
    return {'p':str(p),'q':str(q),'pulse_tan_half_angle':str(tangent),
            'pulse_angle':'2 atan('+str(tangent)+')',
            'original_weights':list(map(str,weights)),
            'coarse_records':['disabled','port2','port3','port10'],
            'initial_family_density':_entries(initial),
            'subnormalized_history_densities':[_entries(x) for x in histories],
            'unread_family_density':_entries(output),
            'original_total_probability':'1',
            'actual_purity_loss':str(loss),
            'proved_purity_loss_floor':str(floor),
            'unhalved_distance_to_every_unitary_channel_floor':str(floor),
            'fine_records_discarded_by_conditioning':False,
            'working_register':'actual X-plus sector of the six-row logical Nambu/family service',
            'physical_particle_vacuum_assigned':False}


def all_lag_parts():
    """Actual native C-duty mean and ordering potential; immutable matrices."""
    row = source.exact_response('C'); projectors=row['projectors']; k=row['K']
    mean = _clean(sum((row['weight'][i]*projectors[i] for i in range(12)),s.zeros(3)))
    potential = _clean(-s.I*sum((k[i,j]*(projectors[i]*projectors[j]-projectors[j]*projectors[i])
                                for i in range(12) for j in range(12) if k[i,j]!=0),s.zeros(3))/2)
    if potential == s.zeros(3) or _clean(potential-potential.H) != s.zeros(3):
        raise ArithmeticError('the actual C-duty ordering potential changed')
    return {'mean':mean,'potential':potential,
            'trace_potential_squared':s.simplify(s.trace(potential*potential))}


def four_mode_pairing(tangent):
    """Exact canonical CAR control for the actual rank-one Higgs pairing."""
    tangent=_exact(tangent)
    co=(1-tangent*tangent)/(1+tangent*tangent);si=2*tangent/(1+tangent*tangent)
    creators=[]
    for mode in range(4):
        op=s.zeros(16)
        for mask in range(16):
            if not mask&(1<<mode):
                op[mask|(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
        creators.append(op)
    pairs=[creators[0]*creators[1],-creators[2]*creators[3]]
    generators=[a+a.H for a in pairs]
    if generators[0]*generators[1] != generators[1]*generators[0]:
        raise ArithmeticError('disjoint even pair generators do not commute')
    pulses=[s.eye(16)+(co-1)*a*a-s.I*si*a for a in generators]
    u=_clean(pulses[0]*pulses[1])
    if _clean(u.H*u)!=s.eye(16) or s.simplify(u[0,0]-co*co)!=0:
        raise ArithmeticError('the complete four-mode pairing control failed')
    return {'unitary':u,'generator':s.ImmutableMatrix(sum(generators,s.zeros(16))),
            'cos_angle':co,'sin_angle':si,'vacuum_amplitude':co*co}


def full_car_vacuum_control(p=s.Rational(1,12),q=s.Rational(1,6),tangent=s.Rational(1,10)):
    """Vacuum witness for the configured full-Fock C-duty kick, not a vacuum selection."""
    purity_floor(p,q,tangent)  # Strict source and nontrivial rational-clock domain.
    weights=c_duty_weights(p,q);control=four_mode_pairing(tangent)
    squared_overlap=control['cos_angle']**4
    floor=2*weights[0]*weights[2]*(1-squared_overlap)
    if floor<=0:raise ArithmeticError('full-CAR vacuum purity floor failed')
    return {'original_weights':list(map(str,weights)),
            'pulse_tan_half_angle':str(_exact(tangent)),
            'canonical_CAR_modes':4,'complete_Fock_input_columns':16,
            'pair_generator_spectrum':{str(k):v for k,v in control['generator'].eigenvals().items()},
            'vacuum_return_amplitude':str(control['vacuum_amplitude']),
            'disabled_P3_overlap_squared':str(squared_overlap),
            'full_CAR_purity_loss_floor':str(floor),
            'unhalved_distance_to_every_unitary_channel_floor':str(floor),
            'other_records_retained':True,
            'full_ninety_mode_four_record_density_replayed':False,
            'physical_vacuum_selection_claimed':False,
            'scope':'pure feedback kick; common free unitaries before/after preserve the witness; no arbitrary simultaneous H0 claim'}


def hamiltonian_trace_control(momentum=s.Rational(1,3), step=s.Rational(1,11)):
    """One-helicity canonical kinetic slice; compulsory diffusion is separate."""
    momentum,step=_exact(momentum),_exact(step)
    if momentum == 0 or step <= 0:
        raise ValueError('nonzero momentum and positive step required')
    parts=all_lag_parts(); y,a=parts['mean'],parts['potential']
    x=s.Matrix([[0,1],[1,0]]);z=s.diag(1,-1)
    h0=_clean(s.kronecker_product(x,y)+step*s.kronecker_product(s.eye(2),a))
    h=_clean(momentum*s.kronecker_product(z,s.eye(3))+h0)
    h2=_clean(h*h);h02=_clean(h0*h0)
    defect=s.simplify(s.trace(h2*h2)-6*momentum**4
                      -2*momentum**2*s.trace(h02)-s.trace(h02*h02))
    expected=8*step**2*momentum**2*parts['trace_potential_squared']
    if s.simplify(defect-expected)!=0 or expected<=0:
        raise ArithmeticError('the source kinetic/mass trace identity failed')
    return {'momentum':str(momentum),'step':str(step),'matrix_dimension':6,
            'mean_Y':_entries(y),'ordering_A':_entries(a),
            'trace_A_squared':str(parts['trace_potential_squared']),
            'same_speed_mass_shell_trace_defect':str(defect),
            'identity':'8 h^2 k^2 Tr(A^2)',
            'ordering_only_dynamics_is_the_complete_source':False,
            'physical_Lorentz_theory_excluded':False}


@lru_cache(maxsize=1)
def _certificate_json():
    value={'native_finite_step':finite_control(),
           'declared_rate_control':finite_control(s.Rational(1,24),s.Rational(1,12),s.Rational(1,20)),
           'full_CAR_vacuum_witness':full_car_vacuum_control(),
           'actual_C_duty_all_lag_hamiltonian':hamiltonian_trace_control(),
           'general_scope':'0<p,q<1; pulse angle outside 2*pi*Z; stationary memory and no postselection',
           'finite_unread_working_register_exact_unitary_replacement_excluded':True,
           'unitary_dilation_or_complete_archive_excluded':False,
           'weak_limit_mean':'H0+g X tensor Z_C, in the specified protocol and pulse clock',
           'kinetic_whitening_generated_by_that_weak_limit':False,
           'balanced_to_P_response_derived':False,
           'physical_selfenergy_computed':False,
           'physical_mass_assigned':None}
    return json.dumps(value,sort_keys=True,allow_nan=False)


def certificate():
    return json.loads(_certificate_json())
