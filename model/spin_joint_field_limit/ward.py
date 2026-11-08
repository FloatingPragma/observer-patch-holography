"""Source Ward identities, static force slacks and original spectral weights."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import sympy as s

from model.native_ordered_ground import source
from .force import _tensor, _simple, _zero, _nonnegative, _quadratic_interval


def _rational(value, name, positive=False):
    if type(value) not in (int, F):
        raise TypeError(name+' must be an exact rational')
    out = F(value)
    if out < 0 or (positive and out == 0):
        raise ValueError(name+' has invalid sign')
    return out


def interval_slack(energy, cross, susceptibility):
    """Conservative interval with a certified strictly positive response floor."""
    def interval(pair, name):
        if len(pair) != 2 or any(type(x) not in (int, F) for x in pair):
            raise TypeError(name+' requires two exact rational endpoints')
        lo, hi = map(F, pair)
        if lo > hi:
            raise ValueError(name+' interval is reversed')
        return lo, hi
    h0, h1 = interval(energy, 'first moment')
    c0, c1 = interval(cross, 'cross moment')
    x0, x1 = interval(susceptibility, 'susceptibility')
    if min(h0, x0) < 0 or x0 == 0:
        raise ValueError('a positive susceptibility floor is required')
    c2_min = F(0) if c0 <= 0 <= c1 else min(c0*c0, c1*c1)
    c2_max = max(c0*c0, c1*c1)
    low, high = max(F(0), h0-c2_max/x0), h1-c2_min/x1
    if high < low:
        raise ValueError('intervals contradict the positive source Ward slack')
    coefficients = (c0/x0, c0/x1, c1/x0, c1/x1)
    return {'force_form_slack_interval': [low, high],
            'source_coefficient_interval': [min(coefficients), max(coefficients)],
            'susceptibility_floor': x0,
            'empirical_coefficient_fit_used': False}


def tube_weight(cutoff, u, b, delta_q, delta_p, radius,
                original_norm_squared, preparation_normalization):
    """Rational form of (10)-(11), using (x+y)^2 <= 2x^2+2y^2."""
    m = _rational(cutoff, 'energy cutoff', True)
    dq = _rational(delta_q, 'Q slack')
    dp = _rational(delta_p, 'P slack')
    eps = _rational(radius, 'energy tube radius', True)
    nq = _rational(original_norm_squared, 'filtered source norm')
    lam = _rational(preparation_normalization, 'preparation normalization', True)
    if type(u) not in (int, F) or type(b) not in (int, F):
        raise TypeError('source coefficients must be exact rationals')
    u, b = F(u), F(b)
    if u*b <= 0 or nq > lam*lam:
        raise ValueError('positive frequency and a valid preparation norm are required')
    outside = 2*m*(m*m*dq+u*u*dp)/(u*b*eps*eps)
    return {'frequency_squared': u*b,
            'original_filtered_preparation_weight': nq/(lam*lam),
            'outside_tube_source_norm_squared_upper': outside,
            'original_tube_probability_lower': max(F(0), nq-outside)/(lam*lam),
            'success_probability_divided_out': False,
            'source_slack_vanishing_assumed_or_proved_here': False}


def moment_error(normalization, score_bound, half_trace_error, setting_probability=F(1)):
    lam = _rational(normalization, 'preparation normalization', True)
    cap = _rational(score_bound, 'score bound')
    eta = _rational(half_trace_error, 'complete output error')
    pi = _rational(setting_probability, 'setting probability', True)
    if pi > 1:
        raise ValueError('setting probability exceeds one')
    return 2*lam*lam*cap*eta/pi


def uniform_variance_slack(first, inverse, third_inverse, field):
    """The uniform Ward pair's exact residual from original energy moments."""
    h = _rational(first, 'first moment', True)
    chi = _rational(inverse, 'inverse moment', True)
    chi3 = _rational(third_inverse, 'third inverse moment', True)
    b = _rational(field, 'selection field', True)
    delta = h-chi*chi/chi3
    if delta < 0:
        raise ValueError('moments contradict the positive uniform Ward variance')
    return {'Q_force_form_slack': delta, 'P_force_form_slack': F(0),
            'source_cross_moment': b*chi,
            'P_inverse_moment': b*b*chi3, 'P_first_moment': b*b*chi,
            'source_coefficients_u_b': [chi/(b*chi3), b],
            'weighted_energy_squared_mean': chi/chi3,
            'weighted_energy_squared_variance': delta/chi3,
            'auxiliary_weight': 'E^-3 dnu/M_-3; original event probabilities unchanged',
            'continuum_variance_vanishing_inferred': False}


def _expect(ground, operator):
    return _simple(s.trace(ground*operator))


def _complete_records(ground, projectors, q, momentum, lam_q, lam_p, kind='energy'):
    """Four declared settings, each with failure and every source energy label."""
    settings = (('Q', q, lam_q), ('P', momentum, lam_p),
                ('Q_minus_iP', q-s.I*momentum, lam_q+lam_p),
                ('Q_plus_iP', q+s.I*momentum, lam_q+lam_p))
    rows = []
    total_effect = s.zeros(4)
    totals = {}
    outcome_moments = {str(label): {} for label, _ in projectors}
    for label, op, normalization in settings:
        effect = (op.H*op/normalization**2).applyfunc(_simple)
        for candidate in (effect, s.eye(4)-effect):
            for size in range(1, 5):
                for indices in combinations(range(4), size):
                    minor = _simple(candidate.extract(indices, indices).det())
                    if not minor.is_Rational or minor < 0:
                        raise ArithmeticError('source creation is not a contraction')
        failure = (s.eye(4)-effect)/4
        total_effect += failure
        failure_probability = _expect(ground, failure)
        if not _nonnegative(failure_probability):
            raise ArithmeticError('negative source failure probability')
        rows.append({'setting': label, 'outcome': 'preparation_failure',
                     'original_probability': str(failure_probability)})
        success = s.Integer(0)
        for energy, projector in projectors:
            full_effect = op.H*projector*op/(4*normalization**2)
            total_effect += full_effect
            probability = _expect(ground, full_effect)
            if not _nonnegative(probability):
                raise ArithmeticError('negative complete energy record')
            success += probability
            outcome_moments[str(energy)][label] = _simple(4*normalization**2*probability)
            rows.append({'setting': label, 'outcome': kind+'_'+str(energy),
                         'original_probability': str(probability)})
        totals[label] = _simple(4*normalization**2*success)
    if not _zero(total_effect-s.eye(4)):
        raise ArithmeticError('the four-setting source instrument is incomplete')
    if _simple(sum(s.sympify(row['original_probability']) for row in rows)-1) != 0:
        raise ArithmeticError('complete original probabilities do not sum to one')
    measured_c = _simple((totals['Q_minus_iP']-totals['Q_plus_iP'])/4)
    expected_c = _simple(s.im(_expect(ground, q*momentum)))
    if measured_c != expected_c:
        raise ArithmeticError('the original interference cross score has the wrong sign')
    outcome_cross = {}
    for label, projector in projectors:
        values = outcome_moments[str(label)]
        cross = _simple((values['Q_minus_iP']-values['Q_plus_iP'])/4)
        if cross != _simple(s.im(_expect(ground, q*projector*momentum))):
            raise ArithmeticError('a filtered original cross score differs from its source record')
        outcome_cross[str(label)] = str(cross)
    return {'setting_probability': '1/4', 'settings': 4, 'records': rows,
            'complete_record_count': len(rows), 'probability_sum': '1',
            'zero_probability_records_retained': True,
            'original_interference_cross_moment': str(measured_c),
            'original_cross_moment_by_retained_outcome': outcome_cross,
            'preparation_success_renormalized': False,
            'finite_effects_specified': True, 'additional_native_tape_executed': False}


@lru_cache(maxsize=1)
def finite_source():
    """Exact joint n=2 source, including a positive-slack sharp-line control."""
    p, field = s.Rational(1, 4), s.Rational(1, 32)
    integer, scale = source.matrix(2, ((0, 1),), F(1), F(1, 32))
    marshall = s.diag(1, 1, -1, -1)
    hb = marshall*s.Matrix(integer.tolist())*marshall/scale
    h0 = sum((_tensor(a*2) for a in 'XYZ'), s.zeros(4))/4-s.eye(4)/4
    if not _zero(hb-h0+field*(_tensor('XI')-_tensor('IX'))/2):
        raise ArithmeticError('actual physical source frame differs')
    e0, e1 = -(16+s.sqrt(257))/32, (-16+s.sqrt(257))/32
    if hb.eigenvals() != {e0: 1, s.Integer(0): 2, e1: 1}:
        raise ArithmeticError('full actual source spectrum differs')
    vector = (hb-e0*s.eye(4)).nullspace()[0]
    ground = (vector*vector.H/(vector.H*vector)[0]).applyfunc(_simple)
    middle = s.eye(4)-hb*(hb+s.eye(4))/field**2
    upper = (s.eye(4)-ground-middle).applyfunc(_simple)
    k = (hb-e0*s.eye(4))/p
    energy = _simple(-e0/p)
    projectors = ((s.Integer(0), ground), (energy, middle), (_simple((e1-e0)/p), upper))
    if not _zero(k*middle-energy*middle):
        raise ArithmeticError('source energy window changed')
    rows = []
    for profile in ((1, 1), (1, 2)):
        f0, f1 = map(s.Integer, profile)
        q = p*p*(f0*_tensor('ZI')-f1*_tensor('IZ'))/2
        momentum = -p*(f0*_tensor('YI')+f1*_tensor('IY'))/2
        current2 = (_tensor('ZX')-_tensor('XZ'))/4
        current3 = (_tensor('XY')-_tensor('YX'))/4
        z = p*p*(f1-f0)*current2
        torque = -p*(f0+f1)*current3
        dq_op = s.I*(k*q-q*k)
        dp_op = s.I*(k*momentum-momentum*k)
        if not _zero(dq_op-torque-field*momentum):
            raise ArithmeticError('actual Q exchange-current or torque sign differs')
        if not _zero(dp_op+z/(p*p)+q/2):
            raise ArithmeticError('actual P exchange-current or torque sign differs')
        if not _zero((s.eye(4)-middle)*q*ground) or not _zero((s.eye(4)-middle)*momentum*ground):
            raise ArithmeticError('finite transverse source preparation left the complete middle level')
        nq, np = _expect(ground, q*q), _expect(ground, momentum*momentum)
        c = _simple(s.im(_expect(ground, q*momentum)))
        c_source = _expect(ground, p**3*(f0*f0*_tensor('XI')-f1*f1*_tensor('IX'))/4)
        if c != c_source:
            raise ArithmeticError('exact source order cross moment differs')
        hq, hp = _expect(ground, q*k*q), _expect(ground, momentum*k*momentum)
        for operator, first_moment in ((q, hq), (momentum, hp)):
            commutator = k*operator-operator*k
            double = (operator*commutator-commutator*operator)/2
            if _expect(ground, double) != first_moment:
                raise ArithmeticError('static source double commutator differs from first moment')
        cq, cp = _simple(nq/energy), _simple(np/energy)
        u, b = _simple(c/cp), _simple(c/cq)
        deltaq, deltap = _simple(hq-c*c/cp), _simple(hp-c*c/cq)
        rq, rp = s.I*k*q-u*momentum, s.I*k*momentum+b*q
        if _expect(ground, rq.H*middle*rq/energy) != deltaq:
            raise ArithmeticError('Q form slack differs from the full source residual')
        if _expect(ground, rp.H*middle*rp/energy) != deltap:
            raise ArithmeticError('P form slack differs from the full source residual')
        if profile == (1, 1):
            if deltaq != 0 or deltap != 0 or _simple(u*b-energy**2) != 0 or b != s.Rational(1, 2):
                raise ArithmeticError('the exact uniform-source Ward equations changed')
            # K_unscaled A0 = i B M3; all four rows are retained.
            a0, m3 = momentum/p, q/(p*p)
            if not _zero((hb-e0*s.eye(4))*a0*ground-s.I*field*m3*ground):
                raise ArithmeticError('finite-field uniform-spin Ward identity changed')
            chi_m = _expect(ground, m3*middle*m3)/_simple(-e0)
            selected_order = _expect(ground, (_tensor('XI')-_tensor('IX'))/2)
            if _simple(chi_m-selected_order/(2*field)) != 0:
                raise ArithmeticError('exact staggered susceptibility Ward identity differs')
            unscaled_energy = _simple(-e0)
            m1 = _expect(ground, m3*(hb-e0*s.eye(4))*m3)
            mm3 = _expect(ground, m3*middle*m3)/unscaled_energy**3
            weighted_mean = _simple(chi_m/mm3)
            uniform_slack = _simple(m1-chi_m**2/mm3)
            variance = _simple(m1/mm3-weighted_mean**2)
            if (uniform_slack != 0 or variance != 0
                    or _simple(weighted_mean-unscaled_energy**2) != 0
                    or _simple(uniform_slack-mm3*variance) != 0):
                raise ArithmeticError('actual uniform source weighted-energy variance differs')
            uniform_variance = {'clock': 'K=H_B-E_B',
                'preparation': 'M3=-S1^3+S0^3; normalization Lambda=1',
                'original_preparation_weight': str(_expect(ground, m3*m3)),
                'first_moment': str(m1), 'inverse_moment': str(chi_m),
                'third_inverse_moment': str(_simple(mm3)),
                'weighted_energy_squared_mean': str(weighted_mean),
                'weighted_energy_squared_variance': str(variance),
                'Q_force_form_slack': str(uniform_slack),
                'P_force_form_slack': '0',
                'auxiliary_weight_changes_original_event_probabilities': False}
        else:
            if _quadratic_interval(deltaq)[0] <= 0 or _quadratic_interval(deltap)[0] <= 0:
                raise ArithmeticError('positive-slack sharp-line control lost its strict margin')
        lamq, lamp = p*p*(abs(f0)+abs(f1))/2, p*(abs(f0)+abs(f1))/2
        records = _complete_records(ground, projectors, q, momentum, lamq, lamp)
        row = {'profile': list(profile), 'source_norms_squared': [str(nq), str(np)],
                     'source_first_moments': [str(hq), str(hp)],
                     'source_inverse_moments': [str(cq), str(cp)], 'order_cross_moment': str(c),
                     'source_coefficients_u_b': [str(u), str(b)],
                     'positive_form_slacks': [str(deltaq), str(deltap)],
                     'one_retained_positive_energy': str(energy),
                     'all_original_energy_weight_on_that_line': True,
                     'positive_pair_slack_excludes_a_sharp_component': False,
                     'original_Q_creation_weight': str(_simple(nq/lamq**2)),
                     'complete_preparation_energy_records': records}
        if profile == (1, 2):
            charge_projections = []
            charged = []
            charge_operator = (_tensor('XI')+_tensor('IX'))/2
            for charge in (-1, 1):
                projection = (s.eye(4)+charge*_tensor('XI'))*(s.eye(4)+charge*_tensor('IX'))/4
                if (not _zero(projection*projection-projection) or s.trace(projection) != 1
                        or not _zero(k*projection-projection*k)
                        or not _zero(charge_operator*projection-charge*projection)):
                    raise ArithmeticError('actual conserved source charge projector differs')
                nqr = _expect(ground, q*projection*q)
                npr = _expect(ground, momentum*projection*momentum)
                cr = _simple(s.im(_expect(ground, q*projection*momentum)))
                ur, br = _simple(cr*energy/npr), _simple(cr*energy/nqr)
                dqr = _simple(energy*nqr-cr*cr*energy/npr)
                dpr = _simple(energy*npr-cr*cr*energy/nqr)
                rq = s.I*k*projection*q-ur*projection*momentum
                rp = s.I*k*projection*momentum+br*projection*q
                if dqr != 0 or dpr != 0 or not _zero(rq*ground) or not _zero(rp*ground):
                    raise ArithmeticError('the actual filtered source force pair failed to close')
                if _simple(ur*br-energy**2) != 0:
                    raise ArithmeticError('filtered source frequency differs from its complete energy')
                original_q, original_p = _simple(nqr/lamq**2), _simple(npr/lamp**2)
                if min(_quadratic_interval(original_q)[0], _quadratic_interval(original_p)[0]) <= F(1, 25):
                    raise ArithmeticError('actual charged source branch lost its original probability floor')
                charged.append({'internal_source_charge': charge,
                    'filtered_cross_moment': str(cr),
                    'source_coefficients_u_b': [str(ur), str(br)],
                    'force_form_slacks': [str(dqr), str(dpr)],
                    'original_Q_branch_probability': str(original_q),
                    'original_P_branch_probability': str(original_p),
                    'original_probability_floor_for_each': '1/25',
                    'source_excitation_energy': str(energy)})
                charge_projections.append((charge, projection))
            dark = s.eye(4)-sum((proj for _, proj in charge_projections), s.zeros(4))
            if not _zero(dark-(ground+upper)):
                raise ArithmeticError('complete charge complement omits a source state')
            charge_projections.append(('complement', dark))
            row['actual_charge_resolved_force_pairs'] = charged
            row['complete_preparation_charge_records'] = _complete_records(
                ground, charge_projections, q, momentum, lamq, lamp, 'charge')
            if _simple(sum(s.sympify(x['filtered_cross_moment']) for x in charged)-c) != 0:
                raise ArithmeticError('filtered source cross moments do not recover the full moment')
        rows.append(row)
    return {'sites': 2, 'coupling': '1', 'field': '1/32', 'p': '1/4',
            'positive_clock': '(H_B-E_B)/(Jp)', 'actual_source_frame_and_full_spectrum_checked': True,
            'source_current_identities_checked_on_all_rows': True,
            'uniform_spin_Ward_susceptibility_exact': True,
            'uniform_source_weighted_energy_variance': uniform_variance,
            'profiles': rows, 'continuum_particle_or_mass_inferred': False}


def certificate():
    return {'actual_source': finite_source(),
            'rational_tube_example': tube_weight(3, 2, 2, F(1, 10000), F(1, 10000),
                                                F(1, 2), F(1, 4), 1),
            'complete_score_error_example': moment_error(2, 5, F(1, 10000), F(1, 4)),
            'susceptibility_interval_example': interval_slack((F(3, 2), F(151, 100)),
                (F(49, 100), F(1, 2)), (F(1, 5), F(21, 100))),
            'uniform_Ward_broad_measure_implication_control': {
                'source': 'abstract K=diag(1,2), q=(1,1)/sqrt(2), a=iB K^-1 q',
                'actual_native_source_or_continuum_measure': False,
                'measure': [['1', '1/2'], ['2', '1/2']],
                'calculation': uniform_variance_slack(F(3, 2), F(3, 4), F(9, 16), F(1, 3))},
            'static_Ward_data_reduce_both_form_force_defects': True,
            'filtered_cross_moment_replaced_by_full_commutator': False,
            'uniform_continuum_slack_vanishing_proved': False,
            'empirical_velocity_fit_or_physical_mass_used': False}
