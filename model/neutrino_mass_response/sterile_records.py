"""Complete active/sterile and conjugate records on the coupled control ground."""
from functools import lru_cache
import sympy as s

from model.coherent_defect_wave.chiral_joint.family_cell.quotient import add, divide, evaluate, mul, outward
from . import sterile as m


def encoded(value):
    return [str(s.Rational(x)) for x in value.to_list()] if hasattr(value, 'to_list') else [str(s.Rational(value))]


def field_polynomial(value):
    return s.Poly.from_list(list(map(s.Rational, encoded(value))), s.Symbol('e'))


def field_eval(poly, field, energy):
    value = field.zero
    for coefficient in poly.all_coeffs():
        value = value*energy+field.convert(coefficient)
    return value


@lru_cache(maxsize=2)
def ground_block(majorana=False):
    data = m.physical_control(majorana)
    h, metric = data['H'], data['gram']
    ground, rank = data['energy_intervals_with_ranks'][0]
    assert rank == 1
    factor = next(p for p, multiplicity in data['characteristic'].factor_list()[1]
        if (p.eval(ground[0]) == 0 if ground[0] == ground[1] else p.count_roots(*ground) == 1))
    annihilator = s.zeros(h.rows)
    for coefficient in factor.all_coeffs():
        annihilator = h*annihilator+coefficient*s.eye(h.rows)
    basis = s.Matrix.hstack(*annihilator.nullspace())
    assert basis.cols == factor.degree()
    gram = basis.H*metric*basis
    left = gram.inv()*basis.H*metric
    source_h = s.ImmutableSparseMatrix(left*h*basis)
    assert basis*source_h == h*basis and gram*source_h == source_h.H*gram
    if factor.degree() == 1:
        field = s.QQ
        energy = field.convert(ground[0])
    else:
        field = s.QQ.algebraic_field(s.CRootOf(factor.as_expr(), 0))
        energy = field.unit
    coefficients = []
    x = factor.gen
    for i in range(factor.degree()):
        polynomial = s.Poly(sum(factor.nth(j)*x**(j-i-1) for j in range(i+1, factor.degree()+1)), x)
        coefficients.append(field_eval(polynomial, field, energy)/field_eval(factor.diff(), field, energy))
    powers, value = [], s.eye(source_h.rows)
    for coefficient in coefficients:
        powers.append(value)
        value = source_h*value
    trace = sum((coefficient*field.convert(s.trace(power)) for coefficient, power in zip(coefficients, powers, strict=True)), field.zero)
    assert trace == field.one
    return {'control': data, 'basis': basis, 'gram': gram, 'H': source_h,
            'ground_interval': tuple(map(s.Rational, ground)), 'ground_factor': factor,
            'field': field, 'energy': energy, 'coefficients': coefficients, 'powers': powers}


def probe_matrix(kind, spin=0):
    if kind == 'active_composite':
        polynomial = m.lift(m.car.neutral_composite(spin))
    elif kind == 'sterile':
        polynomial = m.car.monomial(((m.STERILE[spin], False),))
    else:
        raise ValueError('active composite or sterile probe required')
    data = m.physical_basis()
    raw = m.matrix(m.restricted(polynomial), data['keys'], data['keys'])
    result = data['left']*raw*data['basis']
    assert data['basis']*result == raw*data['basis']
    return s.ImmutableSparseMatrix(result)


@lru_cache(maxsize=4)
def measure(kind, majorana=False):
    ground = ground_block(majorana)
    data, field, energy = ground['control'], ground['field'], ground['energy']
    operator = probe_matrix(kind)
    image = operator*ground['basis']
    adjoint = ground['gram'].inv()*image.H*data['gram']
    q = data['minimal']
    moments, value = [], image
    for _ in range(q.degree()):
        pullback = adjoint*value
        moments.append(sum((coefficient*field.convert(s.trace(power*pullback))
                            for coefficient, power in zip(ground['coefficients'], ground['powers'], strict=True)), field.zero))
        value = data['H']*value
    z = s.Symbol('z')
    entries = {}
    for i, moment in enumerate(moments):
        for j in range(i+1, q.degree()+1):
            entries[j-i-1] = entries.get(j-i-1, field.zero)+field.convert(q.nth(j))*moment
    numerator = s.Poly.from_dict({(j,): value for j, value in entries.items()}, z, domain=field)
    denominator = s.Poly.from_dict({(j,): field.convert(q.nth(j)) for j in range(q.degree()+1)}, z, domain=field)
    assert numerator.rep.nth(q.degree()-1) == moments[0]
    common = s.gcd(numerator, denominator)
    reduced_numerator = numerator.exquo(common)
    reduced_denominator = s.Poly(denominator.exquo(common).as_expr(), z, domain=s.QQ)
    invisible = s.Poly(common.as_expr(), z, domain=s.QQ)
    rows = []
    for bounds, rank in data['energy_intervals_with_ranks']:
        bounds = tuple(map(s.Rational, bounds))
        zero = (invisible.eval(bounds[0]) == 0 if bounds[0] == bounds[1] else invisible.count_roots(*bounds) == 1)
        if zero:
            weight = (s.Integer(0), s.Integer(0))
        else:
            value = (s.Integer(0), s.Integer(0))
            for power in range(reduced_numerator.degree(), -1, -1):
                value = add(mul(value, bounds), evaluate(field_polynomial(reduced_numerator.rep.nth(power)), ground['ground_interval']))
            weight = divide(value, evaluate(reduced_denominator.diff(), bounds))
            assert weight[0] > 0
        original = (weight[0]/1280, weight[1]/1280)  # preparation1/80, setting1/4, Kraus1/4
        rows.append({'source_energy_interval': outward(bounds, 60), 'full_physical_energy_rank': rank,
            'energy_transfer_interval': outward((bounds[0]-ground['ground_interval'][1], bounds[1]-ground['ground_interval'][0]), 60),
            'ground_conditioned_joint_probability_interval': outward((weight[0]/16, weight[1]/16), 80),
            'complete_preparation_joint_probability_interval': outward(original, 88),
            'exact_zero_residue': bool(zero)})
        if weight[0] == weight[1]:
            rows[-1]['exact_unscaled_probe_residue'] = str(weight[0])
            rows[-1]['exact_ground_conditioned_joint_probability'] = str(weight[0]/16)
            rows[-1]['exact_complete_preparation_joint_probability'] = str(weight[0]/1280)
    norm_interval = evaluate(field_polynomial(moments[0]), ground['ground_interval'])
    assert 0 < norm_interval[0] <= norm_interval[1] <= 2
    failure = field.convert(s.Rational(1, 4))-moments[0]/16
    failure_interval = evaluate(field_polynomial(failure), ground['ground_interval'])
    assert failure_interval[0] > 0
    return {'probe': kind, 'four_settings': [{'operation': operation, 'spin': spin,
            'setting_probability': '1/4', 'successful_Kraus': ('J' if kind == 'active_composite' else 'n^c')+('^dagger/2' if operation == 'addition' else '/2'),
            'complete_energy_records': rows,
            'joint_ground_conditioned_failure_probability_interval': outward(failure_interval, 80),
            'joint_complete_preparation_failure_probability_interval': outward((failure_interval[0]/80, failure_interval[1]/80), 88),
            'outside_complete_physical_target_probability': '0',
            'failure_state_and_fermionic_quantum_reference_retained': True}
            for operation in ('removal', 'addition') for spin in (0, 1)],
        'equal_spin_measures_proof': 'full spin Ward identities and unique physical spin-singlet ground',
        'equal_addition_removal_measures_proof': 'exact full-carrier antiunitary maps each probe to i times its adjoint and preserves the unique ground',
        'exact_norm_squared_coefficients_in_ground_energy': encoded(moments[0]),
        'exact_first_moment_coefficients_in_ground_energy': encoded(moments[1]),
        'exact_full_target_annihilator': str(q.as_expr()),
        'exact_visible_support_annihilator': str(reduced_denominator.as_expr()),
        'exact_zero_residue_factor': str(invisible.as_expr()),
        'exact_resolvent_numerator_coefficients_in_ground_energy': [encoded(numerator.rep.nth(i)) for i in range(numerator.degree(), -1, -1)],
        'complete_record_count_per_ground_prepared_experiment': 4*(len(rows)+2),
        'physical_rest_masses_assigned_to_these_energy_lines': False}


@lru_cache(maxsize=2)
def calculate(majorana=False):
    ground = ground_block(majorana)
    data = ground['control']
    records = [{'energy_interval': list(map(str, bounds)), 'rank': rank,
                'original_preparation_probability': str(s.Rational(rank, 80))}
               for bounds, rank in data['energy_intervals_with_ranks']]
    assert sum(row['rank'] for row in records) == 80
    probes = {kind: measure(kind, majorana) for kind in ('active_composite', 'sterile')}
    dirac_selfenergy = None
    if not majorana:
        z = s.Symbol('z')
        values = probes['sterile']['four_settings'][0]['complete_energy_records']
        weights = {s.Rational(row['energy_transfer_interval'][0]): s.Rational(row['exact_unscaled_probe_residue'])
                   for row in values if not row['exact_zero_residue']}
        assert weights == {s.Integer(1): s.Rational(35, 72), s.Integer(5): s.Rational(1, 72)}
        for row in probes['active_composite']['four_settings'][0]['complete_energy_records']:
            if not row['exact_zero_residue']:
                gap = s.Rational(row['energy_transfer_interval'][0])
                assert s.Rational(row['exact_unscaled_probe_residue']) == gap**2*weights[gap]
        retarded = s.cancel(sum(weight*(1/(z-gap)+1/(z+gap)) for gap, weight in weights.items()))
        selfenergy = s.cancel(z-1/retarded)
        assert s.cancel(selfenergy-s.Rational(75, 73)/z-s.Rational(140, 219)*z/(z*z-s.Rational(73, 3))) == 0
        dirac_selfenergy = {'complete_CAR_retarded_propagator': str(retarded),
            'complete_CAR_normalization': '1',
            'exact_Feshbach_selfenergy': str(selfenergy),
            'exact_positive_residue_expansion': '75/(73 z)+(140/219) z/(z^2-73/3)',
            'full_cyclic_dimension': 4, 'full_CAR_Feshbach_complement_dimension': 3,
            'positive_and_negative_poles': ['-5', '-1', '1', '5'],
            'normalized_removal_resolvent': '(z-44/9)/((z-1)(z-5))',
            'normalized_removal_selfenergy': '(35/81)/(z-44/9)',
            'normalized_removal_diagonal_energy': '10/9',
            'low_frequency_mass_squared_like_residue': '75/73',
            'mass_squared_like_residue_identified_with_relativistic_mass': False}
    return {'source': 'full8CAR/5scalar lepton and coupled sterile control, yu=yd=0,ye=y_nu=1',
        'optional_sterile_Majorana_primitive': int(majorana),
        'full_control_carrier_dimension': 1280, 'full_physical_Gauss_dimension': 80,
        'complete_ground_compared_across_all_physical_sectors': True,
        'ground_rank': 1, 'ground_energy_interval': list(map(str, ground['ground_interval'])),
        'ground_minimal_polynomial': str(ground['ground_factor'].as_expr()),
        'full_physical_characteristic_factorization': str(s.factor(data['characteristic'].as_expr())),
        'complete_preparation': 'vec(P_physical)/sqrt80 and full energy measurement; every nonground outcome retained',
        'ground_original_preparation_probability': '1/80',
        'complete_original_preparation_energy_records': records,
        'terminal_preparation_records_before_ground_probe_refinement': records[1:],
        'record_tree_rule': 'the ground preparation record is refined by the probe children; it is not counted a second time',
        'ground_reference': 'unique physical ground with its retained conjugate preparation reference',
        'clock': 'exp(-it H_control), hbar1',
        'complete_neutral_probe_records': probes,
        'unit_Dirac_control_complete_normal_selfenergy': dirac_selfenergy,
        'same_scalar_ground_or_spectrum_transferred_to_full32_all_unit_source': False,
        'intrinsic_tau_neutrino_mass_or_measured_flavor_mixing_derived': False}
