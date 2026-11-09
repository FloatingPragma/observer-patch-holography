"""Vacuum-compressed full-CAR lift of the configured coherent Koide source."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s
from sympy.polys.domains import QQ_I
from sympy.polys.matrices import DomainMatrix

from model.standard_model_source_map import koide_source as source
from model.common_source import native_matter_attachment as matter


def _clean(a):
    return s.ImmutableSparseMatrix(a.applyfunc(s.simplify))


def _rational(value):
    if type(value) not in (int, F, str):
        raise ValueError('an exact rational is required')
    return F(value.numerator, value.denominator) if type(value) is F else F(value)


def archive_bound(steps, total_time, *, higgs_norm=1, round_dilation_error=0):
    """Unhalved complete flagged-output error, with all quantum archives kept."""
    if type(steps) is not int or steps <= 0:
        raise ValueError('positive integer rounds required')
    t, v, error = map(_rational, (total_time, higgs_norm, round_dilation_error))
    if min(t, v, error) < 0:
        raise ValueError('nonnegative time, amplitude and certified error required')
    norm = 6*v
    leakage = norm*norm*t*t/steps
    return {'pairing_operator_norm': str(norm),
            'ideal_any_nonvacuum_probability_upper': str(min(F(1), leakage)),
            'implemented_any_nonvacuum_probability_upper': str(min(F(1), leakage+steps*error)),
            'complete_unhalved_output_error_upper': str(min(F(2), 3*leakage+2*steps*error)),
            'all_nonvacuum_quantum_archives_retained': True,
            'vacuum_record_postselection': False,
            'arbitrary_physical_occupation_and_reference_inputs': True}


@lru_cache(maxsize=2)
def _source_lift(phase=1):
    """Immutable exact one-particle source factors; no native-tape replay."""
    phase = source._phase(phase)
    prep, fourier, padded = source.analyzer_factors(phase)
    blocks = tuple(s.ImmutableMatrix(a[:3, :3]) for a in padded)
    for a in blocks:
        if _clean(a.H*a) != s.eye(3):
            raise ArithmeticError('source controlled family rotation is not unitary')
    # The first output flag row of U=F controlled(R) prep, with all inputs.
    row = _clean(s.Matrix.hstack(*[
        sum((fourier[0, j]*prep[j, flag]*blocks[j] for j in range(8)), s.zeros(3))
        for flag in range(8)]))
    if _clean(row*row.H) != s.eye(3):
        raise ArithmeticError('the complete analyzer row is not a coisometry')
    y = _clean(row[:, :3].H*row[:, :3])
    expected = source.exact_source(phase)['action_coefficient']
    if _clean(y-expected) != s.zeros(3):
        raise ArithmeticError('physical vacuum compression differs from the source coefficient')
    return row, y


def source_lift(phase=1):
    row, y = _source_lift(source._phase(phase))
    return {'analyzer_flag_zero_row': row, 'physical_family_coefficient': y}


def higgs_covariance():
    """Literal lepton tensor supports and all twelve native current blocks."""
    rep = matter.representation()
    selected = list(rep['selected'])
    weak = rep['operators']['weak_C2'].extract(selected, selected)
    candidates = []
    for line in matter.invariant_higgs_lines():
        support = {i for component in line['components'] for (i, j), v in component.todok().items() if v != 0}
        left_ids = [i for i in support if weak[i, i] != 0]
        right_ids = [i for i in support if weak[i, i] == 0]
        if len(right_ids) == 1:
            candidates.append((line, s.diag(*(int(i in left_ids) for i in range(15))),
                               s.diag(*(int(i in right_ids) for i in range(15)))))
    if len(candidates) != 1:
        raise ArithmeticError('one charged-lepton Higgs line required')
    row, left, right = candidates[0]
    if left*right != s.zeros(15) or s.trace(left) != 2 or s.trace(right) != 1:
        raise ArithmeticError('lepton doublet/singlet support changed')
    currents = rep['matter_currents']
    for g in currents:
        if not matter.zero(g*left-left*g) or not matter.zero(g*right-right*g):
            raise ArithmeticError('family-flag analyzer would break a native current')
    for component in row['components']:
        d = left*component*right
        if not matter.zero(d.H*d-right) or not matter.zero(d*d.H*d-d):
            raise ArithmeticError('Higgs component is not the stated unit partial isometry')
    return {'left_internal_rank': 2, 'right_internal_rank': 1,
            'both_Higgs_components_checked': 2, 'source_current_support_checks': 24,
            'neutral_directed_Higgs_rank': 1,
            'left_frame': 'U on the entire lepton doublet',
            'right_frame': 'complex conjugate U on the singlet',
            'fixed_neutral_Higgs_commutes_all_twelve_currents': False,
            'full_current_covariance_transforms_Higgs_too': True}


def _creator(modes, index):
    out = {}
    for mask in range(1 << modes):
        if not mask >> index & 1:
            out[mask | (1 << index), mask] = (-1)**((mask & ((1 << index)-1)).bit_count())
    return s.SparseMatrix(1 << modes, 1 << modes, out)


def _pairing(w):
    n = w.rows
    cs = tuple(_creator(2*n, j) for j in range(2*n))
    g = s.zeros(1 << (2*n), cls=s.SparseMatrix)
    for i in range(n):
        for j in range(n):
            if w[i, j] != 0:
                g += w[i, j]*(cs[i]*cs[n+j])
    return _clean(g+g.H)


@lru_cache(maxsize=1)
def _finite_control():
    """Eight-mode exact control with a genuinely complex off-diagonal mean."""
    c, d = s.Rational(3, 5), s.Rational(4, 5)
    prep = s.Matrix([[c, -d], [d, c]])
    family = s.Matrix([[c, s.I*d], [-d, s.I*c]])
    fourier = s.Matrix([[1, 1], [1, -1]])/s.sqrt(2)
    u = _clean(s.kronecker_product(fourier, s.eye(2))*s.diag(s.eye(2), family)*
               s.kronecker_product(prep, s.eye(2)))
    q = s.diag(1, 1, 0, 0)
    w = _clean(u.H*q*u)
    y = s.ImmutableMatrix(w[:2, :2])
    if _clean(u.H*u) != s.eye(4) or _clean(w*w) != w or s.trace(w) != 2:
        raise ArithmeticError('finite analyzer projection differs')
    delta = s.zeros(8); delta[:4, 4:] = q; delta[4:, :4] = -q
    frame = s.diag(u, u.conjugate())
    transformed = _clean(frame.H*delta*frame.conjugate())
    target = s.zeros(8); target[:4, 4:] = w; target[4:, :4] = -w.T
    if transformed != target:
        raise ArithmeticError('pairing conjugation convention differs')
    wrong = s.diag(u, u)
    if _clean(wrong.H*delta*wrong.conjugate()) == target:
        raise ArithmeticError('the wrong conjugate frame was not detected')
    g = _pairing(w)
    physical = (0, 1, 4, 5)
    addresses = tuple(sum(1 << physical[j] for j in range(4) if mask >> j & 1)
                      for mask in range(16))
    h = _pairing(y)
    if g.extract(addresses, addresses) != h:
        raise ArithmeticError('all-occupation vacuum compression failed')
    # Full 256-state polynomial identity, independent of eigenvalue fitting.
    generator = g.to_DM().convert_to(QQ_I)
    g2 = generator*generator
    identity = DomainMatrix.eye((256, 256), QQ_I)
    if not (generator*(g2-identity)*(g2-identity.scalarmul(QQ_I(4)))).is_zero_matrix:
        raise ArithmeticError('full-Fock pair spectrum polynomial failed')
    if sum(g2.diagonal()) != QQ_I(256) or sum((g2*g2).diagonal()) != QQ_I(640):
        raise ArithmeticError('complete Fock spectral multiplicities differ')
    # Rational unit-circle h=2 atan(1/64); polynomial functional calculus.
    co, si = s.Rational(4095, 4097), s.Rational(128, 4097)
    x = s.Symbol('x')
    values = {0: 1, 1: co-s.I*si, -1: co+s.I*si,
              2: (co-s.I*si)**2, -2: (co+s.I*si)**2}
    polynomial = s.Poly(s.interpolate(list(values.items()), x), x)
    minimal = s.Poly(x*(x*x-1)*(x*x-4), x)
    adjoint_polynomial = s.Poly(sum(s.conjugate(a)*x**j for (j,), a in polynomial.terms()), x)
    if not (polynomial*adjoint_polynomial-s.Poly(1, x)).rem(minimal).is_zero:
        raise ArithmeticError('clock polynomial does not preserve the exact full spectrum')
    embedding = DomainMatrix({address: {column: QQ_I.one} for column, address in enumerate(addresses)},
                            (256, 16), QQ_I)
    evolution = DomainMatrix.zeros((256, 16), QQ_I)
    for coefficient in polynomial.all_coeffs():
        evolution = generator*evolution+embedding.scalarmul(QQ_I.from_sympy(coefficient))
    evolution = s.ImmutableSparseMatrix(evolution.to_Matrix())
    # Unitarity on every row follows from the checked generator polynomial
    # and the exact modular polynomial identity, without dense multiplication.
    ancillary = (2, 3, 6, 7)
    kraus = []
    for label in range(16):
        amask = sum(1 << ancillary[j] for j in range(4) if label >> j & 1)
        # Canonical physical/ancilla Fock tensor order has a wedge shuffle.
        signs = [(-1)**sum((amask & ((1 << p)-1)).bit_count()
                           for p in physical if word >> p & 1) for word in addresses]
        rows = [word | amask for word in addresses]
        kraus.append(s.ImmutableSparseMatrix(s.diag(*signs)*evolution.extract(rows, range(16))))
    gram = sum((a.H*a for a in kraus), s.zeros(16)).applyfunc(s.expand)
    if gram != s.eye(16):
        raise ArithmeticError('a finite ancillary output was discarded')
    return y, h, tuple(kraus)


def finite_control():
    y, h, kraus = _finite_control()
    return {'family_coefficient': y, 'physical_pairing': h, 'kraus': kraus}


def certificate():
    controls = []
    for phase in (s.S.One, (s.Integer(99)+20*s.I)/101):
        _, y = _source_lift(phase)
        controls.append({'phase_exact': str(phase),
            'Y_exact': [[str(x) for x in row] for row in y.tolist()],
            'trace_Y_exact': str(s.trace(y)),
            'physical_charged_pairing_norm_exact_at_unit_Higgs': str(2*s.trace(y)),
            'source_coefficient_equal': True})
    y, h, kraus = _finite_control()
    probabilities = [s.expand(sum(s.conjugate(x)*x for x in k)/16)
                     for k in kraus]
    return {'source': 'coherent Koide analyzer, full CAR ancillary-mode lift',
            'physical_target_modes': 90, 'enlarged_modes': 720,
            'fresh_ancillary_modes_per_round': 630,
            'flag_coordinate_dimension': 8, 'family_coordinate_dimension': 3,
            'flag_is_gauge_neutral_but_ancillary_fermions_keep_source_charges': True,
            'full_Fock_vacuum_compression_scope': 'all physical occupations and untouched references',
            'current_and_Higgs_covariance': higgs_covariance(),
            'source_controls': controls,
            'archive_example': archive_bound(65536, 1),
            'archive_has_classical_vacuum_flags_and_full_quantum_complements': True,
            'common_H0_requires_reducing_vacuum_subspace_and_actual_implementation': True,
            'finite_control': {'physical_modes': 4, 'ancillary_modes': 4,
                'full_Fock_dimension': 256, 'physical_input_columns': 16,
                'complete_ancillary_occupation_labels': 16,
                'original_tracial_probabilities_exact': list(map(str, probabilities)),
                'complex_mean_offdiagonal': str(y[0, 1]),
                'exact_full_clock_unitarity': True,
                'all_occupation_compression': True,
                'fine_occupation_readout_claimed_nonabelian_gauge_invariant': False},
            'nonvacuum_quantum_archives_discarded': False,
            'coherent_ancillary_archive_identified_with_classical_flags': False,
            'native_full_720_mode_tape_executed': False,
            'source_gate_overhead_equals_target_pulse_clock': False,
            'source_global_physical_clock_selected': False,
            'interacting_ground_or_particle_pole_selected': False,
            'phase_or_physical_pole_masses_derived': False,
            'same_as_four_current_gated_actions': False,
            'physical_particle_assignments': []}
