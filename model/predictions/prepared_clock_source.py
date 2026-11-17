"""Source words for the actual finite, successfully prepared edge reference.

The declared sampler and herald select M=sum_n w_n (I-L/120)^(2n).
Its positive real logarithm fixes H=-log(M/m(0)), without spectral fitting or
integer phase branches.  The source clock is exp(-i H), so its Heisenberg
action obeys sigma_s=alpha_{-s} for rho=M/Tr(M).  The separate sqrt(L) clock
is a comparison, not the modular generator of this finite prepared state.

All routines act on the same sixteen-address protected register.  Four
padding addresses are retained.  Matrices returned as ``raw`` are executed
source words or products of those executed routines, never inserted targets
or polar-corrected matrices.  Fixed numerical guards are explicitly not
directed-rounding certificates or physical apparatus error measurements.
"""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from . import inverse_clock
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*inverse_clock.SOURCE_PATHS,
    'model/predictions/prepared_clock_source.py',
    'model/observer/source_scalar/compiler.py'))))
DIMENSION = 16
QUBITS = 4
SYNTHESIS_BITS = 44
BASIS_GUARD = F(1, 10**8)
PHASE_GUARD = F(1, 10**10)
PAIR_GUARD = F(1, 10**9)
BAND_ADDRESSES = ((0,), (1, 3, 5), (2, 6, 7, 8, 9), (4, 10, 11))
PAIR_ADDRESSES = (1, 2, 4)


def _dyadic(expression):
    """Evaluate symbolic inputs before one reproducible source-grid rounding."""
    return F(int(sp.floor(sp.N(expression, 90)*(1 << SYNTHESIS_BITS)
                          + sp.Rational(1, 2))), 1 << SYNTHESIS_BITS)


def _quantize(word):
    return [compiler.pulse(event['kind'],
            round(event['angle']*(1 << SYNTHESIS_BITS))/(1 << SYNTHESIS_BITS),
            event['targets']) for event in word]


def _float_matrix(matrix):
    return np.array([[float(sp.N(x, 80)) for x in row]
                     for row in matrix.tolist()])


def _basis_word(exact_basis):
    """Gray elimination at fixed precision before rounding emitted angles.

    The public source compiler's primitive Gray rotations and Walsh phases
    are reused.  Deriving the elimination angles here from the exact basis
    avoids making source-word identities depend on platform libm residuals
    at the many zero pivots of a degenerate eigenspace.
    """
    gray = [j ^ (j >> 1) for j in range(16)]
    rotations = []
    with mp.workdps(90):
        residual = [[mp.mpf(str(sp.N(exact_basis[a, b], 100)))
                     for b in gray] for a in gray]
        for column in range(15):
            for row in range(15, column, -1):
                x, y = residual[row-1][column], residual[row][column]
                theta = mp.mpf(0) if x == 0 and y == 0 else mp.atan2(y, x)
                c, s = mp.cos(theta), mp.sin(theta)
                first, second = residual[row-1][:], residual[row][:]
                residual[row-1] = [c*a+s*b for a, b in zip(first, second, strict=True)]
                residual[row] = [-s*a+c*b for a, b in zip(first, second, strict=True)]
                angle = F(int(mp.floor(theta*(1 << SYNTHESIS_BITS)+mp.mpf('0.5'))),
                          1 << SYNTHESIS_BITS)
                rotations.append((gray[row-1], gray[row], angle))
        signs = [1 if residual[j][j] >= 0 else -1 for j in range(16)]
        if max(abs(residual[a][b]-(signs[a] if a == b else 0))
               for a in range(16) for b in range(16)) > mp.mpf('1e-75'):
            raise ArithmeticError('fixed-precision Gray elimination lost its complete basis')
    phases = [0.]*16
    for j, sign in enumerate(signs):
        phases[gray[j]] = 0. if sign == 1 else float(_dyadic(sp.pi))
    word = compiler.compile_diagonal(phases)
    for a, b, angle in reversed(rotations):
        word.extend(compiler._gray_rotation(a, b, -float(angle), 4))
    return _quantize(word)


@lru_cache(maxsize=1)
def _spectral_data():
    """Exact, deterministic Gram-Schmidt inside each native eigenspace."""
    laplacian = sp.zeros(12)
    seams = tuple(inverse_clock.continuation.quantum_mixing.seams())
    if len(seams) != 30 or len(set(seams)) != 30:
        raise ArithmeticError('the declared native seam census changed')
    for a, b in seams:
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    bands = (sp.Integer(0), 5-sp.sqrt(5), sp.Integer(6), 5+sp.sqrt(5))
    basis = sp.zeros(16)
    values = [sp.Integer(0)]*16
    for eigenvalue, addresses in zip(bands, BAND_ADDRESSES, strict=True):
        vectors = (laplacian-eigenvalue*sp.eye(12)).nullspace()
        if len(vectors) != len(addresses):
            raise ArithmeticError('the native eigenspace multiplicity changed')
        orthogonal = []
        for address, vector in zip(addresses, vectors, strict=True):
            for previous in orthogonal:
                vector = (vector-previous.dot(vector)/previous.dot(previous)*previous
                          ).applyfunc(sp.simplify)
            first = next(x for x in vector if x != 0)
            if float(sp.N(first, 80)) < 0:
                vector = -vector
            if any(sp.simplify(x) != 0 for x in laplacian*vector-eigenvalue*vector):
                raise ArithmeticError('an exact basis vector left its native eigenspace')
            squared = sp.simplify(vector.dot(vector))
            if squared.is_positive is not True:
                raise ArithmeticError('an exact native mode has no positive norm')
            for previous in orthogonal:
                if sp.simplify(previous.dot(vector)) != 0:
                    raise ArithmeticError('exact native mode orthogonality failed')
            orthogonal.append(vector)
            basis[:12, address] = vector/sp.sqrt(squared)
            values[address] = eigenvalue
    for address in range(12, 16):
        basis[address, address] = 1
    array = _float_matrix(basis)
    if (np.linalg.norm(array.T@array-np.eye(16), 'fro') > 2e-13
            or any(sp.simplify(x-1/sp.sqrt(12)) != 0 for x in basis[:12, 0])):
        raise ArithmeticError('the normalized native basis or constant mode changed')
    return {'laplacian': sp.ImmutableMatrix(laplacian),
        'basis_exact': sp.ImmutableMatrix(basis), 'basis': array,
        'eigenvalues': tuple(values), 'band_eigenvalues': bands,
        'band_addresses': BAND_ADDRESSES, 'pair_addresses': PAIR_ADDRESSES,
        'seams': seams,
        'exact_basis_certificate': {
            'symmetric_source_laplacian': laplacian == laplacian.T,
            'all_twelve_eigenvector_equations_checked_exactly': True,
            'within_each_band_orthogonality_and_positive_norms_checked_exactly': True,
            'cross_band_orthogonality': 'distinct eigenvalues of the exact real symmetric L',
            'full_space_completeness': '1+3+5+3=12 mutually orthogonal nonzero vectors',
            'padding_extension': 'four orthogonal standard coordinate vectors'}}


def native_spectral_data():
    """Detached exact source L and basis; no caller supplies an eigenvalue."""
    return deepcopy(_spectral_data())


@lru_cache(maxsize=1)
def _reference_data():
    spectrum = _spectral_data()
    weights = inverse_clock.subordinated_filter_weights()
    if (any(type(w) is not F or w < 0 for w in weights)
            or not 0 < weights[0] <= sum(weights) <= 1
            or not any(weights[1:])):
        raise ArithmeticError('the declared finite sampler is not positive and nontrivial')
    lazy = sp.eye(12)-spectrum['laplacian']/120
    polynomial = sp.zeros(12)
    for n, weight in enumerate(weights):
        polynomial += sp.Rational(weight)*lazy**(2*n)
    mass0 = sum(weights)
    masses = tuple(sp.simplify(sum((sp.Rational(w)*(1-lam/120)**(2*n)
                                  for n, w in enumerate(weights)), sp.Integer(0)))
                   for lam in spectrum['eigenvalues'][:12])
    gaps = tuple(sp.log(sp.Rational(mass0)/mass) for mass in masses)
    normalization = sp.trace(polynomial)
    density = polynomial/normalization
    if (normalization <= 0 or any(not 0 < float(sp.N(mass, 80)) <= float(mass0)
                                  for mass in masses)
            or sp.simplify(sum(masses)-normalization) != 0):
        raise ArithmeticError('finite prepared reference lost its positive spectral law')
    extended_density = np.zeros((16, 16), complex)
    extended_density[:12, :12] = _float_matrix(density)
    return {'weights': weights, 'lazy_filter': sp.ImmutableMatrix(lazy),
        'polynomial': sp.ImmutableMatrix(polynomial),
        'normalization': normalization, 'm0': mass0,
        'density_exact': sp.ImmutableMatrix(density), 'density': extended_density,
        'masses': masses, 'modular_gaps_exact': gaps,
        'modular_gaps': tuple(float(sp.N(gap, 80)) for gap in gaps),
        'reference_algebra': 'B(C^12); the four retained padding addresses are outside its support',
        'modular_convention': 'U(t)=exp(-i*t*H); sigma_s=alpha_(-s), H=-log(M/m0)',
        'logarithm_choice': 'positive real logarithm of the positive prepared polynomial',
        'temperature_rescaled_to_match_sqrt_L': False,
        'padding_clock_extension': 'identity',
        'selection': 'fixed supplied subordination sampler and all-success herald'}


def prepared_reference():
    """Actual finite sampler reference, not its asymptotic sqrt(L) target."""
    return deepcopy(_reference_data())


def _band(band):
    if type(band) is not int or band not in (1, 2, 3):
        raise ValueError('an original nonconstant native band 1, 2 or 3 is required')
    return PAIR_ADDRESSES[band-1]


def _pair_spectral_word(band, axis, angle):
    """Rank-two X/Y pulse; the other fourteen addresses are left untouched."""
    address = _band(band)
    if type(axis) is not str or axis not in ('X', 'Y'):
        raise ValueError('a declared X or Y pair axis is required')
    target = address.bit_length()-1
    result = []
    for mask in range(16):
        if mask & address:
            continue
        label = [axis if j == target else 'Z' if mask >> j & 1 else 'I'
                 for j in range(4)]
        result.extend(compiler.pauli_rotation(label, float(angle)/8))
    return _quantize(result)


def _pair_target(band, axis, angle):
    address = _band(band)
    result = np.eye(16, dtype=complex)
    c, s = math.cos(float(angle)), math.sin(float(angle))
    result[0, 0] = result[address, address] = c
    result[0, address] = 1j*s if axis == 'X' else s
    result[address, 0] = 1j*s if axis == 'X' else -s
    return result


def _check_matrix(raw, target, guard, name, dimension=16):
    if type(dimension) is not int or dimension not in (16, 32):
        raise ValueError('the complete four- or five-qubit register is required')
    for label, matrix in (('raw source', raw), ('target', target)):
        if (type(matrix) is not np.ndarray or matrix.shape != (dimension, dimension)
                or matrix.dtype.kind not in 'iufc' or not np.all(np.isfinite(matrix))):
            raise ArithmeticError(f'{label} {name} is not a complete finite register matrix')
    if (type(guard) is not F or type(guard.numerator) is not int
            or type(guard.denominator) is not int or guard.numerator < 0
            or guard.denominator <= 0 or math.gcd(guard.numerator, guard.denominator) != 1):
        raise ArithmeticError(f'source {name} requires a finite nonnegative rational guard')
    try:
        error_limit = float(guard)
        defect_limit = float(2*guard+guard*guard)
    except OverflowError as exc:
        raise ArithmeticError(
            f'source {name} requires a finite nonnegative rational guard') from exc
    if not math.isfinite(error_limit) or not math.isfinite(defect_limit):
        raise ArithmeticError(f'source {name} requires a finite nonnegative rational guard')
    with np.errstate(over='ignore', invalid='ignore'):
        # Integer matrix products can wrap to an apparently unitary matrix;
        # lower-precision floating norms can also underflow before comparison.
        numeric_raw = raw.astype(np.result_type(raw.dtype, np.complex128), copy=False)
        numeric_target = target.astype(np.result_type(target.dtype, np.complex128), copy=False)
        difference = float(np.linalg.norm(numeric_raw-numeric_target, 'fro'))
        defect = float(np.linalg.norm(numeric_raw.conj().T@numeric_raw-np.eye(dimension), 'fro'))
    if (not math.isfinite(difference) or not math.isfinite(defect)
            or difference > error_limit or defect > defect_limit):
        raise ArithmeticError(f'raw source {name} exceeds its complete-register numerical guard')
    return {'complete_target_Frobenius_error': difference,
        'raw_unitarity_Frobenius_defect': defect,
        'complete_target_error_guard': guard,
        'interval_or_physical_error_certificate': False}


def _execute(word, target, guard, name):
    raw = compiler.execute(word, 4)
    return raw, _check_matrix(raw, target, guard, name)


def _combined_guard(*guards):
    result = F(1)
    for guard in guards:
        result *= 1+guard
    return result-1


@lru_cache(maxsize=1)
def _compiled_data():
    spectrum, reference = _spectral_data(), _reference_data()
    basis = spectrum['basis']
    words, raw, targets, guards, diagnostics = {}, {}, {}, {}, {}
    words['basis'] = _basis_word(spectrum['basis_exact'])
    words['unbasis'] = compiler.inverse(words['basis'])
    targets.update(basis=basis, unbasis=basis.T)
    for name in ('basis', 'unbasis'):
        guards[name] = BASIS_GUARD
        raw[name], diagnostics[name] = _execute(words[name], targets[name], guards[name], name)
    for name, values in (
            ('modular_phase', reference['modular_gaps_exact']),
            ('sqrt_phase', tuple(sp.sqrt(x) for x in spectrum['eigenvalues'][:12]))):
        # Padding has exactly zero requested phase.  Targets independently use
        # high-precision unrounded symbolic phases and all sixteen columns.
        phases = tuple(-_dyadic(value) for value in values)+(F(0),)*4
        words[name] = _quantize(compiler.compile_diagonal([float(x) for x in phases]))
        wanted = np.array([-float(sp.N(value, 80)) for value in values]+[0.]*4)
        targets[name] = np.diag(np.exp(1j*wanted))
        guards[name] = PHASE_GUARD
        raw[name], diagnostics[name] = _execute(words[name], targets[name], guards[name], name)
    for name, phase in (('clock1', 'modular_phase'), ('sqrt_clock1', 'sqrt_phase')):
        words[name] = words['unbasis']+words[phase]+words['basis']
        raw[name] = raw['basis']@raw[phase]@raw['unbasis']
        targets[name] = basis@targets[phase]@basis.T
        guards[name] = _combined_guard(BASIS_GUARD, PHASE_GUARD, BASIS_GUARD)
        diagnostics[name] = _check_matrix(raw[name], targets[name], guards[name], name)
    for band in (1, 2, 3):
        settings = (
            (f'kick_plus_{band}', 'X', sp.pi/4, True),
            (f'kick_minus_{band}', 'X', -sp.pi/4, True),
            (f'detector_x_{band}', 'Y', sp.pi/4, False),
            (f'detector_y_{band}', 'X', -sp.pi/4, False))
        for name, axis, angle, return_to_site in settings:
            pair_word = _pair_spectral_word(band, axis, _dyadic(angle))
            target = _pair_target(band, axis, float(sp.N(angle, 80)))
            pair_raw, pair_diagnostic = _execute(pair_word, target, PAIR_GUARD, name+'_pair')
            words[name] = words['unbasis']+pair_word
            raw[name] = pair_raw@raw['unbasis']
            targets[name] = target@basis.T
            guards[name] = _combined_guard(BASIS_GUARD, PAIR_GUARD)
            if return_to_site:
                words[name] += words['basis']
                raw[name] = raw['basis']@raw[name]
                targets[name] = basis@targets[name]
                guards[name] = _combined_guard(BASIS_GUARD, PAIR_GUARD, BASIS_GUARD)
            diagnostics[name] = _check_matrix(raw[name], targets[name], guards[name], name)
            diagnostics[name]['pair_routine'] = pair_diagnostic
    digests = {name: hashlib.sha256(json.dumps(word, sort_keys=True,
                separators=(',', ':'), allow_nan=False).encode()).hexdigest()
               for name, word in words.items()}
    # Exact-algebra identity, numerically checked on the actual finite density;
    # the source implementation allowance is retained separately above.
    state = reference['density']
    stationarity = float(np.linalg.norm(targets['clock1']@state
                                       @targets['clock1'].conj().T-state, 'fro'))
    if stationarity > 2e-13:
        raise ArithmeticError('finite reference and its positive-log clock disagree')
    return {'words': words, 'raw': raw, 'targets': targets, 'guards': guards,
        'diagnostics': diagnostics, 'word_sha256': digests,
        'pulse_counts': {name: len(word) for name, word in words.items()},
        'pair_addresses': {band: _band(band) for band in (1, 2, 3)},
        'reference': reference, 'spectrum': spectrum,
        'ideal_reference_stationarity_error': stationarity,
        'source_synthesis_grid': F(1, 1 << SYNTHESIS_BITS),
        'raw_source_routines_used_without_polar_correction': True,
        'complete_sixteen_address_maps_checked_with_numerical_guards': True,
        'modular_clock_is_sqrt_L_clock': False,
        'physical_seconds_or_autonomous_clock_claimed': False}


def compiled_source():
    """Fresh words and raw source maps for the fixed declared experiment.

    ``clock1`` and ``sqrt_clock1`` act in the site basis.  ``kick_plus_b``
    and ``kick_minus_b`` are exp(+/-i*pi*X_0j/4) in the same register.
    ``detector_x_b`` and ``detector_y_b`` end in the spectral address basis;
    address-zero minus address-j has expectation X_0j or Y_0j.  All other
    addresses, including padding, remain outcomes.  Band 2 is lambda=6.

    Reusing the executed raw routines is linear source execution, not a
    replacement of the prepared density.  Untouched reference columns can
    be carried through each returned map without modification.
    """
    return deepcopy(_compiled_data())


def clock_word(kind='modular'):
    if type(kind) is not str or kind not in ('modular', 'sqrt'):
        raise ValueError('the declared modular or sqrt comparison clock is required')
    return deepcopy(_compiled_data()['words']['clock1' if kind == 'modular' else 'sqrt_clock1'])


def pair_word(band, sign=1):
    _band(band)
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('an original signed pair kick is required')
    return deepcopy(_compiled_data()['words'][f'kick_{"plus" if sign == 1 else "minus"}_{band}'])


def detector_word(band, axis='X'):
    _band(band)
    if type(axis) is not str or axis not in ('X', 'Y'):
        raise ValueError('a declared X or Y detector is required')
    return deepcopy(_compiled_data()['words'][f'detector_{axis.lower()}_{band}'])


@lru_cache(maxsize=8)
def _pair_pointer_data(factors, axis):
    compiled = _compiled_data()
    controlled, wanted = [], np.eye(16, dtype=complex)
    for band, direction, angle in factors:
        address = _band(band)
        target = address.bit_length()-1
        for mask in range(16):
            if mask & address:
                continue
            label = [direction if j == target else 'Z' if mask >> j & 1 else 'I'
                     for j in range(4)]
            # Eight commuting Pauli terms give the pair operator; the extra
            # (I-Z_pointer)/2 controls it on the pointer-one branch.
            controlled.extend(compiler.pauli_rotation(label+['I'], float(angle)/16))
            controlled.extend(compiler.pauli_rotation(label+['Z'], -float(angle)/16))
        wanted = _pair_target(band, direction, angle)@wanted
    controlled = _quantize(controlled)
    target_controlled = np.eye(32, dtype=complex)
    target_controlled[16:, 16:] = wanted
    raw_controlled = compiler.execute(controlled, 5)
    control_guard = len(factors)*PAIR_GUARD
    control_check = _check_matrix(raw_controlled, target_controlled, control_guard,
                                 'controlled_pair_factors', 32)
    prepare = _quantize([compiler.pulse('Y', -float(_dyadic(sp.pi/4)), [4])])
    read = _quantize([compiler.pulse('Y' if axis == 'X' else 'X',
                                    (1 if axis == 'X' else -1)*float(_dyadic(sp.pi/4)), [4])])
    root = math.sqrt(.5)
    target_prepare = np.kron(np.array([[root, -root], [root, root]]), np.eye(16))
    target_read = np.kron(np.array([[root, root], [-root, root]]) if axis == 'X'
                         else np.array([[root, -1j*root], [-1j*root, root]]), np.eye(16))
    raw_prepare, raw_read = compiler.execute(prepare, 5), compiler.execute(read, 5)
    _check_matrix(raw_prepare, target_prepare, PHASE_GUARD, 'pointer_prepare', 32)
    _check_matrix(raw_read, target_read, PHASE_GUARD, 'pointer_read', 32)
    raw_basis = np.kron(np.eye(2), compiled['raw']['basis'])
    raw_unbasis = np.kron(np.eye(2), compiled['raw']['unbasis'])
    target_basis = np.kron(np.eye(2), compiled['spectrum']['basis'])
    word = prepare+compiled['words']['unbasis']+controlled+compiled['words']['basis']+read
    raw = raw_read@raw_basis@raw_controlled@raw_unbasis@raw_prepare
    target_map = target_read@target_basis@target_controlled@target_basis.T@target_prepare
    # Tensoring a 16-address Frobenius error with I_2 costs sqrt(2)<2.
    guard = _combined_guard(PHASE_GUARD, 2*BASIS_GUARD, control_guard,
                            2*BASIS_GUARD, PHASE_GUARD)
    diagnostic = _check_matrix(raw, target_map, guard, 'ordered_pair_pointer', 32)
    return {'word': word, 'raw': raw, 'target': target_map, 'guard': guard,
        'word_sha256': hashlib.sha256(json.dumps(word, sort_keys=True,
            separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
        'diagnostic': diagnostic, 'controlled_diagnostic': control_check,
        'pulse_count': len(word), 'pointer_qubit': 4, 'complete_outcome_count': 32,
        'proper_code_load_events': len(word), 'proper_code_export_events': len(word),
        'protected_work_and_pointer_qubits': 5,
        'chronological_factors': factors, 'pointer_axis': axis,
        'spectral_ordered_unitary': wanted,
        'expectation': 'Re Tr(rho U)' if axis == 'X' else 'Im Tr(rho U)',
        'multiplication_order': 'U=last_factor ... first_factor',
        'input_pointer': 'zero', 'reference_columns_untouched': True,
        'growing_copy_register_execution_claimed': False}


def pair_pointer_source(factors, axis='X'):
    """Actual controlled pair words and pointer readout on all 32 addresses.

    ``factors`` is a chronological tuple of (band, X/Y, rational angle).
    Apply ``raw`` to pointer-zero times the unchanged prepared work/reference
    columns. Sum every output probability with sign (-1)^pointer for the
    indicated characteristic-function component; work addresses are retained.
    """
    if type(axis) is not str or axis not in ('X', 'Y'):
        raise ValueError('a declared X or Y pointer readout is required')
    if type(factors) is not tuple or not 1 <= len(factors) <= 16:
        raise ValueError('one through sixteen chronological pair factors are required')
    for factor in factors:
        if type(factor) is not tuple or len(factor) != 3:
            raise ValueError('each factor requires a band, axis and exact rational angle')
        band, direction, angle = factor
        _band(band)
        if type(direction) is not str or direction not in ('X', 'Y'):
            raise ValueError('a declared X or Y pair axis is required')
        if type(angle) is not F or not -64 <= angle <= 64:
            raise ValueError('an exact rational pair angle in [-64,64] is required')
    return deepcopy(_pair_pointer_data(factors, axis))


def _negative_exponential_lower(value):
    """Rational lower bound, using the sampler's odd alternating polynomial."""
    if not value:
        return F(1)
    pieces = math.ceil(value)
    x = value/pieces
    lower = sum(((-x)**j/math.factorial(j) for j in range(14)), F())**pieces
    grid = 1 << 80
    rounded = F((lower*grid).numerator//(lower*grid).denominator, grid)
    # e<3 also gives a strictly positive fallback when the reporting grid
    # cannot resolve this band. It does not require floating exponentiation.
    return max(rounded, F(1, 3**pieces))


def positive_log_clock_bound(rate, heat_step, band, missing_mass, maximum_requests,
                             filter_operator_error=F(0), *, source_preserves_band=False):
    """Bound the prepared polynomial's log, not a separately chosen clock.

    Premises: D>=0; F=I-hD/2 is a contraction on the encoded domain; retained
    nonnegative weights are coefficientwise below the Sibuya/Poisson law,
    with total mass >=1-delta and cap N. On 0<=D<=Lambda, require hLambda<2.
    The infinite law is exp[-b sqrt(D-hD²/4)], b=R sqrt(h). Normalize the
    retained M by its own zero-reference mass m0=sum w, not by one.

    If a<=exp(-b sqrt(Lambda)), positivity follows from delta<a. Writing
    M=exp(-g)-E, 0<=E<=delta, shows
      |log(m0)-log(M)-g| <= -log(1-delta/a) <= delta/(a-delta).
    Both log corrections are nonpositive, so their difference needs no
    extra factor two. A raw successful filter within epsilon of F gives
    ||M_raw-M|| <= (1+epsilon)^(2N)-1 <= 2N epsilon/(1-2N epsilon).
    The integral formula for log bounds noncommuting source errors only if
    the encoded domain lies wholly in this band or the raw implementation
    preserves its spectral subspace. A small off-band operator error alone
    gives no such bound when outside eigenvalues can approach zero.
    These are conditional analytic bounds, not a new growing compiler.
    """
    rate, heat_step, band = (inverse_clock._fraction(x, 0)
                            for x in (rate, heat_step, band))
    if not rate or not heat_step or not band:
        raise ValueError('positive rate, heat step and bounded spectral band are required')
    missing = inverse_clock._fraction(missing_mass, 0, 1)
    maximum = inverse_clock._integer(maximum_requests, 0, 10**100)
    epsilon = inverse_clock._fraction(filter_operator_error, 0)
    if type(source_preserves_band) is not bool:
        raise ValueError('an explicit Boolean spectral-band preservation premise is required')
    if heat_step*band >= 2:
        raise ValueError('the bounded band requires heat_step*band < 2')
    sqrt_interval = inverse_clock.source.exact._sqrt_interval
    root_band = sqrt_interval(heat_step*band)[1]
    beta = tuple(rate*x for x in sqrt_interval(heat_step))
    mass_floor = _negative_exponential_lower(rate*root_band)
    root_correction = sqrt_interval(1-heat_step*band/4)[0]
    bias = rate*root_band*heat_step*band/(4*(1+root_correction))
    truncation = missing/(mass_floor-missing) if missing < mass_floor else None
    ideal_error = truncation+bias if truncation is not None else None
    amplification = 2*maximum*epsilon
    source_error = amplification/(1-amplification) if amplification < 1 else None
    floor = mass_floor-missing
    logarithm_error = (source_error/(floor-source_error)
                       if (not epsilon or source_preserves_band)
                       and source_error is not None and source_error < floor else None)
    total = (ideal_error+logarithm_error
             if ideal_error is not None and logarithm_error is not None else None)
    return {'rate': rate, 'heat_step': heat_step, 'spectral_band_upper': band,
        'inverse_temperature_interval': beta,
        'infinite_reference_mass_lower': mass_floor,
        'missing_sampler_mass_upper': missing,
        'retained_band_mass_lower': max(F(0), floor),
        'normalization': 'H=log(sum w) I-log(sum_n w_n F^(2n))',
        'logarithm_domain': 'the certified spectral band; global faithfulness is not inferred',
        'normalized_log_truncation_upper': truncation,
        'filter_to_square_root_generator_upper': bias,
        'ideal_generator_to_b_sqrt_D_upper': ideal_error,
        'maximum_filter_requests': maximum,
        'one_filter_operator_error': epsilon,
        'raw_positive_mixture_operator_error_upper': source_error,
        'raw_positive_mixture_log_error_upper': logarithm_error,
        'raw_generator_to_b_sqrt_D_upper': total,
        'unit_model_time_unitary_error_upper': total/beta[0] if total is not None else None,
        'unit_model_time_observable_error_upper': 2*total/beta[0] if total is not None else None,
        'source_log_scalar_convention': 'same ideal m0; scalar changes cancel in modular action',
        'raw_implementation_preserves_band_or_entire_encoded_domain_is_in_band': source_preserves_band,
        'off_band_source_log_bound_claimed': False,
        'fixed_cap_generator_upper': maximum*heat_step*band/(1-heat_step*band/2),
        'maximum_encoded_generator_norm_for_filter_contraction': 4/heat_step,
        'full_encoded_domain_filter_contraction_is_a_premise': True,
        'retained_request_algebra_has_this_modular_generator': False,
        'growing_source_filter_or_phase_compiler_certified': False}


def prepared_log_limit_bound(rate=F(1, 2), weight=F(1, 2), jump_cap=2,
                             count_cap=2, order=13, sampling_bits=64,
                             filter_operator_error=F(0)):
    """Apply the log modulus to the existing executable positive sampler.

    The success herald selects precisely rows0:12 of the 32768-by-12 raw
    filter; all other rows, including the four work-padding addresses, are
    retained failures. Each successful S is therefore a map C12->C12 by
    projection, whose whole declared L spectrum is below8. This is not a
    claim that the later16-address clock exactly preserves that subspace.
    """
    inherited = inverse_clock.subordinated_filter_bound(
        rate, weight, jump_cap, count_cap, order, sampling_bits)
    result = positive_log_clock_bound(inherited['rate'], inherited['heat_step'], F(8),
        inherited['sampling_rejection_probability'], inherited['maximum_filter_requests'],
        filter_operator_error, source_preserves_band=True)
    result['inherited_constant_to_band_modular_gap_upper'] = (
        inherited['constant_to_band_modular_gap_error_upper'])
    result['sampler_coefficients_checked_by_existing_positive_evaluator'] = True
    result['successful_filter_domain_dimension'] = 12
    result['successful_filter_domain_selected_by_herald_projection'] = True
    return result


def positive_log_refinement_schedule(level, inverse_temperature=F(1), band=F(8), *,
                                    source_preserves_band=False):
    """Explicit finite cap/precision obligations; not an executed large sampler.

    Put R=2^l, h=b²/R², B=2^(6l), J=6R+2l, odd Taylor order 4l+3.
    The existing tail theorem, central-binomial bound <=1/sqrt(B+1),
    k!>=(k/e)^k with e<3, and factorial >=2^order give missing mass
    <=4*2^(-2l), including a sufficiently fine dyadic sampler. Set one-filter
    error <=2^(-2l)/(8BJ). On every fixed band the log error tends to zero.
    Raw log convergence additionally requires a band-preserving compiler or
    a code whose whole spectrum is in the band. Fixed caps or fixed source
    precision do not implement this limit.
    """
    level = inverse_clock._integer(level, 0, 12)
    beta, band = (inverse_clock._fraction(x, F(1, 10**12))
                  for x in (inverse_temperature, band))
    rate, scale = 1 << level, F(1, 1 << (2*level))
    heat_step = beta*beta/(rate*rate)
    jumps, counts = 1 << (6*level), 6*rate+2*level
    maximum = jumps*counts
    bits = max(48, maximum.bit_length()+2*level)
    result = positive_log_clock_bound(F(rate), heat_step, band, min(F(1), 4*scale),
        maximum, scale/(8*maximum), source_preserves_band=source_preserves_band)
    result.update(level=level, jump_cap=jumps, count_cap=counts,
        lower_Taylor_order=4*level+3, sampling_bits=bits,
        jump_tail_upper=scale, Poisson_count_tail_upper=scale,
        Taylor_remainder_upper=scale, dyadic_rounding_upper=F(maximum+1, 1 << bits),
        required_phase_word_error_upper=scale,
        growing_sampler_evaluated=False,
        existing_finite_evaluator_covers_parameters=(rate <= 32 and jumps <= 128
            and counts <= 16 and 4*level+3 <= 61 and bits <= 128),
        bounded_band_then_spectral_cutoff_limit=True,
        unbounded_generator_operator_norm_convergence_claimed=False)
    return result


def native_growing_word_bounds(vertices, seams):
    """Constructive pulse counts for the existing source-gate identities.

    Binary-tree preparation has B(n)=3*((n-1)*2**n+1) pulses. A complete
    Walsh diagonal has C(n)=(3*n-5)*2**n+6. One Householder reflection
    costs 2B+C, so a real V-dimensional native spectral basis, extended
    identically on padding, needs at most V such reflections. This is a
    finite-word existence bound, not a call beyond the executable API caps.

    The event shift uses the existing clean AND ladders and Gray paths.
    An r-controlled flip uses 2r-3 forty-pulse Toffolis; a seam swaps at
    most 2n-1 adjacent work addresses. Every zero-angle request is counted.
    """
    vertices = inverse_clock._integer(vertices, 3)
    seams = inverse_clock._integer(seams, 2)
    n, r = (vertices-1).bit_length(), (2*seams-1).bit_length()
    tree = lambda bits: 3*((bits-1)*(1 << bits)+1)
    diagonal = lambda bits: (3*bits-5)*(1 << bits)+6
    toffoli = seams*(2*(2*r-3)+(2*n-1)*(2*n-3))
    flips = seams*(2*r+2*(n-1)*(2*n-1))
    shift = 40*toffoli+2*flips
    preparation = tree(n)+3*n  # Uniform valid-address Bell purification.
    filter_pulses = 2*tree(r)+shift
    encoder = vertices*(2*tree(n)+diagonal(n))
    return {'native_addresses': vertices, 'native_seams': seams,
        'work_qubits': n, 'event_qubits': r,
        'clean_flag_and_scratch_qubits': max(n, r)-1,
        'peak_qubits_with_retained_purifying_reference': 2*n+r+max(n, r)-1,
        'retained_work_padding_addresses': (1 << n)-vertices,
        'retained_event_padding_addresses': (1 << r)-2*seams,
        'binary_tree_work_preparation_pulses': tree(n),
        'binary_tree_event_preparation_pulses': tree(r),
        'complete_work_diagonal_pulses': diagonal(n),
        'uniform_input_purification_pulses': preparation,
        'shift_Toffoli_upper': toffoli, 'shift_flip_upper': flips,
        'native_event_shift_pulses_upper': shift,
        'one_native_lazy_filter_pulses_upper': filter_pulses,
        'full_real_spectral_encoder_pulses_upper': encoder,
        'one_same_state_spectral_clock_word_pulses_upper': 2*encoder+diagonal(n),
        'controlled_complex_star_factor_pulses_upper_in_spectral_basis':
            4*tree(n)+4*diagonal(n)+3*(n+1)*(1 << n),
        'ordered_pointer_common_encoder_and_readout_pulses_upper': 2*encoder+2,
        'one_proper_code_load_and_export_per_pulse': True,
        'all_failure_padding_and_reference_outputs_required': True,
        'native_matrix_spectral_synthesis_is_classical_computation': True,
        'classical_synthesis_runtime_or_physical_duration_included': False,
        'unbounded_word_or_register_executed': False}


def native_growing_source_resources(level, inverse_temperature=F(1), band=F(8), *,
                                    ordered_factors=2):
    """One jointly realizable native filter, encoder and precision schedule.

    D_l=4**l L_l, V_l=240*4**l, M_l=630*4**l. The biased native-event
    filter gives I-hD_l/2 iff p=630*h*16**l<=1. Choose r=2l+c with
    4**c>=b²*max(630,Lambda); R=2**r and h=b²/R². This also gives
    h||D_l||<=12/(630*4**l)<2, excluding the high-mode F² alias.

    The existing positive sampler has an n=0 branch. Its odd Taylor
    polynomial at one is >=1/3. For integer R and sampling bits>=2R+2,
    w0>=3**(-R)-2**(-2R-2)>=(3/4)*3**(-R)>=2**(-2R).
    No extra identity component or altered reference state is introduced.

    Write f=2**(-2R), a=4**(-r), K=BJ, P=max(Pinit,Pfilter). A combined
    primitive/load/export error <=f*a/(64*V*K*P), in both operator norm
    and complete-instrument half-diamond norm, gives stopped preparation
    half-trace error E<=f*a/(32V). The unnormalized successful marginal
    S_raw, including imperfect initial purification, gives M_raw=V*S_raw.
    Thus ||M_raw-M||<=2VE<=f*a/16 and the global log error is <=a/(16-a).
    Normalizing the actual marginal changes only the scalar logarithm.
    This removes the reducing-band premise, at explicitly costly precision.
    Arbitrarily accurate physical proper-code services are NOT derived.

    Only integer exponents of very small error budgets are returned; the
    function does not allocate 2**(2R), growing words or sampler arrays.
    """
    level = inverse_clock._integer(level, 0, 12)
    ordered_factors = inverse_clock._integer(ordered_factors, 1, 64)
    beta = inverse_clock._fraction(inverse_temperature, F(1, 10**12), F(16))
    band = inverse_clock._fraction(band, F(1, 10**12), F(10**6))
    offset = 0
    while F(1 << (2*offset)) < beta*beta*max(F(630), band):
        offset += 1
    index = 2*level+offset
    rate, scale = 1 << index, F(1, 1 << (2*index))
    heat = beta*beta*scale
    vertices, seams = 240*4**level, 630*4**level
    lazy_weight = 630*heat*16**level
    if not 0 < lazy_weight <= 1 or not 12*heat*4**level < 2:
        raise ArithmeticError('the growing filter is not an admitted native convex event word')
    jumps, counts = 1 << (6*index), 6*rate+2*index
    maximum = jumps*counts
    sampling_bits = max(48, maximum.bit_length()+2*index, 2*rate+2)
    words = native_growing_word_bounds(vertices, seams)
    initial = words['uniform_input_purification_pulses']
    filtering = words['one_native_lazy_filter_pulses_upper']
    pulse_maximum = max(initial, filtering)
    stopped_pulses = initial+maximum*filtering
    # 2**tail_bits >=64 V K P. Only this moderate exponent is materialized.
    tail_bits = (64*vertices*maximum*pulse_maximum-1).bit_length()
    primitive_bits = 2*rate+2*index+tail_bits
    error_over_floor = F(stopped_pulses, 1 << (2*index+tail_bits))
    relative_mixture_error = 2*vertices*error_over_floor
    log_error = relative_mixture_error/(1-relative_mixture_error)
    ideal = positive_log_clock_bound(F(rate), heat, band, min(F(1), 4*scale), maximum)
    ideal_error = ideal['ideal_generator_to_b_sqrt_D_upper']
    if error_over_floor > scale/(32*vertices) or log_error > scale/(16-scale):
        raise ArithmeticError('complete stopped preparation exceeds the global logarithm budget')
    clock_pulses = words['one_same_state_spectral_clock_word_pulses_upper']
    clock_bits = 2*index+(2*clock_pulses-1).bit_length()
    pointer_pulses = (words['ordered_pointer_common_encoder_and_readout_pulses_upper']
        +ordered_factors*(words['controlled_complex_star_factor_pulses_upper_in_spectral_basis']
                         +2*clock_pulses))
    pointer_bits = 2*index+(8*pointer_pulses-1).bit_length()
    success_lower = max(F(0), (1-4*scale)/vertices)
    complete_error = scale/(32*vertices)
    return {'native_level': level, 'sampler_index': index,
        'inverse_temperature': beta, 'spectral_band_upper': band,
        'rate': rate, 'heat_step': heat, 'native_lazy_event_weight': lazy_weight,
        'whole_native_generator_norm_upper': 12*4**level,
        'whole_native_filter_alias_parameter_upper': 12*heat*4**level,
        'jump_cap': jumps, 'count_cap': counts, 'maximum_filter_requests': maximum,
        'lower_Taylor_order': 4*index+3, 'sampling_coin_bits': sampling_bits,
        'missing_sampler_mass_upper': min(F(1), 4*scale),
        'zero_request_mass_lower_base2_negative_exponent': 2*rate,
        'required_primitive_with_services_error_base2_negative_exponent': primitive_bits,
        'required_clock_word_primitive_error_base2_negative_exponent': clock_bits,
        'complete_stopped_preparation_source_pulses_upper': stopped_pulses,
        'complete_preparation_error_divided_by_zero_request_floor_upper': error_over_floor,
        'complete_preparation_half_trace_error_upper': complete_error,
        'successful_mixture_error_divided_by_zero_request_floor_upper': relative_mixture_error,
        'global_actual_marginal_logarithm_error_upper': log_error,
        'ideal_prepared_generator_error_upper': ideal_error,
        'actual_marginal_generator_error_upper':
            ideal_error+log_error if ideal_error is not None else None,
        'clock_word_operator_error_upper': scale,
        'ordered_factors': ordered_factors,
        'one_pointer_quadrature_source_pulses_upper': pointer_pulses,
        'two_quadrature_source_pulses_per_successful_copy_upper': 2*pointer_pulses,
        'required_pointer_primitive_with_services_error_base2_negative_exponent': pointer_bits,
        'one_copy_complex_Weyl_estimator_error_upper': scale,
        'nominal_clock_to_actual_modular_unitary_error_per_model_time_upper': log_error/beta,
        'nominal_clock_to_actual_modular_flow_error_per_model_time_upper': 2*log_error/beta,
        'one_model_time_actual_generator_and_clock_error_upper':
            (ideal_error+log_error)/beta+scale if ideal_error is not None else None,
        'ideal_success_probability_lower': success_lower,
        'actual_success_probability_lower': max(F(0), success_lower-2*complete_error),
        'resources': words,
        'initial_purification_and_all_stopped_failures_included_in_error': True,
        'raw_log_bound_requires_reducing_spectral_band': False,
        'extra_identity_regularizer_added': False,
        'increasing_source_service_precision_is_an_explicit_premise': True,
        'full_refinement_word_sampler_or_physical_accuracy_executed': False,
        'radial_architecture_interactions_or_lab_units_derived': False}
