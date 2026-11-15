"""One quantum electric link and two charged CAR channels on a physical code.

Logical index b0+2*b1 selects right occupation in the two channels, whose
charges are +1 and -1. The isometry attaches integer flux E=-b0+b1 to the
ordered four-mode Fock state. Hopping and electric energy act on this SAME
four-dimensional Gauss sector; they are not unrelated tensor demonstrations.

The infinite integer link, CAR representation, charge pair, initial state,
action and coherent proper-code services are declared construction choices.
No spatial continuum, physical clock, Standard Model identification, or
laboratory prediction is supplied by this finite example.
"""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import math
from numbers import Real
from pathlib import Path

import numpy as np

from model.observer.source_scalar import compiler


Gaussian = tuple[Q, Q]
Matrix = tuple[tuple[Gaussian, ...], ...]
ZERO: Gaussian = (Q(0), Q(0))
ONE: Gaussian = (Q(1), Q(0))
MINUS_I: Gaussian = (Q(0), Q(-1))
CHARGES = (1, -1)
MODES = ("L0", "R0", "L1", "R1")
ANGLE = math.atan2(4, 3)
ROOT = Path(__file__).resolve().parents[2]


def _integer(value, minimum, maximum, name):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


def _plus(a: Gaussian, b: Gaussian) -> Gaussian:
    return a[0]+b[0], a[1]+b[1]


def _times(a: Gaussian, b: Gaussian) -> Gaussian:
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _conjugate(a: Gaussian) -> Gaussian:
    return a[0], -a[1]


def _matrix(rows) -> Matrix:
    return tuple(tuple(row) for row in rows)


def identity(size: int) -> Matrix:
    _integer(size, 1, 256, "matrix dimension")
    return tuple(tuple(ONE if i == j else ZERO for j in range(size))
                 for i in range(size))


def dagger(matrix: Matrix) -> Matrix:
    return tuple(tuple(_conjugate(row[j]) for row in matrix)
                 for j in range(len(matrix[0])))


def matmul(left: Matrix, right: Matrix) -> Matrix:
    if len(left[0]) != len(right):
        raise ValueError("matrix multiplication dimensions disagree")
    rows = [[ZERO for _ in right[0]] for _ in left]
    for i, row in enumerate(left):
        for k, a in enumerate(row):
            if a == ZERO:
                continue
            for j, b in enumerate(right[k]):
                if b != ZERO:
                    rows[i][j] = _plus(rows[i][j], _times(a, b))
    return _matrix(rows)


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    if len(left) != len(right) or len(left[0]) != len(right[0]):
        raise ValueError("matrix addition dimensions disagree")
    return tuple(tuple(_plus(a, b) for a, b in zip(x, y, strict=True))
                 for x, y in zip(left, right, strict=True))


def matrix_scale(matrix: Matrix, coefficient: Gaussian) -> Matrix:
    return tuple(tuple(_times(coefficient, cell) for cell in row) for row in matrix)


def numerical_matrix(matrix: Matrix) -> np.ndarray:
    return np.array([[complex(float(real), float(imag)) for real, imag in row]
                     for row in matrix], dtype=complex)


def encode_matrix(matrix: Matrix):
    return [[[str(real), str(imag)] for real, imag in row] for row in matrix]


def car_ladder(occupation: int, mode: int, creation: bool):
    """Return (fermionic sign, new occupation), or None for a blocked ladder.

    Occupation mask bit j denotes mode j in the ordered basis
    (c0*)^n0 (c1*)^n1 (c2*)^n2 (c3*)^n3 |vacuum>.
    """
    _integer(occupation, 0, 15, "four-mode occupation")
    _integer(mode, 0, 3, "CAR mode")
    if type(creation) is not bool:
        raise ValueError("creation must be Boolean")
    occupied = bool(occupation & (1 << mode))
    if occupied == creation:
        return None
    sign = (-1)**((occupation & ((1 << mode)-1)).bit_count())
    return sign, occupation ^ (1 << mode)


def _ladder_word(occupation, word):
    """Apply (mode, creation) pairs in chronological order."""
    sign = 1
    for mode, creation in word:
        result = car_ladder(occupation, mode, creation)
        if result is None:
            return None
        factor, occupation = result
        sign *= factor
    return sign, occupation


def gauss_values(occupation: int, flux: int):
    """G_L=E-Q_L and G_R=-E-Q_R, on the unbounded integer link."""
    _integer(occupation, 0, 15, "four-mode occupation")
    if type(flux) is not int:
        raise ValueError("electric flux must be an integer")
    left = sum(q*((occupation >> (2*j)) & 1) for j, q in enumerate(CHARGES))
    right = sum(q*((occupation >> (2*j+1)) & 1) for j, q in enumerate(CHARGES))
    return flux-left, -flux-right


def dressed_hopping_images(occupation: int, flux: int, channel: int):
    """Images under cL* U^q cR + cR* U^-q cL, including CAR signs.

    U|ell>=|ell+1> on ell2(Z). No artificial finite-boundary rule appears.
    """
    gauss_values(occupation, flux)  # validates the full ambient basis input
    _integer(channel, 0, 1, "channel")
    left, right = 2*channel, 2*channel+1
    result = []
    for source, target, delta in ((right, left, CHARGES[channel]),
                                  (left, right, -CHARGES[channel])):
        moved = _ladder_word(occupation, ((source, False), (target, True)))
        if moved is not None:
            sign, output = moved
            result.append((sign, output, flux+delta))
    return tuple(result)


def physical_basis():
    """Four isometry columns, indexed by b0+2*b1; flux is not another qubit."""
    result = []
    for index in range(4):
        bits = tuple((index >> j) & 1 for j in range(2))
        occupation = sum(1 << (2*j+bits[j]) for j in range(2))
        flux = -sum(q*b for q, b in zip(CHARGES, bits, strict=True))
        result.append({"index": index, "right_occupation": list(bits),
                       "fermion_occupation_mask": occupation,
                       "occupied_modes": [2*j+bits[j] for j in range(2)],
                       "electric_flux": flux,
                       "Gauss_values": list(gauss_values(occupation, flux))})
    return result


def _channel_hopping(channel: int) -> Matrix:
    _integer(channel, 0, 1, "channel")
    basis = physical_basis()
    lookup = {(row["fermion_occupation_mask"], row["electric_flux"]): row["index"]
              for row in basis}
    matrix = [[ZERO for _ in basis] for _ in basis]
    for column in basis:
        images = dressed_hopping_images(column["fermion_occupation_mask"],
                                        column["electric_flux"], channel)
        for sign, occupation, flux in images:
            target = lookup.get((occupation, flux))
            if target is None:
                raise ArithmeticError("dressed hopping left the physical sector")
            matrix[target][column["index"]] = (Q(sign), Q(0))
    return _matrix(matrix)


def hopping_hamiltonian() -> Matrix:
    return matrix_add(_channel_hopping(0), _channel_hopping(1))


def electric_operator() -> Matrix:
    fluxes = [row["electric_flux"] for row in physical_basis()]
    return tuple(tuple((Q(fluxes[i]), Q(0)) if i == j else ZERO
                       for j in range(4)) for i in range(4))


def electric_energy() -> Matrix:
    return matmul(electric_operator(), electric_operator())


def hamiltonian() -> Matrix:
    """The declared interacting generator H_hop+E² in dimensionless model time."""
    return matrix_add(hopping_hamiltonian(), electric_energy())


def electric_current_operator() -> Matrix:
    """J=-i[H_hop,E], so i[H_hop+E^2,E]=-J in model time."""
    h, e = hopping_hamiltonian(), electric_operator()
    commutator = matrix_add(matmul(h, e), matrix_scale(matmul(e, h), (-Q(1), Q(0))))
    return matrix_scale(commutator, MINUS_I)


def occupation_projector(bit: int, outcome: int) -> Matrix:
    _integer(bit, 0, 1, "right-occupation bit")
    _integer(outcome, 0, 1, "binary outcome")
    return tuple(tuple(ONE if i == j and ((i >> bit) & 1) == outcome else ZERO
                       for j in range(4)) for i in range(4))


def ideal_hopping() -> Matrix:
    """U=exp(-i theta H_hop), cos(theta)=3/5, sin(theta)=4/5."""
    factors = [matrix_add(matrix_scale(identity(4), (Q(3, 5), Q(0))),
                          matrix_scale(_channel_hopping(j), (Q(0), -Q(4, 5))))
               for j in range(2)]
    return matmul(factors[1], factors[0])


def ideal_electric() -> Matrix:
    """V=exp(-i theta E^2); both nonzero supported fluxes have E²=1."""
    result = []
    for i, row in enumerate(physical_basis()):
        value = ONE if row["electric_flux"] == 0 else (Q(3, 5), -Q(4, 5))
        result.append(tuple(value if i == j else ZERO for j in range(4)))
    return tuple(result)


def ideal_echo() -> Matrix:
    u = ideal_hopping()
    return matmul(dagger(u), matmul(ideal_electric(), u))


def branch_probability(unitary: Matrix, input_index: int, bit: int, outcome: int) -> Q:
    """Unconditional probability of one projector branch on a basis input."""
    _integer(input_index, 0, 3, "physical input index")
    projected = matmul(occupation_projector(bit, outcome), unitary)
    return sum((real*real+imag*imag for real, imag in (row[input_index] for row in projected)), Q(0))


def _duration(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError("model duration must be a finite real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("model duration must be a finite real number")
    return result


def source_hopping_word(model_duration=ANGLE):
    duration = _duration(model_duration)
    return (compiler.pauli_rotation("XI", -duration)
            + compiler.pauli_rotation("IX", -duration))


def source_electric_word(model_duration=ANGLE):
    # Includes the global phase. Dropping it changes the declared unitary.
    duration = _duration(model_duration)
    return compiler.compile_diagonal([0., -duration, -duration, 0.])


def source_step_word(model_duration):
    """Symmetric finite step for H_hop+E²; not its exact summed exponential.

    The parameter is dimensionless model duration, with no calibrated clock
    or native processor service duration identification.
    """
    duration = _duration(model_duration)
    half_hopping = source_hopping_word(duration/2)
    return half_hopping+source_electric_word(duration)+half_hopping


def source_echo_word():
    u = source_hopping_word()
    return u+source_electric_word()+compiler.inverse(u)


def _exact_checks():
    basis = physical_basis()
    if any(row["Gauss_values"] != [0, 0] for row in basis):
        raise ArithmeticError("physical code violates operator Gauss")
    if len({(row["fermion_occupation_mask"], row["electric_flux"]) for row in basis}) != 4:
        raise ArithmeticError("the physical code is not an isometry")
    # CAR on the full four-mode Fock space, not just the chosen two-particle sector.
    car_cases = 0
    for occupation in range(16):
        for i in range(4):
            for j in range(4):
                for create_j in (False, True):
                    coefficients = {}
                    for word in (((j, create_j), (i, False)),
                                 ((i, False), (j, create_j))):
                        image = _ladder_word(occupation, word)
                        if image is not None:
                            sign, output = image
                            coefficients[output] = coefficients.get(output, 0)+sign
                    coefficients = {key: value for key, value in coefficients.items() if value}
                    expected = {occupation: 1} if create_j and i == j else {}
                    if coefficients != expected:
                        raise ArithmeticError("canonical anticommutation relation failed")
                    car_cases += 1
    ambient_gauss_cases = 0
    for occupation in range(16):
        for flux in (-2, -1, 0, 1, 2):
            for channel in range(2):
                for _, target, new_flux in dressed_hopping_images(occupation, flux, channel):
                    if gauss_values(target, new_flux) != gauss_values(occupation, flux):
                        raise ArithmeticError("dressed hop changed a Gauss charge")
                    ambient_gauss_cases += 1
    ident = identity(4)
    for unitary in (ideal_hopping(), ideal_electric(), ideal_echo()):
        if matmul(dagger(unitary), unitary) != ident:
            raise ArithmeticError("exact code operation is not unitary")
    h, v = hopping_hamiltonian(), electric_energy()
    if matmul(h, v) == matmul(v, h):
        raise ArithmeticError("electric and matter factors unexpectedly commute")
    for bit in (0, 1):
        p, q = (occupation_projector(bit, outcome) for outcome in (0, 1))
        if (matrix_add(p, q) != ident or matmul(p, p) != p or matmul(q, q) != q
                or any(cell != ZERO for row in matmul(p, q) for cell in row)):
            raise ArithmeticError("binary readout is not a complete projective instrument")
    if branch_probability(ideal_echo(), 0, 0, 1) != Q(576, 3125):
        raise ArithmeticError("quantum electric feedback witness changed")
    return {"CAR_cases": car_cases, "ambient_Gauss_hopping_images": ambient_gauss_cases,
            "physical_code_columns": 4, "exact_unitaries": 3,
            "hopping_and_electric_energy_do_not_commute": True,
            "binary_projector_completeness_exact": True}


def audit():
    """Replay exact identities and every native unitary column.

    The numerical channel check includes all 16 matrix units; its residual
    is a numerical diagnostic, not an outward floating-point error proof.
    Exact ideal identities and the analytic source-control identities have
    their separate roles. Complete transfer failures remain declared by the
    inherited compiler, even though ideal valid-code input never leaks.
    """
    exact = _exact_checks()
    tape = source_echo_word()
    actual = compiler.execute(tape, 2)
    target = numerical_matrix(ideal_echo())
    channel_residual = 0.
    for i in range(4):
        for j in range(4):
            output = np.outer(actual[:, i], actual[:, j].conj())
            expected = np.outer(target[:, i], target[:, j].conj())
            channel_residual = max(channel_residual, float(np.linalg.norm(output-expected, 2)))
    residual = float(np.linalg.norm(actual-target, 2))
    if max(residual, channel_residual) > 2e-12:
        raise ArithmeticError("native code word disagrees with exact gauge-matter operation")
    source_paths = ("model/observer/source_scalar/compiler.py",
                    "code/source_selection_model/response.json",
                    "code/m1_source_realization/receipt.json")
    return {"physical_basis": physical_basis(), "charges": list(CHARGES),
            "mode_order": list(MODES), "logical_bit_order": "index=b0+2*b1; bit j is right occupation in channel j",
            "ambient_electric_link": "ell2(Z), U|ell>=|ell+1>; no cyclic cutoff",
            "exact_checks": exact, "ideal_echo": encode_matrix(ideal_echo()),
            "right_occupation_probabilities_by_input": [
                [str(branch_probability(ideal_echo(), initial, 0, outcome)) for outcome in (0, 1)]
                for initial in range(4)],
            "electric_disabled_probability_on_input_00": "0",
            "native_full_unitary_residual": residual,
            "native_all_16_matrix_units_channel_residual": channel_residual,
            "resources": compiler.resources(tape),
            "source_sha256": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                              for path in source_paths},
            "finite_interacting_code_executed": True,
            "continuum_limit_or_physical_identification_supplied": False}
