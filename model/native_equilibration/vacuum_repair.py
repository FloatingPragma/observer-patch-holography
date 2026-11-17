"""A native geometry/vacuum repair, its preparation bound and work budget.

Compatibility is specified independently: geometry is uniform on the actual
two-state native seam orbit (one spectator particle), and the conditional
field has small excess energy above the ground state of its own Hamiltonian.
For geometry r>0, the EXISTING native angular matrix supplies

 H_r = r*p.p/2 + q.(K_radial+r**3*L_native).q/(2*r).

Thus density=1/r and angular_scale=r**3. The action-derived proper volume
is constant branchwise, while the two-face area is2/r. A fixed reference
Fock basis is used for every r. The finite Hamiltonian is the actual
compression P_D H_r P_D: squaring an already truncated q or p is incorrect.

One fair native C/E event resamples the geometry bit. Its committed record
selects a source-compiled conditional vacuum preparation and ordinary SWAP
with fresh blank field registers; the displaced input is exported to an
environment. All records and service requirements are counted. This is an
active, nonunital repair apparatus with a fresh-blank premise, not an
autonomous physical cooling law or a selected gravitational dynamics.

For every returned binary preparation, exact rational Rayleigh enclosures
and the known FULL oscillator gap certify its overlap with the unrestricted
Gaussian vacuum:1-F <= (E-E0)/gap. This uses the variational principle, not
an assumed eigensolver accuracy. Native pulse execution is numerical; its
unchanged binary arrays, norm defect and source replay remain explicit.

The resulting joint CQ state is repaired and field-vacuum compatible, but
its regional entropy plus any positive BARE area coefficient is not
stationary along the fixed-volume family. The slab/anisotropic cutoff
limitations of vacuum_geometry remain; this is not a small-ball or
renormalized-Sgen no-go theorem.

There is an additional exact obstruction to making this repair autonomous:
switching geometry and resetting its vacuum consumes positive mean quench
work. For each common normal mode let a_g=mu_g*omega_g. Then

 W_cycle=(omega_++omega_-)/16*(a_+/a_-+a_-/a_+-2) >=0.

The native constant angular mode makes the sum strictly positive for any
distinct branch geometries. No static energy assigned to the geometry bit
pays this stationary cycle cost. Positive branch-dependent clock rescalings
preserve its strict sign. A reciprocal geometric energy/entropy response or
a retained work reservoir is needed; neither is supplied by naming a common
fixed point or by cooling the field. No curvature or stress is defined to
enforce an Einstein equation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
import sympy as s

from tools import strict_json
from ..observer import native
from ..observer.source_scalar import compiler
from ..calibrated_measurement.bounds import root_interval
from ..observer.scalar_response.experiment import canonical
from . import vacuum_geometry as vacuum


ROOT = Path(__file__).resolve().parents[2]
BANDS = ((s.S.Zero, 1), (5-s.sqrt(5), 3), (s.Integer(6), 5), (5+s.sqrt(5), 3))
TRACE_ADMISSION = F(1, 256)


def _geometry(value):
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or not F(1, 2) <= value <= 2):
        raise ValueError("an original Fraction geometry in[1/2,2] required")
    return value


def _domain(cutoff, band, geometry):
    if type(cutoff) is not int or cutoff not in (2, 4, 8):
        raise ValueError("a fixed Fock cutoff2,4 or8 required")
    if type(band) is not int or not 0 <= band < 4:
        raise ValueError("one of the four actual native angular bands required")
    return _geometry(geometry)


def _add(a, b): return a[0]+b[0], a[1]+b[1]
def _mul(a, b):
    corners = [x*y for x in a for y in b]
    return min(corners), max(corners)


@lru_cache(maxsize=128)
def _algebraic_interval(expression):
    """Directed rational arithmetic for the exact quadratic-radical data."""
    value = s.sympify(expression)
    if value.is_Rational:
        point = F(int(value.p), int(value.q))
        return point, point
    if value.is_Add:
        result = F(), F()
        for term in value.args:
            result = _add(result, _algebraic_interval(term))
        return result
    if value.is_Mul:
        result = F(1), F(1)
        for term in value.args:
            result = _mul(result, _algebraic_interval(term))
        return result
    if value.is_Pow:
        base = _algebraic_interval(value.base)
        exponent = value.exp
        if exponent in (s.Rational(1, 2), s.Rational(-1, 2)):
            if base[0] <= 0:
                raise ValueError("positive quadratic-radical base required")
            base = root_interval(base[0])[0], root_interval(base[1])[1]
            if exponent > 0:
                return base
            return 1/base[1], 1/base[0]
        if exponent.is_Integer:
            count = int(exponent)
            if count < 0:
                if base[0] <= 0 <= base[1]:
                    raise ValueError("nonzero exact algebraic denominator required")
                base, count = (1/base[1], 1/base[0]), -count
            result = F(1), F(1)
            for _ in range(count):
                result = _mul(result, base)
            return result
    raise ValueError("unsupported exact vacuum algebraic expression")


@lru_cache(maxsize=3)
def _terms(cutoff):
    """P q² P uses2N+I+a²+a†², retaining the top-level diagonal."""
    a = np.diag(np.sqrt(np.arange(1, cutoff)), 1)
    eye = np.eye(cutoff)
    number = np.diag(2*np.arange(cutoff)+1)
    square, coordinate = a@a+(a@a).T, a+a.T
    kron = lambda a, b, c: np.kron(np.kron(a, b), c)
    diagonal = kron(number, eye, eye)+kron(eye, number, eye)+kron(eye, eye, number)
    squeeze = kron(square, eye, eye)+kron(eye, square, eye)+kron(eye, eye, square)
    cross = kron(coordinate, coordinate, eye)+kron(eye, coordinate, coordinate)
    return diagonal, squeeze, cross


def _coefficient_expressions(band, r):
    lam = BANDS[band][0]
    reference_frequency = s.sqrt(8+lam)
    common = (8+lam*r**3)/(4*r*reference_frequency)
    return (r*reference_frequency/4+common,
            -r*reference_frequency/4+common, -2/(r*reference_frequency))


def coefficients(band, geometry):
    return _coefficient_expressions(band, s.Rational(_domain(2, band, geometry)))


def hamiltonian(cutoff=8, band=2, geometry=F(1)):
    geometry = _domain(cutoff, band, geometry)
    return sum(float(c)*term for c, term in zip(coefficients(band, geometry), _terms(cutoff), strict=True))


def reduced(state, cutoff):
    wave = np.asarray(state).reshape(cutoff, cutoff, cutoff).transpose(1, 0, 2).reshape(cutoff, cutoff**2)
    return wave@wave.conj().T


def entropy(state, cutoff):
    wave = np.asarray(state).reshape(cutoff, cutoff, cutoff).transpose(1, 0, 2).reshape(cutoff, cutoff**2)
    eigenvalues = np.linalg.svd(wave, compute_uv=False)**2
    return -math.fsum(float(p)*math.log(float(p)) for p in eigenvalues if p > 0)


def ground(cutoff=8, band=2, geometry=F(1)):
    H = hamiltonian(cutoff, band, geometry)
    values, vectors = eigh(H, subset_by_index=(0, 1))
    target = vectors[:, 0]
    if target[np.argmax(abs(target))] < 0:
        target = -target
    if values[1]-values[0] <= 0 or np.linalg.norm(H@target-values[0]*target) > 1e-10:
        raise ArithmeticError("finite Hamiltonian ground state is unresolved")
    return H, float(values[0]), target


def response(cutoff=4, band=2, geometry=F(1)):
    """Differentiate compatibility on the SAME Fock algebra, not its basis."""
    geometry = _domain(cutoff, band, geometry)
    H = hamiltonian(cutoff, band, geometry)
    values, vectors = eigh(H)
    psi = vectors[:, 0]
    r = s.Symbol("r", positive=True)
    derivative = [s.diff(c, r).subs(r, s.Rational(geometry)) for c in _coefficient_expressions(band, r)]
    Hp = sum(float(c)*term for c, term in zip(derivative, _terms(cutoff), strict=True))
    force = float(psi@Hp@psi)
    tangent = -vectors[:, 1:]@((vectors[:, 1:].T@Hp@psi)/(values[1:]-values[0]))
    wave = psi.reshape(cutoff, cutoff, cutoff).transpose(1, 0, 2).reshape(cutoff, cutoff**2)
    delta = tangent.reshape(cutoff, cutoff, cutoff).transpose(1, 0, 2).reshape(cutoff, cutoff**2)
    rho, drho = wave@wave.T, delta@wave.T+wave@delta.T
    U, singular, Vh = np.linalg.svd(wave, full_matrices=False)
    rates = np.diag(U.T@delta@Vh.T)
    # This is a two-sided PURE-STATE amplitude path. At a zero Schmidt
    # singular value its eigenvalue changes quadratically, and s*log(s²)
    # tends to zero. No finite first-law claim on arbitrary rank-changing
    # density tangents is inferred from this limit.
    slope = -math.fsum(2*float(x*y)*math.log(float(x*x))
                       for x, y in zip(singular, rates, strict=True) if x > 0)
    residual = np.linalg.norm((H-values[0]*np.eye(len(H)))@tangent+(Hp-force*np.eye(len(H)))@psi)
    if residual > 2e-10 or abs(float(psi@tangent)) > 1e-10 or abs(float(np.trace(drho))) > 1e-10:
        raise ArithmeticError("differentiated zero-excess-energy compatibility failed")
    return {"entropy_derivative": slope, "Hellmann_Feynman_force": force,
        "differentiated_ground_constraint_residual": float(residual),
        "regional_trace_tangent_residual": abs(float(np.trace(drho))),
        "minimum_regional_eigenvalue_diagnostic": float(min(singular)**2),
        "entropy_derivative_domain": "two-sided reduced pure-state amplitude path; s*log(s^2) zero limit"}


def _binary_norm(vector):
    return sum((F(float(z.real))**2+F(float(z.imag))**2 for z in vector), F())


def energy_certificate(vector, cutoff, band, geometry):
    """Exact binary-array Rayleigh enclosure; no floating eigensolver premise."""
    geometry = _domain(cutoff, band, geometry)
    raw = np.asarray(vector)
    if (np.ma.isMaskedArray(vector) or raw.dtype.kind not in "fc"
            or (raw.dtype.kind == "f" and raw.dtype.itemsize > 8)
            or (raw.dtype.kind == "c" and raw.dtype.itemsize > 16)
            or raw.shape != (cutoff**3,) or not np.isfinite(raw).all()):
        raise ValueError("one original finite real/complex source vector required")
    entries = tuple((F(float(z.real)), F(float(z.imag))) for z in raw)
    norm = _binary_norm(raw)
    if norm <= 0 or abs(norm-1) > F(1, 10**8):
        raise ValueError("original source norm must lie within its declared admission")
    diagonal, squeeze, cross = (F(), F()), (F(), F()), (F(), F())
    root_cache = {}
    def root(n):
        if n not in root_cache:
            root_cache[n] = root_interval(F(n))
        return root_cache[n]
    def pairing(i, j):
        return 2*(entries[i][0]*entries[j][0]+entries[i][1]*entries[j][1])
    def index(digits): return (digits[0]*cutoff+digits[1])*cutoff+digits[2]
    for digits in product(range(cutoff), repeat=3):
        i = index(digits)
        weight = (2*sum(digits)+3)*(entries[i][0]**2+entries[i][1]**2)
        diagonal = _add(diagonal, (weight, weight))
        for site in range(3):
            if digits[site]+2 < cutoff:
                other = list(digits); other[site] += 2
                term = pairing(i, index(other))
                squeeze = _add(squeeze, _mul((term, term), root((digits[site]+1)*(digits[site]+2))))
        for left in (0, 1):
            for a, b in product((-1, 1), repeat=2):
                other = list(digits); other[left] += a; other[left+1] += b
                if min(other) < 0 or max(other) >= cutoff:
                    continue
                j = index(other)
                if j <= i:
                    continue
                factor = max(digits[left], other[left])*max(digits[left+1], other[left+1])
                term = pairing(i, j)
                cross = _add(cross, _mul((term, term), root(factor)))
    energy = (F(), F())
    for c, term in zip(coefficients(band, geometry), (diagonal, squeeze, cross), strict=True):
        energy = _add(energy, _mul(_algebraic_interval(c), term))
    energy = energy[0]/norm, energy[1]/norm
    lam, r = BANDS[band][0], s.Rational(geometry)
    frequencies = tuple(s.sqrt(k+lam*r**3) for k in (8-4*s.sqrt(2), s.Integer(8), 8+4*s.sqrt(2)))
    exact_ground = _algebraic_interval(sum(frequencies)/2)
    gap = _algebraic_interval(frequencies[0])
    excess_upper = (energy[1]-exact_ground[0])/gap[0]
    if excess_upper < 0:
        raise ArithmeticError("Rayleigh energy fell below the full oscillator ground bound")
    return {"norm_squared_exact": norm, "Rayleigh_energy_interval": energy,
        "full_Gaussian_ground_energy_interval": exact_ground,
        "full_Hamiltonian_gap_lower": gap[0], "infidelity_upper": excess_upper,
        "normalization_or_probability_repair_performed": False}


def prepare(cutoff=8, band=2, geometry=F(1)):
    _, energy, target = ground(cutoff, band, geometry)
    word = compiler.compile_real_state(target)
    bits = (cutoff**3).bit_length()-1
    blank = np.zeros(cutoff**3); blank[0] = 1
    actual = compiler.execute(word, bits, blank)
    error = float(np.linalg.norm(actual-target))
    if error > 2e-9 or abs(float(np.vdot(actual, actual).real)-1) > 2e-9:
        raise ArithmeticError("native vacuum preparation exceeded its numerical execution admission")
    cert = energy_certificate(actual, cutoff, band, geometry)
    return actual, {"cutoff": cutoff, "band": band, "geometry": str(geometry),
        "native_pulses": compiler.resources(word)["native_pulses"],
        "source_state_error_to_finite_target": error, "finite_ground_energy_diagnostic": energy,
        "source_regional_entropy_diagnostic": entropy(actual, cutoff), "certificate": cert}


def native_geometry():
    """Actual total-two orbit; C and E each supply the same four branches."""
    states = []
    for bit in (0, 1):
        loads = [0]*12; loads[bit] = loads[5] = 1
        states.append(tuple(loads))
    result = []
    for law in ("C", "E"):
        rows, labels = [], set()
        for initial, coin in product(range(2), repeat=2):
            after, (_, current) = native.native(states[initial], (0, 1), coin)
            final = states.index(after)
            label = (coin, current if law == "C" else states[initial][0]-states[initial][1])
            if label in labels or final != 1-coin:
                raise ArithmeticError("native geometry event lost its complete input/output record")
            labels.add(label)
            rows.append({"initial": initial, "final": final, "coin": coin, "label": label,
                         "stationary_branch_probability": "1/4"})
        result.append({"law": law, "branches": rows,
            "exact_population_transition": [["1/2", "1/2"], ["1/2", "1/2"]],
            "complete_record_entropy_nats": math.log(4),
            "reduced_geometry_channel": "rho -> I2*Tr(rho)/2, including offdiagonal inputs"})
    return result


def reset_dilation(geometry=F(1), band=2):
    """Every input column of a literal three-qubit field/reset dilation."""
    _, _, target = ground(2, band, geometry)
    prep = compiler.remap(compiler.compile_real_state(target), (3, 4, 5))
    swap = np.eye(4)[:, [0, 2, 1, 3]]
    swap_word = compiler.compile_real_orthogonal(swap)
    source_swap = compiler.execute(swap_word, 2)
    word = list(prep)
    for bit in range(3):
        word += compiler.remap(swap_word, (bit, bit+3))
    columns = np.zeros((64, 8), complex); columns[:8] = np.eye(8)
    actual = compiler.execute(word, 6, columns)
    desired = np.zeros_like(actual)
    for initial in range(8):
        desired[8*initial:8*(initial+1), initial] = target
    error = float(np.linalg.norm(actual-desired, 2))
    if error > 2e-9:
        raise ArithmeticError("source reset dilation failed its complete-input comparison")
    return {"source_isometry": actual, "target_isometry": desired,
        "all_input_columns": 8, "native_pulses": compiler.resources(word)["native_pulses"],
        "ordinary_swap_native_pulses": compiler.resources(swap_word)["native_pulses"],
        "ordinary_swap_all_column_error": float(np.linalg.norm(source_swap-swap, 2)),
        "all_input_isometry_error": error,
        "isometry_residual": float(np.linalg.norm(actual.conj().T@actual-np.eye(8), 2)),
        "fresh_blank_qubits": 3, "discarded_environment_qubits": 3,
        "input_or_reference_information_destroyed_globally": False}


def joint_reset_control():
    """All16 input columns of native geometry plus conditional source reset."""
    resets = [reset_dilation(r) for r in (F(9, 10), F(11, 10))]
    source, target = np.zeros((512, 16), complex), np.zeros((512, 16), complex)
    for initial, final in product(range(2), repeat=2):
        record = 2*initial+final
        rows = slice((2*record+final)*64, (2*record+final+1)*64)
        columns = slice(8*initial, 8*(initial+1))
        source[rows, columns] = resets[final]["source_isometry"]/math.sqrt(2)
        target[rows, columns] = resets[final]["target_isometry"]/math.sqrt(2)
    error = float(np.linalg.norm(source-target, 2))
    if error > 3e-9:
        raise ArithmeticError("coupled native geometry/vacuum history lost an input column")
    return source, target, {"all_joint_input_columns": 16, "complete_native_record_branches": 4,
        "literal_joint_domain": {"Fock_cutoff": 2, "native_band_index": 2,
            "native_angular_eigenvalue": "6", "field_oscillators": 3,
            "entire36_oscillator_reset_executed": False},
        "source_isometry_error": error,
        "source_isometry_completeness_residual": float(np.linalg.norm(source.conj().T@source-np.eye(16), 2)),
        "reduced_joint_channel_exact_target": "rho_GF -> (sum_g |g><g| tensor |vac_g><vac_g|/2)*Tr(rho_GF)",
        "source_feedback_reads_only_committed_native_geometry_record": True,
        "fresh_field_blanks_required_per_attempt": 3,
        "native_geometry_laws_checked": ["C", "E"]}


def work(central=F(1), contrast=F(1, 10)):
    central = _geometry(central)
    if type(contrast) is not F or not 0 <= contrast < central or central+contrast > 2 or central-contrast < F(1, 2):
        raise ValueError("two distinct admitted positive branch geometries required")
    lo, hi = central-contrast, central+contrast
    rows, total = [], 0.
    for band, (lam, multiplicity) in enumerate(BANDS):
        for radial, kappa in enumerate((8-4*s.sqrt(2), s.Integer(8), 8+4*s.sqrt(2))):
            wm, wp = s.sqrt(kappa+lam*s.Rational(lo)**3), s.sqrt(kappa+lam*s.Rational(hi)**3)
            am, ap = wm/s.Rational(lo), wp/s.Rational(hi)
            expression = (wm+wp)/16*((ap-am)**2/(ap*am))
            value = float(expression)
            rows.append({"band": band, "radial": radial, "multiplicity": multiplicity,
                         "mean_quench_work_diagnostic": value})
            total += multiplicity*value
    if (contrast > 0 and total <= 0) or (contrast == 0 and total != 0):
        raise ArithmeticError("distinct native geometry vacua produced no quench work")
    return {"branch_geometries": [str(lo), str(hi)], "modes": rows,
        "total_mean_quench_work_diagnostic": total,
        "exact_per_mode_nonnegative_form": "(omega_minus+omega_plus)*(a_plus-a_minus)^2/(16*a_plus*a_minus)",
        "constant_native_angular_mode_supplies_strict_positive_lower_bound": contrast > 0,
        "positive_branch_clock_rescalings_preserve_nonnegativity_and_equality_condition": True,
        "any_static_geometry_potential_mean_change_at_stationary_uniform_law": "0",
        "work_parameter": "shared model-coordinate Hamiltonian; no calibrated joule or universal physical clock"}


def source_work_interval(prepared, cutoff=8):
    """Direct exact cycle work of the RETURNED binary state rays.

    Rayleigh quotients are scale-invariant comparisons; source arrays are
    never normalized or repaired. This bound does not infer work or entropy
    derivatives from trace-distance closeness to an ideal vacuum.
    """
    result = F(), F()
    for band, (_, multiplicity) in enumerate(BANDS):
        energies = {}
        for source, target in product((F(9, 10), F(11, 10)), repeat=2):
            vector = prepared[source, band]
            energies[source, target] = energy_certificate(vector, cutoff, band, target)["Rayleigh_energy_interval"]
        plus = _add(energies[F(9, 10), F(11, 10)], energies[F(11, 10), F(9, 10)])
        minus = _add(energies[F(9, 10), F(9, 10)], energies[F(11, 10), F(11, 10)])
        term = plus[0]-minus[1], plus[1]-minus[0]
        result = _add(result, (multiplicity*term[0]/4, multiplicity*term[1]/4))
    return result


def certify_prepared(prepared):
    """Certify the literal supplied arrays, independently of fresh replay."""
    expected = {(r, band) for r in (F(9, 10), F(11, 10)) for band in range(4)}
    if type(prepared) is not dict or set(prepared) != expected:
        raise ValueError("both complete four-band binary preparation families required")
    bounds = {}
    for geometry in (F(9, 10), F(11, 10)):
        infidelity, global_norm = F(), F(1)
        for band, (_, multiplicity) in enumerate(BANDS):
            certificate = energy_certificate(prepared[geometry, band], 8, band, geometry)
            infidelity += multiplicity*certificate["infidelity_upper"]
            global_norm *= certificate["norm_squared_exact"]**multiplicity
        trace_upper = abs(global_norm-1)/2+global_norm*root_interval(infidelity)[1]
        if trace_upper >= TRACE_ADMISSION:
            raise ValueError("literal supplied source field fails global Gaussian vacuum fidelity gate")
        bounds[str(geometry)] = {"global_trace_distance_upper": str(TRACE_ADMISSION),
            "proof": "1-product(F_band^multiplicity)<=sum(multiplicity*(Rayleigh-E0)/gap); raw trace distance<abs(norm-1)/2+norm*sqrt(sum)<1/256",
            "norm_defect_is_included_in_raw_trace_certificate": True}
    work_box = source_work_interval(prepared)
    if work_box[0] <= F(1, 1000):
        raise ValueError("literal supplied source states fail positive cycle-work gate")
    return bounds, work_box


def _json(value):
    if isinstance(value, dict): return {key: _json(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)): return [_json(child) for child in value]
    if type(value) is F: return str(value)
    return value


def run():
    vacuum.native_angular()
    arrays, preparations, response_rows, prepared = {}, [], [], {}
    for geometry in (F(9, 10), F(11, 10)):
        for band, (_, multiplicity) in enumerate(BANDS):
            actual, row = prepare(8, band, geometry)
            arrays[f"prepared_r{geometry.numerator}_{geometry.denominator}_band{band}"] = actual
            prepared[geometry, band] = actual
            # The strict certificate is fixed before execution; detailed
            # binary-ray enclosures are recomputed from retained arrays.
            row.pop("certificate")
            preparations.append(row)
            response_rows.append({"geometry": str(geometry), "band": band,
                                  "multiplicity": multiplicity, **response(8, band, geometry)})
    bounds, work_box = certify_prepared(prepared)
    slope = sum(row["multiplicity"]*row["entropy_derivative"]/2 for row in response_rows)
    area_slope = -float(1/F(9, 10)**2+1/F(11, 10)**2)
    if not slope < 0 or not area_slope < 0:
        raise ArithmeticError("compatible native field repair unexpectedly supplied bare area stationarity")
    dilation = reset_dilation(F(9, 10))
    arrays["reset_isometry"] = dilation.pop("source_isometry")
    arrays["reset_target_isometry"] = dilation.pop("target_isometry")
    joint_source, joint_target, joint_control = joint_reset_control()
    arrays["joint_reset_isometry"], arrays["joint_reset_target_isometry"] = joint_source, joint_target
    sources = set(vacuum.SOURCE_PATHS)|set(compiler.source_summary()["source_sha256"])|{
        "model/native_equilibration/vacuum_repair.py", "model/native_equilibration/test_vacuum_repair.py",
        "model/observer/source_scalar/compiler.py", "model/calibrated_measurement/bounds.py"}
    report = {"schema": "oph.native_equilibration.vacuum_repair.v1",
        "native_geometry": native_geometry(), "preparations": preparations,
        "full36_oscillator_vacuum_preparation_bounds": bounds, "reset_dilation": dilation,
        "joint_native_geometry_and_vacuum_reset": joint_control,
        "compatible_state_geometry_response": response_rows,
        "joint_repaired_family": {"independent_compatibility": "uniform native geometry; conditional small excess field energy for its own H_r",
            "region": "central radial slab times every native angular mode, with all four committed geometry records",
            "exact_fixed_volume_dictionary": "density=1/r,angular_scale=r^3,L=2; two-face area2/r",
            "regional_record_entropy_derivative": 0,
            "finite_cutoff8_regional_field_entropy_derivative": slope,
            "entropy_derivative_is_of_ideal_finite_Hamiltonian_ground_family": True,
            "discrete_source_preparation_fidelity_does_not_bound_entropy_derivatives": True,
            "independent_mean_two_face_area_derivative": area_slope,
            "any_positive_bare_area_coefficient_gives_stationarity": False,
            "state_response_is_computed_on_one_fixed_Fock_algebra": True},
        "stationary_geometry_cycle_work": work(),
        "literal_binary_source_cycle_work": {"exact_certified_lower_bound": "1/1000",
            "interval_midpoint_diagnostic": float(sum(work_box)/2),
            "proof": "exact Rayleigh intervals of the same returned vectors in both branch Hamiltonians",
            "input_arrays_modified_or_renormalized": False},
        "source_resources": {"executed_unique_cutoff8_preparations": 8,
            "executed_cutoff8_preparation_pulses": sum(row["native_pulses"] for row in preparations),
            "full_spectral_field_logical_qubits": 108,
            "required_fresh_blank_logical_qubits_per_full_field_reset": 108,
            "required_full_field_preparation_pulses_per_branch": 12*12291,
            "required_full_field_swap_pulses": 108*dilation["ordinary_swap_native_pulses"],
            "required_full_field_reset_pulses": 12*12291+108*dilation["ordinary_swap_native_pulses"],
            "required_load_export_service_events": 2*(12*12291+108*dilation["ordinary_swap_native_pulses"]),
            "full_field_reset_counts_are_tensor_product_interface_requirements": True,
            "proper_load_export_service_histories_executed": False},
        "scope": {"active_native_preparation_and_swap_dilation_executed": True,
            "binary_ray_Gaussian_fidelity_bound_uses_exact_Rayleigh_enclosures": True,
            "all36_mode_source_wavefunction_materialized": False,
            "full_field_preparation_is_declared_tensor_product_of_executed_band_factors": True,
            "fresh_blanks_and_entropy_export_are_required": True,
            "autonomous_geometry_cooling_or_reciprocal_backreaction_derived": False,
            "geometry_branch_is_comparison_control_not_a_derived_dynamic_native_chart": True,
            "literal_binary_word_family_entropy_derivative_error_bound_certified": False,
            "proper_isotropic_UV_limit_or_renormalized_generalized_entropy": False,
            "small_ball_Einstein_stationarity_or_physical_work_units_derived": False},
        "source_bindings": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(sources)}}
    return _json(report), arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--output", type=Path)
    modes.add_argument("--check", type=Path)
    args = parser.parse_args()
    report, arrays = run()
    path = args.output or args.check
    artifact = path.with_suffix(".npz")
    if args.output:
        if path.exists() or artifact.exists():
            parser.error("fresh receipt/artifact paths required")
        with artifact.open("xb") as stream:
            np.savez_compressed(stream, **arrays)
        report["artifact"] = {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
        path.write_text(canonical(report))
    else:
        retained = strict_json.load(path)
        if retained.pop("artifact", None) != {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}:
            raise ValueError("vacuum repair artifact custody changed")
        with np.load(artifact, allow_pickle=False) as saved:
            if len(saved.files) != len(set(saved.files)) or set(saved.files) != set(arrays):
                raise ValueError("vacuum repair array menu changed")
            for key, fresh in arrays.items():
                old = saved[key]
                if old.dtype != fresh.dtype or old.shape != fresh.shape or not np.isfinite(old).all() or not np.allclose(old, fresh, atol=2e-10, rtol=2e-10):
                    raise ValueError("vacuum repair fresh source output changed: "+key)
            prepared = {(r, band): saved[f"prepared_r{r.numerator}_{r.denominator}_band{band}"]
                        for r in (F(9, 10), F(11, 10)) for band in range(4)}
            certify_prepared(prepared)
            for row in report["preparations"]:
                r, band = F(row["geometry"]), row["band"]
                fresh = arrays[f"prepared_r{r.numerator}_{r.denominator}_band{band}"]
                if (np.linalg.norm(prepared[r, band]-fresh)+row["source_state_error_to_finite_target"] > 2e-9
                        or abs(_binary_norm(prepared[r, band])-1) > F(1, 500000000)):
                    raise ValueError("retained source preparation exceeds its unchanged numerical admission")
            for prefix, admission in (("reset", 2e-9), ("joint_reset", 3e-9)):
                if np.linalg.norm(saved[prefix+"_isometry"]-arrays[prefix+"_target_isometry"], 2) > admission:
                    raise ValueError("retained source reset exceeds all-column admission")
        def compare(a, b, location="root"):
            if type(a) is not type(b): raise ValueError("vacuum repair receipt type changed")
            if isinstance(a, dict):
                if a.keys() != b.keys(): raise ValueError("vacuum repair receipt fields changed")
                for key in a: compare(a[key], b[key], location+"."+key)
            elif isinstance(a, list):
                if len(a) != len(b): raise ValueError("vacuum repair receipt list changed")
                for index, (x, y) in enumerate(zip(a, b, strict=True)): compare(x, y, location+f"[{index}]")
            elif isinstance(a, float):
                if "minimum_regional_eigenvalue" in location and (a < 0 or b < 0):
                    raise ValueError("vacuum repair regional spectrum diagnostic changed sign")
                if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a, b, rel_tol=2e-8, abs_tol=2e-8):
                    raise ValueError("vacuum repair numerical replay changed")
            elif a != b: raise ValueError("vacuum repair exact replay changed")
        compare(retained, report)
        print("native vacuum preparation, joint compatibility and geometry work replay: PASS")


if __name__ == "__main__":
    main()
