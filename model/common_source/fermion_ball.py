"""Regional vacuum states of the same calibrated free fermion action.

Each charged pair of all-left Weyl fields is mapped locally, by a canonical
particle-hole change, to h_D=[[-sigma.p,m],[m,sigma.p]]. The retained mass
packet and CKM frame are unchanged. An even, antiperiodic Fourier grid in a
three-dimensional box regulates this supplied free action. It is a spectral
cutoff, not a local chiral gauge regulator or a derived physical geometry.

The full filled negative-energy projector is restricted to position sites
inside a ball. For this finite quasifree CAR state, K=log((1-C_A)/C_A) and
S=-tr(C_A log C_A+(1-C_A)log(1-C_A)); see Peschel--Eisler,
https://arxiv.org/abs/0906.1663, section3.3. Restricting h_D before choosing
its ground state would give a different state. A positive-energy orbital
with occupation epsilon supplies a legitimate global state variation and
the regional first law dS/depsilon=psi_A^dagger K psi_A.

The independently fixed massless continuum comparison is
2*pi*integral_B (R^2-r^2)/(2R) T00, Casini--Huerta--Myers,
https://arxiv.org/abs/1102.0440, eq.(2.20)/(22). Its use for massive fields
is only a comparison, not an asserted modular identity. One common packet,
box, region, clock-energy convention and source grammar serve every species.
Neither masses nor an area/entropy coefficient are fit in this calculation.
Vacuum preparation, cutoff/region choice and physical service accuracy remain
premises. Native selection of this vacuum, an interacting continuum and joint
state--geometry stationarity are not established.
The regional algebra is the full free-matter CAR algebra. Dynamical gauge
constraints, gauge-invariant dressings and gauge-field edge contributions
to an interacting region are not supplied by this partial trace.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
import sympy as s

from tools import strict_json
from . import calibrated_matter as matter
from .chiral_candidate import gaussian_source as source
from .chiral_candidate import kinetic
from ..native_equilibration import ball_modular


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("fermion_ball_results.json")
MIN_GAP = 1e-11
NUMERIC_ADMISSION = 2e-9
SOURCE_PATHS = tuple(sorted(set((*matter.source_bindings(), *ball_modular.SOURCE_PATHS,
    *source.compiler.source_summary()["source_sha256"],
    "model/common_source/fermion_ball.py",
    "model/common_source/test_fermion_ball.py",
    "model/common_source/chiral_candidate/gaussian_source.py",
    "model/common_source/mass_response/model.py",
    "model/observer/source_scalar/compiler.py"))))


def _positive(value, name):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError(f"positive exact rational {name} required")
    if not F(1, 10**8) <= value <= 10**8:
        raise ValueError(f"{name} outside declared numerical domain")
    return F(value)


def _mass(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, F)):
        raise ValueError("nonnegative finite comparison mass required")
    try:
        value = float(value)
    except OverflowError as error:
        raise ValueError("mass outside declared numerical domain") from error
    if not math.isfinite(value) or not 0 <= value <= 1e8:
        raise ValueError("mass outside declared numerical domain")
    return value


def matrices():
    zero, eye = np.zeros((2, 2)), np.eye(2)
    alpha = np.array([np.block([[-sigma, zero], [zero, sigma]])
                      for sigma in kinetic.pauli_matrices()])
    beta = np.block([[zero, eye], [eye, zero]])
    return alpha, beta


def exact_attachment():
    """Local canonical embedding, with no fitted spinor phases or masses."""
    px, py, pz, mass = s.symbols("px py pz mass", real=True)
    sigma = s.Matrix([[pz, px-s.I*py], [px+s.I*py, -pz]])
    eps = s.Matrix([[0, 1], [-1, 0]])
    k = s.Matrix([[0, mass], [mass, 0]])
    h = -s.kronecker_product(sigma, s.eye(2))
    delta = s.kronecker_product(eps, k)
    bdg = h.row_join(delta).col_join(delta.conjugate().T.row_join(h.T))
    embed = s.zeros(4, 8)
    for row, col, sign in ((0, 0, 1), (1, 2, 1), (2, 7, 1), (3, 5, -1)):
        embed[row, col] = sign
    nambu_exchange = s.zeros(4).row_join(s.eye(4)).col_join(
        s.eye(4).row_join(s.zeros(4)))
    dirac = (-sigma).row_join(mass*s.eye(2)).col_join(
        (mass*s.eye(2)).row_join(sigma))
    if (embed*embed.T != s.eye(4) or embed*nambu_exchange*embed.T != s.zeros(4)
            or embed*bdg != dirac*embed):
        raise ArithmeticError("local all-left/Dirac canonical attachment failed")
    charge = s.diag(1, -1, 1, -1, -1, 1, -1, 1)
    if embed*charge != embed:
        raise ArithmeticError("particle-hole map does not preserve the physical charge convention")
    if dirac**2 != (px*px+py*py+pz*pz+mass*mass)*s.eye(4):
        if s.simplify(dirac**2-(px*px+py*py+pz*pz+mass*mass)*s.eye(4)) != s.zeros(4):
            raise ArithmeticError("same-action Dirac shell failed")
    return {"local_Dirac_components_in_all_left_Nambu_order": ["chi_L_up", "chi_L_down", "chi_Rc_dagger_down", "-chi_Rc_dagger_up"],
        "intertwining_and_CAR_identity_exact": True,
        "both_CAR_anticommutators_checked": True,
        "Dirac_components_have_same_physical_charge_exactly": True,
        "symbol_square_exact": "(|p|^2+m^2) I4",
        "particle_hole_change_is_region_local": True,
        "CKM_and_singlet_frames_are_region_local": True,
        "regional_entropies_unchanged_by_same_weak_mass_frame": True,
        "charged_Dirac_species_including_color": 21,
        "neutral_massless_Weyl_species": 3,
        "total_complex_CAR_modes_per_position_site": 90,
        "optional_Weinberg_masses_included": False}


def spinors(momentum, mass):
    """Deterministic analytic negative/positive frames; no eigensolver gauge."""
    p, m = kinetic._momentum(momentum), _mass(mass)
    energy = math.sqrt(sum(x*x for x in p)+m*m)
    if energy == 0:
        raise ValueError("zero-energy vacuum selection is unresolved")
    sigma = sum((x*v for x, v in zip(p, kinetic.pauli_matrices(), strict=True)),
                np.zeros((2, 2), complex))
    a = math.sqrt((energy+m)/(2*energy))
    b = sigma/math.sqrt(2*energy*(energy+m))
    eye = np.eye(2)
    change = np.block([[eye, eye], [eye, -eye]])/math.sqrt(2)
    negative = change @ np.vstack((b, a*eye))
    positive = change @ np.vstack((a*eye, -b))
    return energy, negative, positive


def full_action_attachment():
    """Check every charged block of the actual complex180-dimensional action."""
    momentum = (0.2, 0.3, 0.4)
    weak = matter.weyl_symbol(momentum)
    w = np.asarray(matter.mass_basis_change().evalf(20), complex)
    spin_frame = np.kron(np.eye(2), w)
    zeros = np.zeros_like(spin_frame)
    frame = np.block([[spin_frame, zeros], [zeros, spin_frame.conj()]])
    mass_frame = frame.conj().T@weak@frame
    alpha, beta = matrices()
    errors = []
    for channel in matter.CHANNELS:
        coordinates = matter.mass_coordinates(channel)
        for left, right in matter._oriented_components(channel)[1].todok():
            for family, mass in enumerate(coordinates):
                a, b = 3*left+family, 3*right+family
                embed = np.zeros((4, 180), complex)
                embed[0, a] = embed[1, 45+a] = embed[2, 135+b] = 1.
                embed[3, 90+b] = -1.
                target = np.einsum("i,ijk->jk", momentum, alpha)+float(mass)*beta
                errors.append(float(np.linalg.norm(embed@mass_frame-target@embed, 2)))
    if len(errors) != 21 or max(errors) > NUMERIC_ADMISSION:
        raise ArithmeticError("regional Dirac blocks do not attach to the actual calibrated action")
    return {"actual_weak_BdG_dimension": 180, "charged_color_family_blocks_checked": len(errors),
        "maximum_intertwining_error": max(errors), "nonaxial_momentum": list(momentum),
        "same_empirical_packet_and_exact_CKM_frame": True}


def vacuum_grid(n, length, mass):
    if type(n) is not int or n % 2 or not 2 <= n <= 40:
        raise ValueError("even grid size2..40 required")
    length, mass = _positive(length, "box length"), _mass(mass)
    return _vacuum_grid(n, length, mass)


@lru_cache(maxsize=32)
def _vacuum_grid(n, length, mass):
    wave = 2*math.pi*(np.fft.fftfreq(n)*n+.5)/float(length)
    momenta = np.stack(np.meshgrid(wave, wave, wave, indexing="ij"), axis=-1)
    alpha, beta = matrices()
    energies = np.sqrt(np.sum(momenta*momenta, axis=-1)+mass*mass)
    h = np.einsum("...i,ijk->...jk", momenta, alpha)+mass*beta
    projector = (np.eye(4)-h/energies[..., None, None])/2
    residual = float(np.max(abs(projector@projector-projector)))
    if residual > NUMERIC_ADMISSION:
        raise ArithmeticError("full vacuum projector failed before restriction")
    kernel = np.fft.ifftn(projector, axes=(0, 1, 2))
    kernel = np.frombuffer(kernel.tobytes(), complex).reshape(kernel.shape)
    return kernel, residual


def region(n, length, radius):
    if type(n) is not int or n % 2 or not 2 <= n <= 40:
        raise ValueError("even grid size2..40 required")
    length, radius = _positive(length, "box length"), _positive(radius, "radius")
    if 2*radius >= length:
        raise ValueError("ball must fit strictly inside the box")
    # The strict geometric test is rational, not a floating boundary choice.
    labels = np.array(list(np.ndindex((n, n, n))), dtype=np.int64)
    doubled = 2*labels-(n-1)
    square = np.sum(doubled*doubled, axis=1)
    limit = (2*n*radius/length)**2
    keep = square*limit.denominator < limit.numerator
    labels = labels[keep]
    if not len(labels):
        raise ValueError("the chosen finite ball contains no grid sites")
    positions = float(length)/(2*n)*(2*labels-(n-1))
    return labels, positions


def restricted_covariance(n, length, mass, labels):
    kernel, _ = vacuum_grid(n, length, mass)
    if np.ma.isMaskedArray(labels):
        raise ValueError("masked site menus are not a region")
    labels = np.asarray(labels)
    if (labels.dtype.kind not in "iu" or labels.ndim != 2 or labels.shape[1] != 3
            or not len(labels) or np.any(labels < 0) or np.any(labels >= n)
            or len(np.unique(labels, axis=0)) != len(labels)):
        raise ValueError("distinct original grid sites required")
    labels = labels.astype(np.int64, copy=False)
    delta = labels[:, None, :]-labels[None, :, :]
    raw = kernel[delta[:, :, 0] % n, delta[:, :, 1] % n, delta[:, :, 2] % n]
    raw = raw*np.exp(1j*math.pi*np.sum(delta, axis=-1)/n)[:, :, None, None]
    return raw.transpose(0, 2, 1, 3).reshape(4*len(labels), 4*len(labels))


def modular(covariance):
    value = np.asarray(covariance)
    if (np.ma.isMaskedArray(covariance) or value.dtype.kind not in "fc"
            or value.ndim != 2 or value.shape[0] != value.shape[1] or not len(value)
            or not np.isfinite(value).all()
            or np.linalg.norm(value-value.conj().T, ord=np.inf) > NUMERIC_ADMISSION):
        raise ValueError("finite Hermitian covariance required without repair")
    eigenvalues, vectors = eigh(value, driver="evr")
    residual = float(np.linalg.norm(value@vectors-vectors*eigenvalues, ord=np.inf))
    gap = float(min(eigenvalues[0], 1-eigenvalues[-1]))
    if gap <= MIN_GAP or residual > min(NUMERIC_ADMISSION, gap/20):
        raise ValueError("regional state is unresolved near purity; no clipping allowed")
    logs = np.log1p(-eigenvalues)-np.log(eigenvalues)
    entropy = float(np.sum(-eigenvalues*np.log(eigenvalues)
                           -(1-eigenvalues)*np.log1p(-eigenvalues)))
    return entropy, (vectors*logs)@vectors.conj().T, gap, residual


def positive_orbital(n, length, mass, positions):
    p = np.full(3, math.pi/float(length))
    energy, _, positive = spinors(p, mass)
    psi = (np.exp(1j*positions@p)[:, None]*positive[:, 0]/n**1.5).ravel()
    return energy, psi


def ball(n=20, length=F(6), radius=F(1), mass=0., *, secant=False):
    if type(secant) is not bool:
        raise ValueError("explicit Boolean entropy control required")
    mass = _mass(mass)
    labels, positions = region(n, length, radius)
    covariance = restricted_covariance(n, length, mass, labels)
    entropy, modular_form, gap, residual = modular(covariance)
    energy, psi = positive_orbital(n, length, mass, positions)
    derivative = float(np.vdot(psi, modular_form@psi).real)
    stress = float(2*math.pi*energy*np.sum((float(radius)**2
                          -np.sum(positions*positions, axis=1))/(2*float(radius)))/n**3)
    result = {"n": n, "box_length": str(length), "ball_radius": str(radius),
        "dimensionless_action_mass": mass, "regional_sites": len(labels),
        "regional_complex_CAR_modes_per_Dirac_species": len(covariance),
        "purity_gap": gap, "eigensystem_residual": residual,
        "entropy": entropy, "positive_orbital_energy": energy,
        "regional_particle_weight": float(np.vdot(psi, psi).real),
        "actual_modular_entropy_derivative": derivative,
        "conformal_weighted_stress_comparison": stress,
        "relative_modular_minus_conformal_comparison": derivative/stress-1,
        "massive_conformal_identity_asserted": False}
    if secant:
        step = 1e-3
        direction = np.outer(psi, psi.conj())
        perturbed = covariance+step*direction
        after, _, perturbed_gap, _ = modular(perturbed)
        measured = (after-entropy)/step
        lower_gap = min(gap, perturbed_gap)-NUMERIC_ADMISSION
        if lower_gap <= 0:
            raise ValueError("entropy secant gap unresolved")
        curvature = float(np.vdot(psi, psi).real**2/(lower_gap*(1-lower_gap)))
        if not derivative-step*curvature/2-2e-7 <= measured <= derivative+2e-7:
            raise ArithmeticError("independent entropy secant violates first-law enclosure")
        result["entropy_first_law_control"] = {"occupation_step": step,
            "independent_entropy_difference_over_step": measured,
            "analytic_secant_lower": derivative-step*curvature/2,
            "analytic_secant_upper": derivative,
            "floating_secant_admission": 2e-7,
            "bound_conditional_on_numerical_spectral_admission": True}
    return result


def nested_relative_entropy_control():
    """Two nested restrictions of one global vacuum/excited-orbital state."""
    n, length, mass, occupation = 12, F(6), 0.7, 0.2
    rows = []
    for radius in (F(3, 4), F(1)):
        labels, positions = region(n, length, radius)
        reference = restricted_covariance(n, length, mass, labels)
        entropy, modular_form, _, _ = modular(reference)
        _, psi = positive_orbital(n, length, mass, positions)
        changed = reference+occupation*np.outer(psi, psi.conj())
        changed_entropy, _, _, _ = modular(changed)
        delta_k = occupation*float(np.vdot(psi, modular_form@psi).real)
        relative = delta_k-(changed_entropy-entropy)
        rows.append({"radius": str(radius), "sites": len(labels),
                     "relative_entropy": relative, "modular_expectation_change": delta_k,
                     "entropy_change": changed_entropy-entropy})
    if not -2e-10 <= rows[0]["relative_entropy"] <= rows[1]["relative_entropy"]+2e-10:
        raise ArithmeticError("nested restrictions violate relative-entropy monotonicity")
    return {"same_global_state_and_orbital": True, "excited_orbital_occupation": occupation,
            "regions": rows, "no_separate_regional_state_preparation": True}


def direct_nambu_covariance(n, length, mass, labels):
    """Independent full all-left BdG vacuum, for a small source control."""
    wave = 2*math.pi*(np.fft.fftfreq(n)*n+.5)/float(length)
    target = np.empty((n, n, n, 8, 8), complex)
    for index in np.ndindex((n, n, n)):
        p = wave[list(index)]
        h = kinetic.bdg_symbol(p, np.array([[0., mass], [mass, 0.]]), speed=1.)
        energy = math.sqrt(float(np.dot(p, p))+mass*mass)
        target[index] = (np.eye(8)+h/energy)/2
    kernel = np.fft.ifftn(target, axes=(0, 1, 2))
    delta = labels[:, None, :]-labels[None, :, :]
    blocks = kernel[delta[:, :, 0] % n, delta[:, :, 1] % n, delta[:, :, 2] % n]
    blocks = blocks*np.exp(1j*math.pi*np.sum(delta, axis=-1)/n)[:, :, None, None]
    nambu = blocks.transpose(0, 2, 1, 3).reshape(8*len(labels), -1)
    local = np.zeros((8, 8), complex)
    for j in range(4):
        local[2*j, j] = local[2*j, 4+j] = 1.
        local[2*j+1, j], local[2*j+1, 4+j] = -1j, 1j
    change = np.kron(np.eye(len(labels)), local)
    gamma = 1j*(change@nambu@change.conj().T-np.eye(8*len(labels)))
    if np.max(abs(gamma.imag)) > NUMERIC_ADMISSION:
        raise ArithmeticError("all-left Majorana covariance is not real")
    return gamma.real


def source_preparation(n=2, length=F(6), mass=1.):
    """Prepare a finite Slater vacuum by actual four-column source gates.

    Adjacent complex Givens rotations factor an analytic Fourier/spinor
    frame. Its angles follow the supplied action, not regional entropy data.
    The full Fock operator is Gaussian; no dense2**(4*n**3) state is used.
    """
    if n != 2:
        raise ValueError("literal preparation control fixes the2x2x2 grid")
    length, mass = _positive(length, "box length"), _mass(mass)
    labels = np.array(list(np.ndindex((n, n, n))))
    wave = 2*math.pi*(np.fft.fftfreq(n)*n+.5)/float(length)
    momenta = np.array(list(np.ndindex((n, n, n))))
    momenta = wave[momenta]
    modes = 4*n**3
    unitary = np.empty((modes, modes), complex)
    for j, p in enumerate(momenta):
        _, negative, positive = spinors(p, mass)
        phase = np.exp(1j*(float(length)/n*labels)@p)/math.sqrt(n**3)
        for sign, frame in enumerate((negative, positive)):
            unitary[:, sign*modes//2+2*j:sign*modes//2+2*j+2] = (
                phase[:, None, None]*frame).reshape(modes, 2)
    work = unitary.copy()
    calls = []
    for column in range(modes):
        for row in range(modes-1, column, -1):
            a, b = work[row-1, column], work[row, column]
            norm = math.hypot(abs(a), abs(b))
            cosine = abs(a)/norm if norm else 1.
            sine = (a/abs(a) if a else 1.)*b.conjugate()/norm if norm else 0.
            rotation = np.array([[cosine, sine], [-sine.conjugate(), cosine]])
            work[[row-1, row]] = rotation@work[[row-1, row]]
            coefficient = 1j*sine/abs(sine) if sine else 1j
            angle = math.atan2(abs(b), abs(a))
            calls.append((row-1, row, coefficient, angle))
    qr_error = float(np.max(abs(work-np.diag(np.diag(work)))))
    if qr_error > NUMERIC_ADMISSION:
        raise ArithmeticError("analytic vacuum frame Givens factorization failed")
    transform = np.eye(2*modes)
    onebody = np.eye(modes, dtype=complex)
    resources = {key: 0 for key in source._ADDITIVE_RESOURCES}
    maximum_local = 0.
    blocks = 0

    def apply(gate, sites, *, number_conserving=False):
        nonlocal maximum_local, blocks
        axes = [2*j+k for j in sites for k in (0, 1)]
        local = gate.orthogonal[:len(axes), :len(axes)]
        transform[axes] = local@transform[axes]
        if number_conserving:
            orbital = (np.array([[gate.fock_unitary[1, 1]]]) if len(sites) == 1
                       else gate.fock_unitary[np.ix_((1, 2), (1, 2))])
            onebody[list(sites)] = orbital@onebody[list(sites)]
        for name in resources:
            resources[name] += gate.resources[name]
        maximum_local = max(maximum_local, gate.diagnostics["majorana_linearity_residual"])
        blocks += 1

    for mode in range(0, modes//2, 2):
        apply(source.local_gate("pairing", 1., math.pi/2), (mode, mode+1))
    for mode, phase in enumerate(np.angle(np.diag(work))):
        apply(source.local_gate("onsite", 1., -float(phase)), (mode,), number_conserving=True)
    for left, right, coefficient, angle in reversed(calls):
        apply(source.local_gate("hopping", coefficient, -angle), (left, right), number_conserving=True)
    gamma = transform@source.vacuum_covariance(modes)@transform.T
    actual = ((gamma[0::2, 1::2]+np.eye(modes))/2+1j*gamma[0::2, 0::2]/2)
    expected = restricted_covariance(n, length, mass, labels)
    dirac_orthogonality = float(np.linalg.norm(transform@transform.T-np.eye(2*modes), 2))
    # Execute the local particle-hole map, including both spin phases and
    # the fermionic swap into the original spin-outer all-left ordering.
    for cell in range(n**3):
        first = 4*cell
        apply(source.local_gate("pairing", 1., math.pi/2), (first+2, first+3))
        for mode in (first+2, first+3):
            apply(source.local_gate("onsite", 1., math.pi/2), (mode,))
        apply(source.local_gate("fswap"), (first+1, first+2))
    all_left = transform@source.vacuum_covariance(modes)@transform.T
    independent = direct_nambu_covariance(n, length, mass, labels)
    errors = {"orbital_frame_operator_error": float(np.linalg.norm(onebody-unitary, 2)),
        "vacuum_covariance_operator_error": float(np.linalg.norm(actual-expected, 2)),
        "Dirac_Majorana_orthogonality_error": dirac_orthogonality,
        "all_left_Majorana_orthogonality_error": float(np.linalg.norm(transform@transform.T-np.eye(2*modes), 2)),
        "independent_all_left_BdG_vacuum_error": float(np.linalg.norm(all_left-independent, 2)),
        "maximum_local_gaussian_linearity_error": maximum_local}
    if max(errors.values()) > NUMERIC_ADMISSION:
        raise ArithmeticError(f"actual source vacuum preparation failed: {errors}")
    # Every emitted elementary angle is at most pi in absolute value.
    # A fixed envelope avoids attaching physical content to the arbitrary
    # QR phases of numerically empty elimination entries on different hosts.
    if resources["integrated_processor_norm_upper"] > math.pi*resources["native_pulses"]:
        raise ArithmeticError("uniform source-angle envelope failed")
    costs = {key: value for key, value in resources.items()
             if key not in ("zero_angle_requests", "integrated_processor_norm_upper")}
    return {"grid": [n, n, n], "complex_CAR_modes": modes, "filled_negative_orbitals": modes//2,
        "mass": float(mass), "all_four_local_Fock_input_columns_executed": True,
        "local_particle_hole_phases_and_FSWAP_executed": True,
        "adjacent_Givens_calls_including_zero_angles": len(calls),
        "source_blocks_including_blank_to_filled_pairs": blocks,
        "numerical_errors": errors, "actual_resource_counts": costs,
        "integrated_processor_norm_upper_conservative": math.pi*resources["native_pulses"],
        "zero_angle_calls_charged_at_full_scheduled_cost": True,
        "reference_global_Fock_dimension": str(1 << modes),
        "dense_global_Fock_execution": False,
        "native_passive_vacuum_selection_proved": False,
        "physical_service_error_certified": False}


def build():
    # One length conversion for every field; no species-dependent cutoff.
    radius_in_inverse_GeV = F(1, 100)
    rows = []
    for channel in matter.CHANNELS:
        for family, exact in enumerate(matter.mass_coordinates(channel)):
            mass = float(exact)
            rows.append({"channel": channel, "family": family,
                "color_multiplicity": 1 if channel == "charged_lepton" else 3,
                "supplied_running_mass_GeV": str(exact),
                "ball": ball(mass=mass*float(radius_in_inverse_GeV),
                             secant=channel == "up" and family == 2)})
    massless = ball(secant=True)
    total = sum(row["color_multiplicity"]*row["ball"]["entropy"] for row in rows)+1.5*massless["entropy"]
    return {"schema": "oph.fermion_ball.v1", "classification": "same-calibrated-free-action regional vacuum calculation",
        "exact_attachment": exact_attachment(), "full_calibrated_action_check": full_action_attachment(),
        "constant_native_comparison_frame": ball_modular.proper_frame(),
        "common_regulator": {"grid": [20, 20, 20], "boundary": "antiperiodic in all three spatial directions",
            "ball_radius_in_inverse_GeV": str(radius_in_inverse_GeV), "box_side_in_ball_radii": "6",
            "same_cutoff_and_radius_for_all_species": True, "local_chiral_gauge_regulator": False},
        "charged_species": rows, "massless_Dirac_control": massless,
        "full_90_component_regional_entropy": total,
        "massless_Weyl_entropy_is_half_Dirac_by_particle_hole_pairing": True,
        "grid_controls": [ball(n=n, mass=mass) for n in (12, 16, 24, 28) for mass in (0., float(matter.mass_coordinates("up")[2])/100)],
        "fixed_spacing_box_controls": [ball(n=n, length=F(3*n, 10), mass=mass)
            for n in (16, 24) for mass in (0., float(matter.mass_coordinates("up")[2])/100)],
        "same_state_nested_region_control": nested_relative_entropy_control(),
        "source_vacuum_preparation": source_preparation(mass=float(matter.mass_coordinates("up")[2])/100),
        "scope": {"same_calibrated_action_used_for_spectrum_and_entropy": True,
            "state_preparation_and_spectral_cutoff_supplied": True,
            "large_grid_native_vacuum_preparation_executed": False,
            "regional_entropy_estimated_from_detector_records": False,
            "regional_algebra": "full free-matter CAR; no interacting gauge constraints or edge algebra supplied",
            "same_native_continuum_geometry_identified": False,
            "interacting_gauge_Higgs_vacuum_executed": False,
            "massive_conformal_modular_formula_assumed": False,
            "renormalized_area_density_or_Newton_constant_derived": False,
            "joint_state_geometry_entropy_stationarity_derived": False,
            "numerical_diagnostics_are_outward_certificates": False,
            "mass_inputs_removed_or_derived": False},
        "source_bindings": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}}


def compare(old, new, path="root"):
    if isinstance(old, dict):
        if not isinstance(new, dict) or old.keys() != new.keys():
            raise ValueError(f"receipt keys differ at {path}")
        for key in old:
            compare(old[key], new[key], path+"."+key)
    elif isinstance(old, list):
        if not isinstance(new, list) or len(old) != len(new):
            raise ValueError(f"receipt menu differs at {path}")
        for i, (a, b) in enumerate(zip(old, new, strict=True)):
            compare(a, b, path+f"[{i}]")
    elif type(old) is float:
        exact = path.rsplit(".", 1)[-1] in {"mass", "dimensionless_action_mass", "occupation_step", "floating_secant_admission"}
        tolerance = (2e-7 if "independent_entropy_difference_over_step" in path else
                     5e-12 if path.endswith("purity_gap") else 2e-10)
        if type(new) is not float or not math.isfinite(old) or not math.isfinite(new) or (old != new if exact else not math.isclose(old, new, rel_tol=2e-9, abs_tol=tolerance)):
            raise ValueError(f"numerical replay differs at {path}")
    elif type(old) is not type(new) or old != new:
        raise ValueError(f"exact replay differs at {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", type=Path)
    group.add_argument("--write", type=Path)
    args = parser.parse_args()
    result = build()
    if args.check:
        compare(strict_json.load(args.check), result)
        print("Same-action fermion ball vacuum, source preparation and modular response replay passed")
    else:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
        print(f"Wrote {args.write}")


if __name__ == "__main__":
    main()
