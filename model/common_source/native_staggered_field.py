"""Native exchange echo and a controlled free staggered field in3+1 dimensions.

The existing correlation/exchange architecture supplies B=J(I-Swap) and
native occupation gaps. A signed coherent echo cancels its density interaction
on ALL Fock inputs. Added parity transport, quarter-turn phases, cubic geometry
and staggered onsite signs specify a new free architecture. The real-space
operator itself derives an eight-corner Clifford representation with two
Dirac tastes; no logical Dirac matrix or measured mass packet is supplied.

The negative-energy sea is half filled. Empty encoded occupation is not its
ground state. Geometry, filling preparation, inverse controls, coherent source
accuracy and the programmed coordinate remain construction premises. This is
not the unechoed interacting Hamiltonian or a chiral Standard Model theory.

Primary staggered reference: Susskind, Phys.Rev.D16,3031(1977), and
https://arxiv.org/html/2604.02078v1 sections2.1--2.2. Finite product-formula
bound: https://arxiv.org/abs/1912.08854 Proposition16.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import repair_field_attachment as attachment
from . import record_code
from .chiral_candidate import gaussian_source as gaussian
from .mass_response import history as instruments
from .halfline import pi_interval
from ..observer.source_scalar import compiler


ROOT = Path(__file__).resolve().parents[2]
SOURCE_GATE = 1e-8
MACRO_GATE = 1e-10
REPLAY_TOLERANCE = 2e-9
TIME = F(1)
STEPS = 32


def _freeze(value):
    value = np.asarray(value, complex)
    return np.frombuffer(value.tobytes(), complex).reshape(value.shape)


def _rational(value, *, positive=False):
    if not (type(value) in (int, F) or isinstance(value, s.Rational)) or isinstance(value, bool):
        raise ValueError("an original exact rational required")
    result = F(value)
    if result < 0 or (positive and not result):
        raise ValueError("a nonnegative exact model parameter required")
    return result


def native_parameters(sector=0):
    if type(sector) is not int or sector not in (0, 1, 2):
        raise ValueError("one of the three native carrier sectors required")
    p = attachment.logical_parameters()
    return {"mass": F(p["positive_occupation_gaps"][sector]),
            "exchange": F(p["exchange_strength"]),
            "coordinate_speed": 2*F(p["exchange_strength"]),
            "onsite_baseline": F(p["onsite_baselines"][sector])}


def _phase_word(angle, mode=0):
    return [compiler.pulse("I", angle/2, [mode]), compiler.pulse("Z", -angle/2, [mode])]


def echo_word(duration):
    """Z U_B(-t/2) Z U_B(t/2)=exp(+i J t*(XX+YY)/2)."""
    duration = _rational(duration)
    first = attachment.exchange_word(duration/2)
    z = _phase_word(math.pi)
    return list(first)+z+compiler.inverse(first)+z


def echo_certificate(duration=F(3, 5)):
    duration = _rational(duration)
    J = native_parameters()["exchange"]
    B = J*(s.eye(4)-attachment.carrier_swap(0, 1, 2))
    Z = s.diag(1, -1, 1, -1)
    difference = (B-Z*B*Z)/2
    hop = s.zeros(4); hop[1, 2] = hop[2, 1] = -s.Rational(J)
    if difference != hop or B*(Z*B*Z) != (Z*B*Z)*B:
        raise ArithmeticError("native exchange echo did not exactly cancel the density interaction")
    actual = _freeze(compiler.execute(echo_word(duration), 2))
    target = expm(-1j*float(duration)*np.asarray(hop, complex))
    error = float(np.linalg.norm(actual-target, 2))
    if error > MACRO_GATE:
        raise ArithmeticError("source exchange echo exceeded its all-Fock operator gate")
    return {"unitary": actual, "target": _freeze(target), "word": tuple(echo_word(duration)),
            "resources": compiler.resources(echo_word(duration)),
            "operator_error_diagnostic": error,
            "density_term_cancelled_on_double_occupation": True,
            "inverse_coherent_exchange_and_Z_controls_supplied": True}


def _coordinates(shape):
    if (type(shape) is not tuple or len(shape) != 3
            or any(type(n) is not int or not 2 <= n <= 8 for n in shape)):
        raise ValueError("three finite lattice sides from2 through8 required")
    return tuple(product(*(range(n) for n in shape)))


def _eta(x, axis):
    return (-1)**sum(x[:axis])


def _epsilon(x):
    return (-1)**sum(x)


def lattice_operator(shape=(4, 4, 4), *, spacing=F(1), sector=0, periodic=True):
    """Native-derived staggered one-body operator, oriented bonds retained."""
    coordinates, spacing = _coordinates(shape), _rational(spacing, positive=True)
    if type(periodic) is not bool or (periodic and any(n % 2 or n < 4 for n in shape)):
        raise ValueError("periodic central differences require even sides at least4;2-site periodic hopping cancels")
    index = {x: j for j, x in enumerate(coordinates)}
    p = native_parameters(sector)
    h = s.zeros(len(coordinates))
    for j, x in enumerate(coordinates):
        h[j, j] = s.Rational(p["mass"]*_epsilon(x))
        for axis in range(3):
            y = list(x); y[axis] += 1
            if y[axis] == shape[axis]:
                if not periodic:
                    continue
                y[axis] = 0
            k = index[tuple(y)]
            z = -s.I*s.Rational(p["exchange"]*_eta(x, axis)/spacing)
            h[j, k] += z
            h[k, j] += s.conjugate(z)
    return s.ImmutableMatrix(h)


@lru_cache(maxsize=1)
def clifford_matrices():
    """Derived from Fourier corner shifts by the three real-space signs."""
    corners = tuple(product((0, 1), repeat=3))
    index = {r: j for j, r in enumerate(corners)}
    alphas, beta = [], s.zeros(8)
    for axis in range(3):
        matrix = s.zeros(8)
        for j, r in enumerate(corners):
            out = tuple(1-v if k < axis else v for k, v in enumerate(r))
            matrix[index[out], j] = (-1)**r[axis]
        alphas.append(s.ImmutableMatrix(matrix))
    for j, r in enumerate(corners):
        beta[index[tuple(1-v for v in r)], j] = 1
    generators = (*alphas, s.ImmutableMatrix(beta))
    for i, left in enumerate(generators):
        for j, right in enumerate(generators):
            if left*right+right*left != (2*s.eye(8) if i == j else s.zeros(8)):
                raise ArithmeticError("the source staggered signs do not derive a Clifford representation")
    return tuple(generators)


def corner_symbol(momentum, *, spacing=F(1), sector=0):
    if type(momentum) is not tuple or len(momentum) != 3:
        raise ValueError("three exact rational momenta required")
    spacing = _rational(spacing, positive=True)
    momentum = tuple(_rational_signed(p) for p in momentum)
    p, matrices = native_parameters(sector), clifford_matrices()
    return s.ImmutableMatrix(s.Rational(p["mass"])*matrices[3]+sum(
        (s.Rational(p["coordinate_speed"]/spacing)*s.sin(s.Rational(spacing*k))*matrices[j]
         for j, k in enumerate(momentum)), s.zeros(8)))


def _rational_signed(value):
    if not (type(value) in (int, F) or isinstance(value, s.Rational)) or isinstance(value, bool):
        raise ValueError("an original exact rational momentum required")
    return F(value)


def continuum_contract(*, spacing=F(1, 4), momentum_cutoff=F(1), time=F(1), excitations=1):
    """Uniform corner-band error, with the full Fock scaling stated.

    |sin(ap)/a-p|<=a²|p|³/6 and Clifford anticommutation give
    ||h_a-h_D||<=c*a²*sqrt(sum|p_j|^6)/6. Unitary Duhamel and exterior
    products bound a fixed finite number of particles and holes. The whole
    growing Fock space is not covered by a fixed-excitation bound.
    """
    a, cutoff, t = (_rational(spacing, positive=True), _rational(momentum_cutoff), _rational(time))
    if type(excitations) is not int or not 0 <= excitations <= 1000000:
        raise ValueError("a declared finite number of particle/hole excitations required")
    if 2*a*cutoff >= pi_interval()[0]:
        raise ValueError("the declared compact band must lie strictly inside the first reduced Brillouin zone")
    # sqrt3<2 yields a simple directed rational componentwise-cutoff bound.
    delta = native_parameters()["coordinate_speed"]*a*a*cutoff**3/3
    smallest_mass = min(native_parameters(j)["mass"] for j in range(3))
    projector = min(F(1), delta/smallest_mass)
    return {"spacing": str(a), "componentwise_momentum_cutoff": str(cutoff),
            "source_coordinate_duration": str(t), "particle_plus_hole_excitations": excitations,
            "one_body_operator_error_upper": str(delta),
            "propagator_operator_error_upper": str(min(F(2), t*delta)),
            "fixed_excitation_Fock_error_upper": str(min(F(2), excitations*t*delta)),
            "excitation_identification": "match positive/negative spectral fibers of the two normal-ordered filled-sea representations; raw source vacuum phase is retained separately",
            "vacuum_projector_operator_error_upper_all_native_sectors": str(projector),
            "vacuum_time_ordered_two_point_band_error_upper": str(min(F(2), projector+t*delta)),
            "unit_norm_smearing_four_point_Wick_error_upper": str(min(F(2), 4*(projector+t*delta))),
            "trace_zero_identity_term_or_wavefunction_fit_introduced": False,
            "source_product_formula_and_native_service_errors_additional": True,
            "continuum_contract": "fixed reduced-zone compact momentum band and finite excitation number; a->0 at declared c=2J",
            "whole_growing_Fock_or_interacting_continuum_claimed": False,
            "infinite_volume_vacuum_vector_norm_convergence_claimed": False,
            "native_hopping_control_duration_scales_as": "1/a",
            "arbitrarily_fast_coherent_service_and_physical_clock_dictionary_derived": False}


def clifford_certificate():
    matrices = clifford_matrices()
    I, X, Y, Z = s.eye(2), s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)
    tastes = (s.kronecker_product(I, Z, Y), s.kronecker_product(Z, X, Y), s.kronecker_product(Z, Y, I))
    for i, left in enumerate(tastes):
        if any(left*g != g*left for g in matrices):
            raise ArithmeticError("the proposed taste operators do not commute with the source Clifford action")
        for j, right in enumerate(tastes):
            if left*right+right*left != (2*s.eye(8) if i == j else s.zeros(8)):
                raise ArithmeticError("two-fold taste algebra failed")
    chirality = s.I*matrices[0]*matrices[1]*matrices[2]
    if (chirality != chirality.H or chirality**2 != s.eye(8)
            or any(chirality*g != g*chirality for g in matrices[:3])
            or chirality*matrices[3]+matrices[3]*chirality != s.zeros(8)):
        raise ArithmeticError("massless chirality and mass-breaking identities failed")
    u, v, w, mass = s.symbols("u v w m", real=True)
    h = mass*matrices[3]+u*matrices[0]+v*matrices[1]+w*matrices[2]
    if s.expand(h*h) != (mass**2+u**2+v**2+w**2)*s.eye(8):
        raise ArithmeticError("source staggered dispersion failed")
    gamma = (matrices[3], *(matrices[3]*alpha for alpha in matrices[:3]))
    metric = (1, -1, -1, -1)
    for i, left in enumerate(gamma):
        for j, right in enumerate(gamma):
            if left*right+right*left != (2*metric[i]*s.eye(8) if i == j else s.zeros(8)):
                raise ArithmeticError("the derived continuum action lacks its Lorentzian Clifford signature")
    spin = tuple(-s.I*matrices[(i+1) % 3]*matrices[(i+2) % 3]/2 for i in range(3))
    if sum((S*S for S in spin), s.zeros(8)) != s.Rational(3, 4)*s.eye(8):
        raise ArithmeticError("the source continuum lift does not have spin one-half")
    for i in range(3):
        if spin[i]*spin[(i+1) % 3]-spin[(i+1) % 3]*spin[i] != s.I*spin[(i+2) % 3]:
            raise ArithmeticError("the derived continuum rotation algebra failed")
    return {"corner_dimension": 8, "Dirac_tastes_per_native_sector": 2,
            "positive_energy_multiplicity": 4, "negative_energy_multiplicity": 4,
            "derived_alpha_labels": ["ZII", "XZI", "XXZ"], "derived_beta_label": "XXX",
            "derived_continuum_gamma_signature": [1, -1, -1, -1], "spin_Casimir": "3/4",
            "commuting_taste_Pauli_labels": ["IZY", "ZXY", "ZYI"],
            "massless_chirality_eigenvalue_multiplicities": [[-1, 4], [1, 4]],
            "single_Weyl_or_SM_chiral_representation_selected": False,
            "dispersion_exact": "E²=m²+(2J/a)²*sum_j sin²(a*p_j)",
            "negative_modes_discarded_or_taste_root_taken": False}


def periodic_fourier_certificate():
    coordinates = np.array(_coordinates((4, 4, 4)), int)
    h = np.asarray(lattice_operator(), complex)
    p = native_parameters()
    matrices = tuple(np.asarray(g, complex) for g in clifford_matrices())
    rows = []
    for bits in product((0, 1), repeat=3):
        momentum = -math.pi/2*np.array(bits)
        corners = np.array(list(product((0, 1), repeat=3)))
        wave = np.exp(1j*coordinates@(momentum+math.pi*corners).T)/8
        target = float(p["mass"])*matrices[3]+sum((float(p["coordinate_speed"])*math.sin(k)*g
                     for k, g in zip(momentum, matrices[:3], strict=True)), np.zeros((8, 8), complex))
        error = float(np.linalg.norm(h@wave-wave@target, 2))
        if error > SOURCE_GATE:
            raise ArithmeticError("actual cubic lattice does not intertwine its derived Clifford corner symbol")
        rows.append({"reduced_momentum_pi_over2_units": [-j for j in bits], "operator_error_diagnostic": error})
    return {"periodic_sites": 64, "corner_blocks": rows,
            "zero_momentum_poles": [str(native_parameters(j)["mass"]) for j in range(3)],
            "one_site_staggered_components": 1, "momentum_corner_components": 8}


def _cz_word():
    # Native J(-pi)=-CZ. Preserve the scalar phase explicitly.
    return [compiler.pulse("I", math.pi, [0]), compiler.pulse("J", -math.pi, [0, 1])]


def bond_factors(left, right, eta, duration, *, spacing=F(1)):
    """Actual source macros for a parity-dressed oriented staggered bond."""
    if type(left) is not int or type(right) is not int or not 0 <= left < right < 512 or type(eta) is not int or eta not in (-1, 1):
        raise ValueError("ordered CAR endpoints and a signed staggered phase required")
    duration, spacing = _rational(duration), _rational(spacing, positive=True)
    phi = -eta*math.pi/2
    result = []
    for middle in range(left+1, right):
        result.append((_cz_word(), (middle, right)))
    result.append((_phase_word(-phi), (right,)))
    result.append((echo_word(duration/spacing), (left, right)))
    result.append((_phase_word(phi), (right,)))
    for middle in range(left+1, right):
        result.append((_cz_word(), (middle, right)))
    return tuple(result)


def _fock_hamiltonian(one_body):
    """Independent bitwise CAR lift, retaining all occupied sectors."""
    h = np.asarray(one_body, complex)
    modes = len(h)
    if modes > 8:
        raise ValueError("explicit full-Fock control is limited to eight modes")
    matrix = np.zeros((1 << modes,)*2, complex)
    for before in range(1 << modes):
        for right in range(modes):
            if not before >> right & 1:
                continue
            after_remove = before ^ (1 << right)
            first_sign = (-1)**((before & ((1 << right)-1)).bit_count())
            for left in range(modes):
                if h[left, right] == 0 or after_remove >> left & 1:
                    continue
                sign = first_sign*(-1)**((after_remove & ((1 << left)-1)).bit_count())
                matrix[after_remove | (1 << left), before] += sign*h[left, right]
    return matrix


def decoded_occupation_word(word, modes=8):
    """Literal native decoded-wire map; every other code bit stays identity."""
    literal = compiler.remap(word, tuple(9*j for j in range(modes)))
    checked = compiler._validate(literal, 9*modes)
    if any(target % 9 for event in checked for target in event["targets"]):
        raise ArithmeticError("source request touches a code syndrome or spectator wire")
    lowered = [compiler.pulse(event["kind"], event["angle"], [target//9 for target in event["targets"]]) for event in checked]
    if lowered != list(word):
        raise ArithmeticError("source factor order changed under native occupation wiring")
    return literal


def _ideal_local_word(word):
    """Independent intended Pauli/relative-phase action for at most two bits."""
    checked = compiler._validate(word)
    bits = max(target for event in checked for target in event["targets"])+1
    if bits > 2:
        raise ValueError("local source macro target is limited to two bits")
    dim = 1 << bits
    result = np.eye(dim, dtype=complex)
    for event in checked:
        kind, targets, angle = event["kind"], event["targets"], event["angle"]
        generator = np.zeros((dim, dim), complex)
        for mask in range(dim):
            bit = (mask >> targets[0]) & 1
            if kind == "J":
                generator[mask, mask] = int(not all(mask >> target & 1 for target in targets))
            elif kind == "I":
                generator[mask, mask] = 1
            elif kind == "Z":
                generator[mask, mask] = 1-2*bit
            elif kind in ("X", "Y"):
                generator[mask ^ (1 << targets[0]), mask] = 1 if kind == "X" else 1j*(1-2*bit)
            elif kind == "H":
                generator[mask, mask] = (1-2*bit)/math.sqrt(2)
                generator[mask ^ (1 << targets[0]), mask] = 1/math.sqrt(2)
            else:
                raise ValueError("unrecognized intended primitive")
        result = expm(1j*angle*generator)@result
    return result


def _schedule(shape, spacing, sector, time, steps):
    coordinates = _coordinates(shape)
    index = {x: j for j, x in enumerate(coordinates)}
    pieces = []
    dt = time/(2*steps)
    m = native_parameters(sector)["mass"]
    for j, x in enumerate(coordinates):
        # Signed native onsite baselines cancel between the balanced even
        # and odd sublattices. No unknown global phase is subtracted.
        pieces.append(((_phase_word(-float(dt*m*_epsilon(x))), (j,)),))
    for axis in range(3):
        for j, x in enumerate(coordinates):
            y = list(x); y[axis] += 1
            if y[axis] >= shape[axis]:
                continue
            pieces.append(bond_factors(j, index[tuple(y)], _eta(x, axis), dt, spacing=spacing))
    return tuple(pieces)


def source_cube(*, sector=0, spacing=F(1), time=TIME, steps=STEPS):
    """All256 Fock columns of an open2³ source cube, with full parity words."""
    spacing, time = _rational(spacing, positive=True), _rational(time)
    if type(steps) is not int or not 1 <= steps <= 1024 or time > 16:
        raise ValueError("one through1024 source split steps and duration at most16 required")
    pieces = _schedule((2, 2, 2), spacing, sector, time, steps)
    catalogue, order, resources = [], [], {key: 0 for key in gaussian._ADDITIVE_RESOURCES}
    value = np.eye(256, dtype=complex)
    cache = {}
    for _ in range(steps):
        for piece in (*pieces, *reversed(pieces)):
            for local_word, sites in piece:
                key = tuple((event["kind"], event["angle"], tuple(event["targets"])) for event in local_word)
                if key not in cache:
                    matrix = _freeze(compiler.execute(local_word, max(max(event["targets"]) for event in local_word)+1))
                    if np.linalg.norm(matrix-_ideal_local_word(local_word), 2) > MACRO_GATE:
                        raise ArithmeticError("local source macro exceeded its all-input phase/operator admission")
                    cache[key] = len(catalogue)
                    catalogue.append({"matrix": matrix, "word": tuple(local_word), "resources": compiler.resources(local_word)})
                block = cache[key]
                order.append((block, sites))
                # Literal decoded occupation requests are validated before
                # lowering. Other syndrome and spectator bits remain identity.
                decoded_occupation_word(compiler.remap(local_word, sites))
                value = attachment._little_apply(catalogue[block]["matrix"], sites, 8, value)
                for name in resources:
                    resources[name] += catalogue[block]["resources"][name]
    h = np.asarray(lattice_operator((2, 2, 2), spacing=spacing, sector=sector, periodic=False), complex)
    target = expm(-1j*float(time)*_fock_hamiltonian(h))
    error = float(np.linalg.norm(value-target, 2))
    q = (native_parameters(sector)["mass"],)+(native_parameters()["exchange"]/spacing,)*3
    # Four commuting subgroups (mass and each spatial direction). On the
    # open cube their one-body generators anticommute, square to q_i² I8,
    # and dGamma of a traceless Clifford generator has norm4*|q_i|.
    # Prop16 then gives (4/3)*q_i*sum_tail(q²)+(2/3)*q_i²*sqrt(sum_tail(q²)).
    # Replacing that root by sum_tail(q) is directed and purely rational.
    constant = sum((F(4, 3)*v*sum(x*x for x in q[j+1:])+F(2, 3)*v*v*sum(q[j+1:])
                    for j, v in enumerate(q)), F())
    splitting = F(time)**3*constant/steps**2
    summed_macro_error = F(len(order), 10**10)
    macro_bound = summed_macro_error/(1-summed_macro_error)
    if error > float(splitting+macro_bound):
        raise ArithmeticError("full-Fock source cube exceeded its finite splitting/source budget")
    vacuum_error = abs(value[0, 0]-1)
    if vacuum_error > len(order)*MACRO_GATE:
        raise ArithmeticError("source cube lost its absolute empty-state phase")
    return {"unitary": _freeze(value), "target": _freeze(target), "one_body_hamiltonian": _freeze(h),
            "catalogue": catalogue, "factor_order": order, "resources": resources,
            "full_Fock_columns": 256, "native_CAR_modes": 8, "steps": steps,
            "source_coordinate_duration": str(time), "lattice_spacing": str(spacing),
            "splitting_operator_upper_exact": str(splitting),
            "composed_macro_operator_upper_conditional": str(macro_bound),
            "source_local_macro_operator_admission": "1/10000000000",
            "target_operator_error_diagnostic": error,
            "empty_Fock_scalar_phase_error_diagnostic": float(vacuum_error),
            "complete_source_output_normalized_or_polar_corrected": False,
            "full_native_many_carrier_service_history_executed": False}


def fiber_transport_certificate(momentum=(F(1, 3), F(1, 5), F(1, 7)), *, spacing=F(1, 4)):
    """Explicit identification used ONLY for normal-ordered excitation bounds."""
    momentum, spacing = tuple(map(_rational_signed, momentum)), _rational(spacing, positive=True)
    if len(momentum) != 3:
        raise ValueError("three momentum components required")
    matrices = tuple(np.asarray(g, complex) for g in clifford_matrices())
    m, c = map(float, (native_parameters()["mass"], native_parameters()["coordinate_speed"]))
    beta = matrices[3]
    kinetic = [sum((c*(math.sin(float(spacing*k))/float(spacing) if lattice else float(k))*g
                   for k, g in zip(momentum, matrices[:3], strict=True)), np.zeros((8, 8), complex))
               for lattice in (True, False)]
    rotations, generators, energies = [], [], []
    for K in kinetic:
        energy = math.sqrt(m*m+float(np.trace(K@K).real)/8)
        R = ((energy+m)*np.eye(8)+K@beta)/math.sqrt(2*energy*(energy+m))
        h = K+m*beta
        if max(np.linalg.norm(R.conj().T@R-np.eye(8), 2), np.linalg.norm(R.conj().T@h@R-energy*beta, 2)) > SOURCE_GATE:
            raise ArithmeticError("declared spectral-fiber identification failed")
        rotations.append(R); generators.append(h); energies.append(energy)
    transported = rotations[0].conj().T@generators[0]@rotations[0]-rotations[1].conj().T@generators[1]@rotations[1]
    difference = abs(energies[0]-energies[1])
    if np.linalg.norm(transported-(energies[0]-energies[1])*beta, 2) > SOURCE_GATE:
        raise ArithmeticError("normal-ordered generator comparison hid an extra conjugation error")
    return {"matching": "R_a†h_aR_a=E_a beta and R_0†h_0R_0=E_0 beta; identify those spectral coordinates",
            "transported_generator_norm_diagnostic": float(np.linalg.norm(transported, 2)),
            "energy_difference_diagnostic": difference,
            "unchanged_local_detectors_covered_by_excitation_bound": False,
            "local_CAR_comparison_uses_separate_projector_and_two_point_bounds": True}


def hierarchy_certificate():
    gaps = [native_parameters(j)["mass"] for j in range(3)]
    if tuple(gap/min(gaps) for gap in reversed(gaps)) != (1, 3, 12):
        raise ArithmeticError("native free-field pole ratio changed")
    return {"zero_momentum_poles": list(map(str, gaps)), "ascending_mass_ratio": [1, 3, 12],
            "free_self_energy": "0", "kinetic_wavefunction_coefficient": "1",
            "constant_positive_field_rescaling_changes_pole_zeros": False,
            "reason": "det[Z_s*(omega-h_s)]=Z_s**8*det(omega-h_s); canonical normalization removes Z_s",
            "common_spatial_scale_or_speed_changes_zero_momentum_mass": False,
            "common_clock_or_lapse_changes_mass_ratios": False,
            "species_dependent_lapse_or_mass_renormalization_derived": False,
            "un_echoed_interacting_exchange_is_this_Hamiltonian": False,
            "measured_mass_inputs_or_target_fitting_used": False}


def open_cube_component_certificate():
    h = lattice_operator((2, 2, 2), periodic=False)
    coordinates = _coordinates((2, 2, 2))
    components = [s.diag(*h.diagonal())]
    for axis in range(3):
        matrix = s.zeros(8)
        for i, x in enumerate(coordinates):
            for j, y in enumerate(coordinates):
                if x[axis] != y[axis] and all(x[k] == y[k] for k in range(3) if k != axis):
                    matrix[i, j] = h[i, j]
        components.append(matrix)
    q = (s.Rational(native_parameters()["mass"]),)+(s.Rational(native_parameters()["exchange"]),)*3
    if sum(components, s.zeros(8)) != h:
        raise ArithmeticError("source direction grouping changed the native lattice operator")
    for i, left in enumerate(components):
        for j, right in enumerate(components):
            if left*right+right*left != (2*q[i]**2*s.eye(8) if i == j else s.zeros(8)):
                raise ArithmeticError("source cube direction groups do not have the declared commutator bound")
    return {"anticommuting_one_body_groups": 4, "group_amplitudes": list(map(str, q)),
            "full_Fock_group_norms": [str(4*v) for v in q],
            "Strang_group_order": ["onsite", "x", "y", "z", "z", "y", "x", "onsite"],
            "commuting_disjoint_bonds_grouped_before_splitting": True}


def _exterior(one_body):
    modes = len(one_body)
    if modes > 8:
        raise ValueError("explicit exterior target limited to eight CAR modes")
    occupied = [[j for j in range(modes) if mask >> j & 1] for mask in range(1 << modes)]
    target = np.zeros((1 << modes,)*2, complex)
    for before, columns in enumerate(occupied):
        for after, rows in enumerate(occupied):
            if len(rows) == len(columns):
                target[after, before] = np.linalg.det(one_body[np.ix_(rows, columns)]) if rows else 1.
    return target


def _vacuum_frame():
    h = np.asarray(lattice_operator((2, 2, 2), periodic=False), complex)
    beta = np.diag([_epsilon(x) for x in _coordinates((2, 2, 2))])
    m, J = native_parameters()["mass"], native_parameters()["exchange"]
    energy = math.sqrt(float(m*m+3*J*J))
    kinetic = h-float(m)*beta
    frame = ((energy+float(m))*np.eye(8)+kinetic@beta)/math.sqrt(2*energy*(energy+float(m)))
    negative = [j for j in range(8) if beta[j, j] == -1]
    target = (np.eye(8)-h/energy)/2
    if max(np.linalg.norm(frame.conj().T@frame-np.eye(8), 2),
           np.linalg.norm(frame[:, negative]@frame[:, negative].conj().T-target, 2)) > SOURCE_GATE:
        raise ArithmeticError("analytic half-filled native ground frame failed")
    return frame, target, negative, -4*energy


def source_vacuum():
    """Prepare the negative-energy sea by actual number-conserving source words.

    The declared seed has four encoded charges. It is not an autonomous
    cooling result or a claim that C/E repair selects this many-body vacuum.
    All256 columns of its Gaussian frame and all output amplitudes remain.
    """
    frame, projector, negative, ground_energy = _vacuum_frame()
    work, calls = frame.copy(), []
    for column in range(8):
        for row in range(7, column, -1):
            a, b = work[row-1, column], work[row, column]
            norm = math.hypot(abs(a), abs(b))
            cosine = abs(a)/norm if norm else 1.
            sine = (a/abs(a) if a else 1.)*b.conjugate()/norm if norm else 0.
            rotation = np.array([[cosine, sine], [-sine.conjugate(), cosine]])
            work[[row-1, row]] = rotation@work[[row-1, row]]
            coefficient = 1j*sine/abs(sine) if sine else 1j
            calls.append((row-1, row, coefficient, math.atan2(abs(b), abs(a))))
    if np.linalg.norm(work-np.diag(np.diag(work)), 2) > SOURCE_GATE:
        raise ArithmeticError("native vacuum frame source elimination failed")
    value = np.eye(256, dtype=complex)
    resources = {key: 0 for key in gaussian._ADDITIVE_RESOURCES}
    for mode, phase in enumerate(np.angle(np.diag(work))):
        word = _phase_word(float(phase))
        actual = compiler.execute(word, 1)
        value = attachment._little_apply(actual, (mode,), 8, value)
        for key in resources:
            resources[key] += compiler.resources(word)[key]
    for left, right, coefficient, angle in reversed(calls):
        gate = gaussian.local_gate("hopping", coefficient, -angle)
        value = attachment._little_apply(gate.fock_unitary, (left, right), 8, value)
        for key in resources:
            resources[key] += gate.resources[key]
    target = _exterior(frame)
    error = float(np.linalg.norm(value-target, 2))
    if error > SOURCE_GATE:
        raise ArithmeticError("actual all-Fock source vacuum frame exceeds its phase/operator gate")
    seed = sum(1 << j for j in negative)
    state = _freeze(value[:, seed])
    probability = float(np.vdot(state, state).real)
    wrong_charge = math.sqrt(sum(abs(state[mask])**2 for mask in range(256) if mask.bit_count() != 4))
    H = _fock_hamiltonian(np.asarray(lattice_operator((2, 2, 2), periodic=False), complex))
    residual = float(np.linalg.norm(H@state-ground_energy*state))
    if max(abs(probability-1), wrong_charge, residual) > SOURCE_GATE:
        raise ArithmeticError("source prepared state is not the declared half-filled ground state")
    # Costs retain a fixed angle envelope instead of arbitrary QR phases
    # of numerical empty pivots. Actual output amplitudes are still checked.
    costs = {key: count for key, count in resources.items() if key not in ("zero_angle_requests", "integrated_processor_norm_upper")}
    return {"frame_unitary": _freeze(value), "state": state, "target_projector": _freeze(projector),
            "analytic_orbital_frame": _freeze(frame), "resources": costs,
            "source_phase_angle_envelope": "pi per primitive request",
            "source_ground_energy": ground_energy, "encoded_charge": 4,
            "empty_encoded_state_is_ground": False,
            "native_CE_selects_or_cools_this_state": False,
            "frame_operator_error_diagnostic": error,
            "ground_eigenvector_error_diagnostic": residual}


def history_instrument(unitary, state):
    """A continuing3D density read/feedback/read instrument, all input columns."""
    first = instruments.pauli_read(unitary, 0, axis="Z")
    matrices = []
    pointer_cost = first.resources["native_pulses"]
    for record, branch in enumerate(first.outputs):
        if record:
            feedback = compiler.execute(_phase_word(math.pi), 1)
            branch = attachment._little_apply(feedback, (7,), 8, branch)
        second = instruments.pauli_read(unitary@branch, 7, axis="Z")
        pointer_cost += second.resources["native_pulses"]
        matrices.extend(second.outputs)
    outputs = np.stack(matrices)
    probabilities = [float(np.linalg.norm(K@state)**2) for K in outputs]
    defect = float(np.linalg.norm(sum((K.conj().T@K for K in outputs), np.zeros((256, 256), complex))-np.eye(256), 2))
    if defect > SOURCE_GATE or abs(sum(probabilities)-1) > SOURCE_GATE:
        raise ArithmeticError("complete source history instrument lost its Gram identity")
    return {"outputs": _freeze(outputs), "record_probabilities": probabilities,
            "all_unknown_Fock_input_columns": 256, "density_record_addresses": [[0, 0], [0, 1], [1, 0], [1, 1]],
            "first_read_backaction_and_record_feedback_retained": True,
            "shared_observers_are_marginals_of_one_retained_instrument": True,
            "arbitrary_untouched_quantum_reference_supported": True,
            "physical_clock_or_native_service_flags_authenticated": False,
            "pointer_native_pulses_per_complete_schedule": pointer_cost,
            "trace_preservation_error_diagnostic": defect}


def execute():
    echo = echo_certificate()
    cube, vacuum = source_cube(), source_vacuum()
    history = history_instrument(cube["unitary"], vacuum["state"])
    arrays = {"native_exchange_echo": echo["unitary"], "source_cube_unitary": cube["unitary"],
              "source_cube_one_body_H": cube["one_body_hamiltonian"],
              "source_vacuum_frame": vacuum["frame_unitary"], "source_filled_vacuum": vacuum["state"],
              "source_history_kraus": history["outputs"]}
    catalogue = []
    for j, block in enumerate(cube["catalogue"]):
        name = "source_macro_"+str(j)
        arrays[name] = block["matrix"]
        catalogue.append({"array": name, "word": list(block["word"]), "resources": block["resources"]})
    report = {"schema": "oph.native_staggered_field.v1",
        "native_parameters": {"masses": [str(native_parameters(j)["mass"]) for j in range(3)],
            "J": str(native_parameters()["exchange"]), "coordinate_speed": str(native_parameters()["coordinate_speed"]),
            "native_empty_occupied_code_labels": [[0, 1], [3, 2], [5, 4]],
            "other_sector_chart_rule": "field occupation is after the declared fixed hole-basis recoding in sectors1/2; not the unrelabelled code parity bit"},
        "echo": {"generator_identity": "(B-ZBZ)/2=-J*(a†b+b†a), [B,ZBZ]=0",
            "full_Fock_columns": 4, "double_occupation_and_scalar_phase_retained": True,
            "word": list(echo["word"]), "resources": echo["resources"],
            "positive_B_and_its_inverse_and_matching_Z_are_chosen_controls": True},
        "periodic_operator_and_corner_check": periodic_fourier_certificate(),
        "spin_and_taste": clifford_certificate(), "finite_cube_commutator_certificate": open_cube_component_certificate(),
        "source_cube": {"shape": [2, 2, 2], "boundary": "open; not a periodic2-site derivative",
            "literal_executed_native_sector": 0,
            "all_Fock_columns": 256, "CAR_modes": 8, "time": str(TIME), "spacing": "1", "steps": STEPS,
            "catalogue": catalogue,
            "factor_order": [[block, list(sites)] for block, sites in cube["factor_order"]],
            "resources": cube["resources"], "splitting_operator_upper_exact": cube["splitting_operator_upper_exact"],
            "composed_macro_operator_upper_conditional": cube["composed_macro_operator_upper_conditional"],
            "full_native_carrier_CE_or_physical_service_history_executed": False,
            "parity_strings_and_syndrome_free_native_occupation_wiring_checked": True},
        "vacuum": {"encoded_charge": 4, "meaning": "all four negative-energy orbitals filled",
            "relative_particle_hole_charge": "Q=N-4 on the eight-site source control",
            "empty_encoded_state_is_ground": False, "chosen_seed_and_coherent_frame_preparation": True,
            "native_repair_selects_this_vacuum": False, "resources": vacuum["resources"],
            "raw_empty_state_or_filled_vacuum_global_phase_removed": False},
        "continuing_observer": {key: value for key, value in history.items()
            if key not in ("outputs", "trace_preservation_error_diagnostic")},
        "continuum_band_contracts": [continuum_contract(spacing=F(1, 2**j), excitations=4) for j in range(2, 9)],
        "spectral_fiber_identification": fiber_transport_certificate(), "hierarchy": hierarchy_certificate(),
        "source_contract": {"local_complete_macro_operator_admission": "1/10000000000",
            "full_preparation_pointer_and_Gram_operator_admission": "1/100000000",
            "complete_array_operator_replay_tolerance": "1/500000000",
            "all_raw_complex_outputs_retained": True,
            "norm_projection_polar_correction_or_success_postselection": False,
            "numerical_admissions_are_outward_or_physical_certificates": False,
            "product_formula_bound_primary_source": "https://arxiv.org/abs/1912.08854; Proposition16",
            "full_growing_field_error_needs_additional_source_precision_and_service_scaling": True},
        "scope": {"supplied_logical_Dirac_matrix_or_empirical_mass_packet": False,
            "chosen_cubic_geometry_staggered_signs_and_inverse_echo_controls": True,
            "continuum_is_free_vectorlike_two_taste_Dirac_per_native_sector": True,
            "source_open_cube_is_a_continuum_or_a_periodic64_mode_history": False,
            "physical_electric_weak_color_charge_or_Higgs_attachment_derived": False,
            "SM_chiral_interacting_continuum_or_mass_hierarchy_derived": False,
            "vacuum_energy_renormalization_or_gravitational_counterterm_derived": False,
            "programmed_coordinate_is_a_calibrated_physical_or_native_attempt_clock": False},
        "primary_staggered_sources": ["https://doi.org/10.1103/PhysRevD.16.3031",
            "https://arxiv.org/html/2604.02078v1; sections2.1--2.2"]}
    report["backend_diagnostics"] = verify_retained_arrays(report, arrays)
    return report, arrays


def _stack_norm(value):
    value = np.asarray(value)
    return float(np.linalg.norm(value.reshape(-1, value.shape[-1]), 2)) if value.ndim > 2 else float(np.linalg.norm(value, 2))


def verify_retained_arrays(report, arrays):
    """Gate actual retained phase/operator arrays before portable replay."""
    if any(type(value) is not np.ndarray or value.dtype != np.dtype("complex128") or not np.isfinite(value).all()
           for value in arrays.values()):
        raise ValueError("complete finite complex128 source arrays required")
    shapes = {"native_exchange_echo": (4, 4), "source_cube_unitary": (256, 256),
        "source_cube_one_body_H": (8, 8), "source_vacuum_frame": (256, 256),
        "source_filled_vacuum": (256,), "source_history_kraus": (4, 256, 256)}
    for name, shape in shapes.items():
        if name not in arrays or arrays[name].shape != shape:
            raise ValueError("complete source array schema differs: "+name)
    echo_target = expm(-1j*float(F(3, 5))*np.array([[0, 0, 0, 0], [0, 0, -float(native_parameters()["exchange"]), 0],
        [0, -float(native_parameters()["exchange"]), 0, 0], [0, 0, 0, 0]], complex))
    echo_error = _stack_norm(arrays["native_exchange_echo"]-echo_target)
    if echo_error > MACRO_GATE:
        raise ValueError("retained exchange echo exceeds its all-Fock phase gate")
    value = np.eye(256, dtype=complex)
    maximum_macro = 0.
    for block in report["source_cube"]["catalogue"]:
        matrix = arrays[block["array"]]
        target = _ideal_local_word(block["word"])
        if matrix.shape != target.shape:
            raise ValueError("local source macro schema differs")
        error = _stack_norm(matrix-target)
        maximum_macro = max(maximum_macro, error)
        if error > MACRO_GATE:
            raise ValueError("retained local source macro exceeds its operator gate")
    for block, sites in report["source_cube"]["factor_order"]:
        item = report["source_cube"]["catalogue"][block]
        decoded_occupation_word(compiler.remap(item["word"], sites))
        value = attachment._little_apply(arrays[item["array"]], tuple(sites), 8, value)
    U = arrays["source_cube_unitary"]
    if _stack_norm(value-U) > REPLAY_TOLERANCE:
        raise ValueError("retained source cube does not equal its complete factor word")
    h = np.asarray(lattice_operator((2, 2, 2), periodic=False), complex)
    if not np.array_equal(h, arrays["source_cube_one_body_H"]):
        raise ValueError("retained Hamiltonian does not equal the native-derived local law")
    target = expm(-1j*float(TIME)*_fock_hamiltonian(h))
    target_error = _stack_norm(U-target)
    admitted = F(report["source_cube"]["splitting_operator_upper_exact"])+F(report["source_cube"]["composed_macro_operator_upper_conditional"])
    if target_error > float(admitted)+REPLAY_TOLERANCE:
        raise ValueError("retained full-Fock source exceeds its splitting/operator budget")
    number = np.diag([mask.bit_count() for mask in range(256)])
    if max(_stack_norm(U.conj().T@U-np.eye(256)), _stack_norm(U@number-number@U)) > SOURCE_GATE:
        raise ValueError("retained source is not unitary and number preserving within its admission")
    frame, _, negative, energy = _vacuum_frame()
    prepared = arrays["source_vacuum_frame"]
    preparation_error = _stack_norm(prepared-_exterior(frame))
    if preparation_error > SOURCE_GATE:
        raise ValueError("retained vacuum source frame exceeds its all-Fock phase gate")
    state = arrays["source_filled_vacuum"]
    if not np.array_equal(state, prepared[:, sum(1 << j for j in negative)]):
        raise ValueError("retained vacuum state lost its exact in-archive source ancestry")
    eigen_error = float(np.linalg.norm(_fock_hamiltonian(h)@state-energy*state))
    if eigen_error > SOURCE_GATE:
        raise ValueError("retained filled vacuum is not the source-derived ground state")
    kraus = arrays["source_history_kraus"]
    expectations = []
    P0 = np.diag([1-int(mask & 1) for mask in range(256)])
    P7 = np.diag([1-int((mask >> 7) & 1) for mask in range(256)])
    feedback = np.diag([(-1)**((mask >> 7) & 1) for mask in range(256)])
    for first in (0, 1):
        for second in (0, 1):
            expectations.append((P7 if not second else np.eye(256)-P7)@U@(feedback if first else np.eye(256))@(P0 if not first else np.eye(256)-P0)@U)
    pointer_error = _stack_norm(kraus-np.stack(expectations))
    gram_error = _stack_norm(sum((K.conj().T@K for K in kraus), np.zeros((256, 256), complex))-np.eye(256))
    if max(pointer_error, gram_error) > SOURCE_GATE:
        raise ValueError("retained continuing instrument exceeds its complete operator/Gram gate")
    probabilities = [float(np.linalg.norm(K@state)**2) for K in kraus]
    if any(abs(a-b) > REPLAY_TOLERANCE for a, b in zip(probabilities, report["continuing_observer"]["record_probabilities"], strict=True)):
        raise ValueError("retained observer probabilities differ from complete source outputs")
    return {"exchange_echo_operator_residual": echo_error, "maximum_macro_operator_residual": maximum_macro,
        "cube_Fock_target_residual": target_error, "vacuum_preparation_operator_residual": preparation_error,
        "vacuum_ground_residual": eigen_error, "history_pointer_operator_residual": pointer_error,
        "history_Gram_residual": gram_error}


def source_bindings():
    paths = set(attachment.source_bindings())|{
        "model/common_source/native_staggered_field.py", "model/common_source/test_native_staggered_field.py"}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def _same(expected, actual, path="report"):
    if type(expected) is not type(actual):
        raise ValueError("original JSON type differs: "+path)
    if type(expected) is dict:
        if expected.keys() != actual.keys():
            raise ValueError("complete report keys differ: "+path)
        for key in expected:
            _same(expected[key], actual[key], path+"."+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError("complete report list differs: "+path)
        for j, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _same(left, right, path+"["+str(j)+"]")
    elif type(expected) is float:
        if not math.isfinite(actual) or (expected != actual if path.endswith(".angle") else not math.isclose(expected, actual, abs_tol=REPLAY_TOLERANCE, rel_tol=2e-10)):
            raise ValueError("portable numerical observable differs: "+path)
    elif expected != actual:
        raise ValueError("exact source/architecture value differs: "+path)


def compare_execution(report, arrays, fresh, fresh_arrays):
    left, right = deepcopy(report), deepcopy(fresh)
    left_metrics, right_metrics = left.pop("backend_diagnostics"), right.pop("backend_diagnostics")
    _same(right, left)
    if arrays.keys() != fresh_arrays.keys():
        raise ValueError("complete native staggered source array keys differ")
    for name, expected in fresh_arrays.items():
        actual = arrays[name]
        if type(actual) is not np.ndarray or actual.dtype != expected.dtype or actual.shape != expected.shape or not np.isfinite(actual).all():
            raise ValueError("source array schema differs: "+name)
        if _stack_norm(actual-expected) > REPLAY_TOLERANCE:
            raise ValueError("complete source output operator differs: "+name)
    for metadata, original, values in ((left_metrics, report, arrays), (right_metrics, fresh, fresh_arrays)):
        metrics = verify_retained_arrays(original, values)
        if metadata.keys() != metrics.keys() or any(type(metadata[key]) is not float or not math.isfinite(metadata[key])
                or abs(metadata[key]-metrics[key]) > 2e-12 for key in metrics):
            raise ValueError("reported residual is not a diagnostic of its own original source arrays")


def check(path):
    path = Path(path)
    report = strict_json.load(path)
    attachment._finite_json(report)
    pins = source_bindings()
    if report.get("source_bindings") != pins:
        raise ValueError("native staggered closed source bindings differ")
    artifact = path.with_suffix(".npz")
    if report.get("artifact") != {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError("native staggered complete artifact custody differs")
    with np.load(artifact, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError("duplicate complete source array names")
        arrays = {name: archive[name] for name in archive.files}
    fresh, fresh_arrays = execute()
    scientific = {key: value for key, value in report.items() if key not in ("source_bindings", "artifact")}
    compare_execution(scientific, arrays, fresh, fresh_arrays)
    if source_bindings() != pins:
        raise ValueError("native staggered sources changed during replay")
    return scientific


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--output", type=Path)
    choice.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print("native exchange echo, staggered field, filled vacuum and complete observer replay: PASS")
        return
    artifact = args.output.with_suffix(".npz")
    if args.output.exists() or artifact.exists():
        parser.error("fresh JSON/NPZ output paths required")
    pins = source_bindings()
    report, arrays = execute()
    if source_bindings() != pins:
        raise RuntimeError("native staggered sources changed during execution")
    report["source_bindings"] = pins
    with artifact.open("xb") as stream:
        np.savez_compressed(stream, **arrays)
    report["artifact"] = {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False)+"\n")
    print(json.dumps(report["continuing_observer"]["record_probabilities"]))


if __name__ == "__main__":
    main()
