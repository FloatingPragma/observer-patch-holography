"""Finite continuum Weyl symbols for the declared M15 x F3 Higgs candidate.

These matrices attach a canonical spatial Weyl operator to the exact internal
action. They are a momentum-cutoff diagnostic, not a local chiral gauge lattice
regulator. The source adapter executes the corresponding CAR Hamiltonian
separately; diagonalizing this symbol is not itself a native source execution.
"""
from __future__ import annotations

import itertools
import math
from functools import lru_cache
from numbers import Real

import numpy as np
import sympy as s
from sympy.polys.matrices import DomainMatrix

from . import action
from ...observer_particles.spectral_response import charge_resolved_spectrum, spectrum_summary


def _numeric(matrix):
    return np.asarray(matrix.evalf(17) if hasattr(matrix, "evalf") else matrix,
                      dtype=complex)


def _momentum(momentum):
    try:
        values = tuple(itertools.islice(iter(momentum), 4))
    except TypeError as exc:
        raise ValueError("three finite real momentum components required") from exc
    if len(values) != 3 or any(isinstance(x, (bool, np.bool_))
                              or not isinstance(x, Real) for x in values):
        raise ValueError("three finite real momentum components required")
    try:
        result = tuple(float(x) for x in values)
    except (OverflowError, ValueError) as exc:
        raise ValueError("finite execution momenta required") from exc
    if not all(math.isfinite(x) for x in result):
        raise ValueError("three finite real momentum components required")
    return result


def pauli_matrices():
    return (np.array([[0, 1], [1, 0]], complex),
            np.array([[0, -1j], [1j, 0]], complex),
            np.diag([1., -1.]).astype(complex))


def spin_epsilon():
    return np.array([[0, 1], [-1, 0]], complex)


def coordinate_speed():
    """c=Omega/2 in the declared spatial coordinate/modular-parameter units."""
    return float(action.source_coefficients()["Omega"] / 2)


def _speed(value):
    value = coordinate_speed() if value is None else value
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError("a finite positive coordinate speed required")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("a finite positive coordinate speed required")
    return value


def _dimension(value):
    if type(value) is not int or value < 1:
        raise ValueError("a positive internal dimension required")
    return value


def _kernel(kernel):
    result = _numeric(action.mass_kernel() if kernel is None else kernel)
    if result.ndim != 2 or result.shape[0] != result.shape[1] or not len(result):
        raise ValueError("a nonempty symmetric internal mass kernel required")
    if not np.isfinite(result).all() or not np.allclose(result, result.T, atol=1e-13, rtol=1e-13):
        raise ValueError("a finite symmetric internal mass kernel required")
    return result


def weyl_kinetic(momentum, internal_dimension=45, *, speed=None):
    """All-left h(p)=-c sigma.p tensor I, with spin outermost."""
    p, n, c = _momentum(momentum), _dimension(internal_dimension), _speed(speed)
    sigma_p = sum((x * sigma for x, sigma in zip(p, pauli_matrices(), strict=True)),
                  np.zeros((2, 2), complex))
    return -c * np.kron(sigma_p, np.eye(n))


def hamiltonian_pairing(kernel=None):
    """Creation Delta=epsilon tensor K* for L_mass=+1/2 chi^T(epsilon K)chi+h.c.

    The conjugation is essential for a complex symmetric K. The default native
    neutral-Higgs kernel is real, so this agrees with action.pairing_kernel().
    """
    return np.kron(spin_epsilon(), _kernel(kernel).conj())


def bdg_symbol(momentum, kernel=None, *, speed=None):
    """Hermitian symbol on (chi(p), chi†(-p)), dimension4*dim(K).

    The lower block is -h(-p)^T. Opposite momentum and Nambu redundancy must
    both be respected when turning this symbol into a many-body CAR source.
    """
    p, k = _momentum(momentum), _kernel(kernel)
    h = weyl_kinetic(p, len(k), speed=speed)
    hminus = weyl_kinetic(tuple(-x for x in p), len(k), speed=speed)
    delta = hamiltonian_pairing(k)
    return np.block([[h, delta], [delta.conj().T, -hminus.T]])


def spin_conjugated_symbol(momentum, kernel=None, *, speed=None):
    """Same symbol after chi†(-p) -> epsilon chi†(-p); displays L/R mixing."""
    k = _kernel(kernel)
    d = 2 * len(k)
    identity, zero = np.eye(d), np.zeros((d, d))
    change = np.block([[identity, zero], [zero, np.kron(spin_epsilon(), np.eye(len(k)))]])
    return change @ bdg_symbol(momentum, k, speed=speed) @ change.conj().T


def momentum_pair_source(momentum, kernel=None, *, speed=None):
    """h,Delta for actual CAR modes (chi(p),chi(-p)) at distinct momenta.

    At p=0 use weyl_kinetic(0) and hamiltonian_pairing directly instead. The
    pair construction doubles the original CAR modes, not the species count.
    """
    p, k = _momentum(momentum), _kernel(kernel)
    if not any(x != 0 for x in p):
        raise ValueError("p=0 is self-conjugate; do not duplicate its CAR modes")
    hp = weyl_kinetic(p, len(k), speed=speed)
    hm = weyl_kinetic(tuple(-x for x in p), len(k), speed=speed)
    delta, zero = hamiltonian_pairing(k), np.zeros_like(hp)
    return (np.block([[hp, zero], [zero, hm]]),
            np.block([[zero, delta], [-delta.T, zero]]))


def nambu_spin(internal_dimension=45):
    """Physical spin generators before conjugating the hole spin basis."""
    n = _dimension(internal_dimension)
    result = []
    for sigma in pauli_matrices():
        s = np.kron(sigma / 2, np.eye(n))
        z = np.zeros_like(s)
        result.append(np.block([[s, z], [z, -s.T]]))
    return tuple(result)


def helicity(momentum, internal_dimension=45):
    p = _momentum(momentum)
    norm = math.hypot(*p)
    if norm == 0:
        raise ValueError("helicity is undefined at zero momentum")
    return sum((x / norm * s for x, s in zip(p, nambu_spin(internal_dimension), strict=True)),
               np.zeros((4 * internal_dimension, 4 * internal_dimension), complex))


def nambu_charge(charge=None):
    q = _numeric(action.electric_charge() if charge is None else charge)
    if q.ndim != 2 or q.shape[0] != q.shape[1] or not np.isfinite(q).all():
        raise ValueError("a finite square Hermitian internal charge required")
    if not np.allclose(q, q.conj().T, atol=1e-13, rtol=1e-13):
        raise ValueError("a Hermitian internal charge required")
    qs = np.kron(np.eye(2), q)
    zero = np.zeros_like(qs)
    return np.block([[qs, zero], [zero, -qs.T]])


def naive_weyl_nodes():
    """All eight exact zeros of -c sigma.sin(k), with Jacobian chirality."""
    nodes = []
    for bits in itertools.product((0, 1), repeat=3):
        nodes.append({"pi_coordinates": list(bits),
                      "jacobian_chirality": -1 if sum(bits) % 2 == 0 else 1})
    return nodes


def naive_lattice_symbol(momentum, internal_dimension=45, *, speed=None):
    """Negative control: this local sin(k) choice has eight Weyl nodes."""
    return weyl_kinetic(tuple(math.sin(x) for x in _momentum(momentum)),
                        internal_dimension, speed=speed)


@lru_cache(maxsize=1)
def _weak_raising_exact():
    """Recover E34 in the actual 12-current source basis, then lift it to M15F3."""
    matrices = action.matrices()
    fundamental = s.zeros(5)
    fundamental[3, 4] = 1
    basis = s.Matrix.hstack(*(s.Matrix(g).reshape(25, 1) for g in matrices["generators5"]))
    # The Hermitian Gram matrix is real over Q(sqrt5). Use that exact field to
    # avoid enormous unreduced expressions from generic symbolic elimination.
    gram = (basis.H * basis).applyfunc(s.expand)
    rhs = (basis.H * fundamental.reshape(25, 1)).applyfunc(s.expand)
    domain = s.QQ.algebraic_field(s.sqrt(5))
    inverse = DomainMatrix.from_Matrix(gram).convert_to(domain).inv()
    real = DomainMatrix.from_Matrix(s.re(rhs)).convert_to(domain)
    imag = DomainMatrix.from_Matrix(s.im(rhs)).convert_to(domain)
    coefficients = inverse.matmul(real).to_Matrix() + s.I*inverse.matmul(imag).to_Matrix()
    if (basis*coefficients-fundamental.reshape(25, 1)).applyfunc(s.expand) != s.zeros(25, 1):
        raise ArithmeticError("the weak raising current did not reconstruct from actual source currents")
    current = sum((coefficient*g for coefficient, g in
                   zip(coefficients, matrices["generators45"], strict=True)), s.zeros(45))
    return tuple(coefficients), s.ImmutableSparseMatrix(current.applyfunc(s.expand))


def weak_raising_current():
    """Internal normalized weak J+ from the actual source current matrices."""
    return _numeric(_weak_raising_exact()[1])


def weak_chirality_report(momentum=(.23, -.31, .17)):
    """Total J†J weight on charged-lepton energy/helicity branches.

    J†J sums the complete final-state space before projecting the incoming
    electron branch. It therefore avoids a spurious zero from selecting only
    one same-momentum neutrino spin/frequency as the final state.
    """
    p = _momentum(momentum)
    momentum_energy = coordinate_speed()*math.hypot(*p)
    if momentum_energy == 0:
        raise ValueError("the weak-helicity discriminator requires nonzero momentum")
    operators = action.identification_operators(families=False)
    charge, casimir = operators["Q"], operators["weak_C2"]
    left = [i for i in range(15) if charge[i, i] == -1 and casimir[i, i] == s.Rational(3, 4)]
    conjugate = [i for i in range(15) if charge[i, i] == 1 and casimir[i, i] == 0]
    neutrino = [i for i in range(15) if charge[i, i] == 0 and casimir[i, i] == s.Rational(3, 4)]
    if len(left) != 1 or len(conjugate) != 1 or len(neutrino) != 1:
        raise ArithmeticError("the actual charge/weak-Casimir projectors did not isolate L/e^c")
    current = np.kron(np.eye(2), weak_raising_current())
    zero = np.zeros_like(current)
    nambu_current = np.block([[current, zero], [zero, -current.T]])
    total_strength = nambu_current.conj().T @ nambu_current
    hmass = bdg_symbol(p)
    hfree = bdg_symbol(p, np.zeros((45, 45)))
    helicity_matrix = helicity(p)
    kernel = _kernel(None)
    rows = []
    massless = []
    left_trace, hole_trace = 0., 0.
    for family in range(3):
        il, ic = 3*left[0]+family, 3*conjugate[0]+family
        indices = [spin*45+il for spin in range(2)] + [90+spin*45+ic for spin in range(2)]
        hm = hmass[np.ix_(indices, indices)]
        hf = hfree[np.ix_(indices, indices)]
        lam = helicity_matrix[np.ix_(indices, indices)]
        strength = total_strength[np.ix_(indices, indices)]
        left_trace += float(np.trace(strength[:2, :2]).real)
        hole_trace += float(np.trace(strength[2:, 2:]).real)
        mass = float(abs(kernel[il, ic]))
        energy = math.hypot(momentum_energy, mass)
        weights, free_weights, expected = {}, {}, {}
        for sign, name in ((-1, "negative_helicity"), (1, "positive_helicity")):
            helicity_projector = (np.eye(4)+sign*2*lam)/2
            positive_mass = (np.eye(4)+hm/energy)/2
            positive_free = (np.eye(4)+hf/momentum_energy)/2
            branch = positive_mass@helicity_projector
            free_branch = positive_free@helicity_projector
            weights[name] = float(np.trace(branch@strength).real)
            free_weights[name] = float(np.trace(free_branch@strength).real)
            expected[name] = (1-sign*momentum_energy/energy)/2
            if not np.isclose(np.trace(branch), 1, atol=1e-12):
                raise ArithmeticError("energy/helicity projection did not isolate one incoming mode")
        rows.append({"native_family_axis": family, "mass_parameter": mass,
                     "positive_energy": energy, "total_weak_current_weight": weights,
                     "predicted_chiral_mixing_weight": expected})
        massless.append({"native_family_axis": family, "total_weak_current_weight": free_weights})
    return {
        "current": "normalized E34 reconstructed from the actual twelve source generators, then lifted to the retained exterior representation",
        "source_current_coefficients": [str(s.expand(x)) for x in _weak_raising_exact()[0]],
        "decoded_one_generation_indices": {"charged_left_doublet": left[0],
                                            "charge_conjugated_singlet": conjugate[0],
                                            "neutral_left_doublet": neutrino[0]},
        "incoming_electric_charge": -1,
        "momentum": list(p), "coordinate_momentum_energy": momentum_energy,
        "unprojected_spin_family_strength_trace": {"charged_left_doublet": left_trace,
                                                    "charge_conjugated_singlet_holes": hole_trace},
        "zero_Higgs_control": massless,
        "massive_native_trial": rows,
        "wrong_helicity_weight_formula": "(1-c|p|/sqrt(c^2|p|^2+m^2))/2; approaches m^2/(4c^2|p|^2) for m/(c|p|)->0",
        "complete_final_space_used_in_JdaggerJ": True,
        "scope": "a normalized current-vertex discriminator of the declared quadratic candidate, not a scattering rate, propagating weak gauge field, or observed-particle identification",
    }


def _norm(matrix):
    return float(np.linalg.norm(matrix, ord="fro"))


def _signature(matrix, tolerance=1e-10):
    values = np.linalg.eigvalsh(matrix)
    return {"positive": int(np.count_nonzero(values > tolerance)),
            "negative": int(np.count_nonzero(values < -tolerance)),
            "zero": int(np.count_nonzero(abs(values) <= tolerance))}


def report():
    """JSON-ready independent kinematic checks; floating residuals are diagnostics."""
    k = _kernel(None)
    p = (0.23, -0.31, 0.17)
    h, h0 = bdg_symbol(p, k), bdg_symbol((0, 0, 0), k)
    hminus = bdg_symbol(tuple(-x for x in p), k)
    n, d = len(k), 2 * len(k)
    zero, identity = np.zeros((d, d)), np.eye(d)
    exchange = np.block([[zero, identity], [identity, zero]])
    scale = coordinate_speed() ** 2 * math.fsum(x*x for x in p)
    expected_square = np.block([
        [scale*np.eye(d) + np.kron(np.eye(2), k.conj() @ k), zero],
        [zero, scale*np.eye(d) + np.kron(np.eye(2), k @ k.conj())]])
    currents = tuple(_numeric(t) for t in action.generators())
    current_residuals = []
    stabilizer_columns = []
    for t in currents:
        t_spin = np.kron(np.eye(2), t)
        kinetic = weyl_kinetic(p)
        current_residuals.append(_norm(t_spin @ kinetic - kinetic @ t_spin))
        derivative = t.T @ k + k @ t
        stabilizer_columns.append(np.r_[derivative.real.ravel(), derivative.imag.ravel()])
    stabilizer_rank = int(np.linalg.matrix_rank(np.array(stabilizer_columns).T, tol=1e-10))
    charge = nambu_charge()
    # Fixed coordinate interpolators, not selected energy eigenvectors. Their
    # dark sectors remain in the full spectral census and are not normalized
    # away. This is a supplied-action diagnostic, not detector tomography.
    interpolator_indices = [0, n, d, d+n]
    spectrum = charge_resolved_spectrum(h, charge, np.eye(2*d)[:, interpolator_indices])
    counts = {}
    for row in spectrum["clusters"]:
        key = str(round(sum(row["charge_interval"])/2, 12))
        counts[key] = counts.get(key, 0)+int(np.count_nonzero(row["frequencies"] > 1e-10))
    spectral_report = spectrum_summary(spectrum)
    spectral_report["fixed_interpolator_indices"] = interpolator_indices
    spectral_report["interpolator_design"] = (
        "first internal coordinate, both spin components and both Nambu blocks; "
        "canonical field coefficients with no energy-eigenvector selection; "
        "charged preparation/readout and a phase-reference apparatus are not supplied here")
    spin = nambu_spin()
    nodes = naive_weyl_nodes()
    return {
        "candidate": "declared all-left M15 x native F3 continuum-cutoff Weyl/Higgs symbol",
        "coordinate_speed": coordinate_speed(),
        "units": "declared spatial coordinate and programmed modular parameter",
        "dimensions": {"internal_left_Weyl_components": n, "spin_CAR_modes_at_one_momentum": d,
                       "Nambu_symbol": 2*d, "charged_Dirac_fields_counting_color": 21,
                       "massless_left_Weyl_neutrino_fields": 3, "independent_right_neutrinos": 0},
        "spectrum_at_nonzero_momentum": _signature(h),
        "spectrum_at_zero_momentum": _signature(h0),
        "positive_energy_electric_charge_multiplicities": counts,
        "charge_resolved_spectral_response": spectral_report,
        "massless_neutrino_internal_indices": action.ordering()["massless_neutrino_internal_family_indices"],
        "zero_momentum_caveat": "12 Nambu zeros encode six complex CAR zero modes of three left Weyl neutrinos; helicity at p=0 is undefined",
        "redundancy": "negative Nambu energies are conjugate partners at opposite momentum, not extra physical species or sterile nuR",
        "vacuum_caveat": "particle/antiparticle readout requires the quasifree ground-state reference, not an empty encoded Fock vacuum",
        "identities": {"Hermiticity_Frobenius": _norm(h-h.conj().T),
                       "dispersion_square_Frobenius": _norm(h@h-expected_square),
                       "particle_hole_Frobenius": _norm(exchange@h.conj()@exchange+hminus),
                       "electric_charge_commutator_Frobenius": _norm(charge@h-h@charge),
                       "max_kinetic_current_commutator_Frobenius": max(current_residuals),
                       "spin_Casimir_Frobenius": _norm(sum((s@s for s in spin), np.zeros_like(h))-.75*np.eye(2*d)),
                       "helicity_commutator_Frobenius": _norm(helicity(p)@h-h@helicity(p))},
        "fixed_neutral_Higgs_current_span": {"dimension": 12,
                                            "broken_rank_numeric": stabilizer_rank,
                                            "unbroken_dimension_numeric": 12-stabilizer_rank,
                                            "interpretation": "SU3 plus electromagnetic U1; full action covariance also transforms the Higgs"},
        "weak_chirality_discriminator": weak_chirality_report(p),
        "naive_lattice_negative_control": {"nodes": nodes, "left_nodes": sum(x["jacobian_chirality"]<0 for x in nodes),
                                           "right_nodes": sum(x["jacobian_chirality"]>0 for x in nodes),
                                           "net_chirality": sum(x["jacobian_chirality"] for x in nodes)},
        "spatial_contract": "finite menu of continuum momenta with a declared spectral cutoff; spatial attachment is supplied",
        "local_gauge_covariant_chiral_regulator_established": False,
        "measured_Standard_Model_particles_identified": False,
        "native_trial_masses_match_observed_hierarchy": False,
        "physical_seconds_GeV_or_electromagnetic_charge_normalization_selected": False,
    }
