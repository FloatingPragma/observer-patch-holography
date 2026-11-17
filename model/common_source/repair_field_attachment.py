"""Finite CAR attachment of a chosen native correlation/exchange architecture.

The native corrected cycle itself is still logically identity. This separate
architecture orbit-averages its existing two-hop correlation observable and
adds positive collective exchange between identical native carriers. Neither
control rule is selected by C/E repair. Their exact encoded action, rather
than a measured mass packet or prescribed logical hopping matrix, determines
the finite field considered here. Native total load remains two per carrier;
the conserved encoded occupation is a different, operationally read quantity.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import parameter_selection as selection, record_code as code
from .chiral_candidate import gaussian_source as gaussian
from .mass_response import history as instruments
from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..native_modular_reference.phase_inference import atan_interval
from .halfline import pi_interval
from ..scalar_observations.experimental_readout.stats import concentration_certificate
from ..calibrated_measurement.bounds import root_interval


ROOT = Path(__file__).resolve().parents[2]
SOURCE_GATE = 1e-8
LOCAL_MACRO_GATE = 1e-10
REPLAY_TOLERANCE = 2e-9
# Same-array diagnostic provenance only: target construction/SVD can round
# differently across backends. This is not an operator or physical allowance.
RETAINED_DIAGNOSTIC_TOLERANCE = 64*np.finfo(np.float64).eps
CELLS = 3
PROBE_TIME = F(1)
STEPS = 128
FIELD_TIME = F(8)
READOUT_RADIUS = F(1, 1000)
GAP_TIME = F(512)
TRIALS = 2_000_000
RECORD_RADIUS = F(1, 200)
RECORD_SOURCE_ALLOWANCE = F(1, 10**6)
PHYSICAL_SCORE_ALLOWANCE = F(1, 5000)
SAMPLER_ALLOWANCE = F(1, 10**7)


def _freeze(value):
    raw = np.asarray(value, dtype=np.complex128)
    return np.frombuffer(raw.tobytes(), np.complex128).reshape(raw.shape)


def _sector(value):
    if type(value) is not int or value not in (0, 1, 2):
        raise ValueError("one of the three declared carrier sectors required")
    return value


def _rational(value, *, positive=False):
    if (not (type(value) in (int, F) or isinstance(value, s.Rational)) or isinstance(value, bool)):
        raise ValueError("an original exact rational model parameter required")
    result = s.Rational(value)
    if result < 0 or (positive and result == 0):
        raise ValueError("a nonnegative exact parameter required")
    return result


@lru_cache(maxsize=1)
def _parameters():
    probe = selection.diagonal_cost(2)
    energies = tuple(probe["compression"][j, j] for j in range(6))
    vacuum, particle, gaps = [], [], []
    for sector in range(3):
        labels = 2*sector, 2*sector+1
        low, high = sorted(labels, key=lambda j: energies[j])
        vacuum.append(low); particle.append(high); gaps.append(energies[high]-energies[low])
    strength = sum(probe["coefficients"][a, b] for a in range(12) for b in range(a+1, 12))
    if tuple(gaps) != (s.Rational(1, 256), s.Rational(1, 1024), s.Rational(1, 3072)) or strength != s.Rational(125, 1024):
        raise ArithmeticError("the fixed native correlation normalization changed")
    return energies, tuple(vacuum), tuple(particle), tuple(gaps), strength


def logical_parameters():
    energies, vacuum, particle, gaps, strength = _parameters()
    return {"six_native_compression_energies": energies,
            "vacuum_labels": vacuum, "occupied_labels": particle,
            "onsite_baselines": tuple(energies[j] for j in vacuum),
            "positive_occupation_gaps": gaps, "exchange_strength": strength}


def occupation(cells=CELLS):
    if type(cells) is not int or not 1 <= cells <= 6:
        raise ValueError("one through six CAR sites required")
    return s.ImmutableMatrix(s.diag(*(mask.bit_count() for mask in range(1 << cells))))


def annihilator(site, cells=CELLS):
    occupation(cells)
    if type(site) is not int or not 0 <= site < cells:
        raise ValueError("a site in the declared CAR ordering required")
    value = s.zeros(1 << cells)
    for mask in range(1 << cells):
        if mask >> site & 1:
            value[mask ^ (1 << site), mask] = (-1)**((mask & ((1 << site)-1)).bit_count())
    return s.ImmutableMatrix(value)


def carrier_swap(a, b, cells=CELLS):
    """Ordinary exchange of entire distinguishable native carriers, not FSWAP."""
    occupation(cells)
    if type(a) is not int or type(b) is not int or not 0 <= a < b < cells:
        raise ValueError("two ordered distinct carrier sites required")
    matrix = s.zeros(1 << cells)
    for before in range(1 << cells):
        after = before ^ (((before >> a) ^ (before >> b)) & 1)*((1 << a)+(1 << b))
        matrix[after, before] = 1
    return s.ImmutableMatrix(matrix)


def chain_hamiltonian(sector=0, *, strength=F(1), cells=CELLS):
    """Exact open-chain8x8 at three sites, including its absolute baseline.

    ``strength`` multiplies all exchange bonds together. It is an optional
    declared geometry/control parameter, not a fitted species coefficient.
    The default J is fixed once from the native W² off-diagonal sum.
    """
    sector = _sector(sector)
    strength = _rational(strength)
    p = logical_parameters()
    dim = 1 << cells
    h = cells*p["onsite_baselines"][sector]*s.eye(dim)+p["positive_occupation_gaps"][sector]*occupation(cells)
    for site in range(cells-1):
        h += strength*p["exchange_strength"]*(s.eye(dim)-carrier_swap(site, site+1, cells))
    return s.ImmutableMatrix(h)


def current(bond, *, strength=F(1), cells=CELLS):
    if type(bond) is not int or not 0 <= bond < cells-1:
        raise ValueError("an oriented nearest-neighbor open-chain bond required")
    J = _rational(strength)*logical_parameters()["exchange_strength"]
    a, b = annihilator(bond, cells), annihilator(bond+1, cells)
    return s.ImmutableMatrix(-s.I*J*(a.T*b-b.T*a))


def _matrix_json(value):
    return [[str(item) for item in row] for row in s.Matrix(value).tolist()]


def native_obstruction():
    """Exact compression of every cross-copy two-port diagonal/swap term."""
    geometry = selection.code_geometry()
    for port in range(12):
        for support in code.SUPPORTS:
            if sum(state[port] for state in support) != 1:
                raise ArithmeticError("the native one-port scalar compression changed")
    # A cross-copy port swap changes the per-copy load2 unless its two
    # bits agree. In that case it fixes the entire native configuration.
    equal = (s.Rational(1, 6)**2+s.Rational(5, 6)**2)
    difference = 2*s.Rational(1, 6)-2*s.Rational(1, 6)**2
    if equal != s.Rational(13, 18) or difference != s.Rational(5, 18):
        raise ArithmeticError("cross-carrier scalar obstruction changed")
    return {"native_ports_checked": 12, "carrier_tensor_code_dimension": 36,
            "arbitrary_diagonal_cross_port_bilinear_compression": "I36/36",
            "arbitrary_cross_port_squared_difference_compression": str(difference)+"*I36",
            "arbitrary_single_cross_port_swap_compression": str(equal)+"*I36",
            "single_native_particle_transfer_compression": "0_36x36",
            "conclusion": "two-port cross-copy terms cannot provide a first-order encoded hopping generator on this code",
            "higher_order_elimination_or_collective_controls_excluded": False}


@lru_cache(maxsize=1)
def _orbit_data():
    """Actual full512 native observable and chosen matching-support cycle."""
    original = tuple(F(value) for value in selection._probe_energies())
    permutation = list(range(512))
    for support in code.SUPPORTS:
        addresses = [circuit.encode(state) for state in support]
        for a, b in zip(addresses, addresses[1:]+addresses[:1], strict=True):
            permutation[a] = b
    if sorted(permutation) != list(range(512)):
        raise ArithmeticError("native support twirl is not a complete permutation")
    average = []
    for before in range(512):
        address, total = before, F()
        for _ in range(6):
            total += original[address]
            address = permutation[address]
        if address != before:
            raise ArithmeticError("support-cycle period does not divide six")
        average.append(total/6)
    expected = _parameters()[0]
    for logical, support in enumerate(code.SUPPORTS):
        if {average[circuit.encode(state)] for state in support} != {F(expected[logical])}:
            raise ArithmeticError("finite orbit averaging did not preserve the native compression")
    return original, tuple(permutation), tuple(average)


def orbit_average():
    original, permutation, averaged = _orbit_data()
    changed = sum(a != b for a, b in zip(original, averaged, strict=True))
    return {"cycle": list(permutation), "native_original_energy": list(map(str, original)),
            "averaged_native_energy": list(map(str, averaged)),
            "changed_addresses": changed,
            "code_preserving_exactly": True,
            "code_six_level_compression": list(map(str, _parameters()[0])),
            "construction": "average(T^-r H2 T^r,r=0..5); T cycles each fixed matching support and fixes every other address",
            "actual_emitted_protocol": "direct Walsh compilation of the derived averaged full512 diagonal",
            "six_leg_equivalence": "equal complete endpoint unitary because all conjugated native potentials are diagonal and commute",
            "six_permutation_legs_or_intermediate_record_equivalence_executed": False,
            "unmodified_native_repair_selects_this_control": False}


def _walsh_diagonal_word(values, time):
    if len(values) != 512:
        raise ValueError("the entire existing nine-bit native carrier required")
    coefficients = list(values)
    stride = 1
    while stride < 512:
        for block in range(0, 512, 2*stride):
            for offset in range(stride):
                i, j = block+offset, block+offset+stride
                a, b = coefficients[i], coefficients[j]
                coefficients[i], coefficients[j] = a+b, a-b
        stride *= 2
    word = []
    for mask, coefficient in enumerate(coefficients):
        label = ["Z" if mask >> bit & 1 else "I" for bit in range(9)]
        word += compiler.pauli_rotation(label, -float(time*coefficient/512))
    return word


def onsite_source(time=PROBE_TIME):
    time = _rational(time)
    actual_code = code.code_matrix()
    word = _walsh_diagonal_word(_orbit_data()[2], time)
    output = compiler.execute(word, 9, actual_code)
    ideal = code.ideal_code_matrix()
    phases = np.exp(-1j*float(time)*np.asarray(list(map(float, _parameters()[0]))))
    target = ideal*phases[None, :]
    error = float(np.linalg.norm(output-target, 2))
    if error > SOURCE_GATE:
        raise ArithmeticError("native averaged correlation source exceeded its all-code-column gate")
    return {"output": _freeze(output), "encoded_input": actual_code,
            "word": tuple(word), "resources": compiler.resources(word),
            "operator_error_diagnostic": error, "code_columns_checked": 6,
            "complete_native_rows_retained": 512}


def _little_apply(matrix, sites, bits, value):
    """Tensor contraction of an already source-executed small complete macro.

    The compiler itself admits at most16 register bits. This factor method
    does not ask it to execute an18/19-bit tape, or replace local amplitudes
    by an ideal target. Matrix bitj is sites[j], including all scratch rows.
    """
    if (type(bits) is not int or not 1 <= bits <= 19 or len(set(sites)) != len(sites)
            or any(type(j) is not int or not 0 <= j < bits for j in sites)):
        raise ValueError("ordered distinct local macro sites in the finite register required")
    matrix, value = np.asarray(matrix), np.asarray(value)
    if matrix.shape != (1 << len(sites),)*2 or value.ndim != 2 or value.shape[0] != 1 << bits:
        raise ValueError("complete local macro and retained reference columns required")
    return compiler._apply(value, matrix, list(reversed(sites)), bits)


def _fredkin_target():
    value = np.zeros((8, 8), float)
    for before in range(8):
        after = before
        if before & 4 and ((before & 1) != ((before >> 1) & 1)):
            after ^= 3
        value[after, before] = 1
    return value


def fredkin_source():
    word = compiler.compile_real_orthogonal(_fredkin_target())
    matrix = compiler.execute(word, 3)
    error = float(np.linalg.norm(matrix-_fredkin_target(), 2))
    if error > LOCAL_MACRO_GATE:
        raise ArithmeticError("actual full-column controlled carrier-bit swap failed")
    return {"matrix": _freeze(matrix), "word": tuple(word), "resources": compiler.resources(word),
            "operator_error_diagnostic": error}


def _collective_outputs(encoded, macro, H, phase, scalar):
    incoming = np.kron(encoded, encoded)
    value = np.zeros((1 << 19, 4), complex)
    value[:1 << 18] = incoming
    for stage in (0, 1):
        value = _little_apply(H, (18,), 19, value)
        for bit in range(9):
            value = _little_apply(macro, (bit, 9+bit, 18), 19, value)
        value = _little_apply(H, (18,), 19, value)
        if stage == 0:
            value = _little_apply(phase, (18,), 19, value)
    value = _little_apply(scalar, (18,), 19, value)
    return incoming, value[:1 << 18], value[1 << 18:]


def _collective_local_gates(macro, H, phase, scalar, theta=None):
    if theta is None:
        theta = float(PROBE_TIME*logical_parameters()["exchange_strength"])
    targets = (_fredkin_target(), 1j*np.array([[1, 1], [1, -1]])/math.sqrt(2),
               np.diag([np.exp(1j*theta), np.exp(-1j*theta)]), np.exp(-1j*theta)*np.eye(2))
    for actual, expected in zip((macro, H, phase, scalar), targets, strict=True):
        if (actual.shape != expected.shape or not np.isfinite(actual).all()
                or np.linalg.norm(actual-expected, 2) > LOCAL_MACRO_GATE):
            raise ValueError("retained collective local macro exceeds its all-input admission")


def native_collective_exchange(time=PROBE_TIME, *, sector=0):
    """Four unknown two-mode Fock inputs through the actual19-bit factor word.

    Nine source-executed Fredkin blocks form a controlled swap of the two
    full nine-bit native carriers. A retained ancilla implements the exact
    involution exponential; its nonzero failed-return output is not removed.
    The complete outputs are represented by input/local-macro factors and
    this fixed chronological grammar, rather than an archived dense vector.
    """
    sector, time = _sector(sector), _rational(time)
    parameters = logical_parameters()
    labels = parameters["vacuum_labels"][sector], parameters["occupied_labels"][sector]
    encoded = code.code_matrix()[:, labels]
    # Little endian: carrier0 bits0..8, carrier1 bits9..17, ancilla18.
    macro = fredkin_source()
    Hword = [compiler.pulse("H", math.pi/2, [0])]
    H = compiler.execute(Hword, 1)  # iH, whose four occurrences give phase1.
    theta = time*parameters["exchange_strength"]
    phase_word = [compiler.pulse("Z", float(theta), [0])]
    scalar_word = [compiler.pulse("I", -float(theta), [0])]
    phase, scalar = compiler.execute(phase_word, 1), compiler.execute(scalar_word, 1)
    _collective_local_gates(macro["matrix"], H, phase, scalar, float(theta))
    incoming, returned, failed = _collective_outputs(encoded, macro["matrix"], H, phase, scalar)
    swapped = incoming.reshape(512, 512, 4).transpose(1, 0, 2).reshape(1 << 18, 4)
    target = np.exp(-1j*float(theta))*(math.cos(float(theta))*incoming+1j*math.sin(float(theta))*swapped)
    operator = float(np.linalg.norm(returned-target, 2))
    failure = float(np.linalg.norm(failed, 2))
    if max(operator, failure) > SOURCE_GATE:
        raise ArithmeticError("source collective exchange or retained ancilla return exceeded gate")
    projection = incoming.conj().T@returned
    target_logical = expm(-1j*float(theta)*np.asarray(s.eye(4)-carrier_swap(0, 1, 2), complex))
    if np.linalg.norm(projection-target_logical, 2) > SOURCE_GATE:
        raise ArithmeticError("whole-carrier source does not restrict to ordinary logical swap")
    pulses = 18*macro["resources"]["native_pulses"]+6
    return {"encoded_two_columns": _freeze(encoded), "fredkin_matrix": macro["matrix"],
            "ancilla_H": _freeze(H), "ancilla_Z": _freeze(phase), "scalar_I": _freeze(scalar),
            "logical_projection": _freeze(projection),
            "source_operator_error_diagnostic": operator, "failed_return_operator_error_diagnostic": failure,
            "chronological_factor_rule": "H18;product_jFredkin(j,9+j,18);H18;Z18(theta);H18;product_jFredkin(j,9+j,18);H18;I18(-theta)",
            "retained_full_native_output_rows": 1 << 19, "unknown_Fock_columns": 4,
            "native_pulses_per_collective_bond": pulses,
            "two_carrier_encoder_pulses": 2*code.encoding_resources()["native_pulses"],
            "proper_load_export_events_per_collective_bond": 2*pulses,
            "local_all_input_numerical_admission": "1/10000000000",
            "complete_macro_operator_bound_conditional": str(24*F(1, 10**10)/(1-24*F(1, 10**10))),
            "local_admission_is_physical_or_outward_interval_certificate": False,
            "factorized_contraction_not_direct19qubit_compiler_call": True,
            "unnormalized_failed_ancilla_output_retained_by_factor_rule": True}


def exact_field_checks():
    """All Fock sectors, continuity, stable vacuum and independent species."""
    n = occupation()
    identity = s.eye(8)
    a = [annihilator(j) for j in range(3)]
    for j in range(3):
        for k in range(3):
            if a[j]*a[k]+a[k]*a[j] != s.zeros(8) or a[j]*a[k].T+a[k].T*a[j] != (identity if j == k else s.zeros(8)):
                raise ArithmeticError("independently occupiable CAR identities failed")
    report = []
    p = logical_parameters()
    for sector in range(3):
        h = chain_hamiltonian(sector)
        if h*n != n*h:
            raise ArithmeticError("derived native collective action does not conserve encoded charge")
        reconstructed = 3*p["onsite_baselines"][sector]*identity+p["positive_occupation_gaps"][sector]*n
        for j in range(2):
            x, y = a[j].T*a[j], a[j+1].T*a[j+1]
            reconstructed += p["exchange_strength"]*(x+y-2*x*y-a[j].T*a[j+1]-a[j+1].T*a[j])
        if h != reconstructed:
            raise ArithmeticError("native carrier exchange is not the full interacting CAR operator")
        for site in range(3):
            density = a[site].T*a[site]
            expected = (current(site-1) if site else s.zeros(8))-(current(site) if site < 2 else s.zeros(8))
            if s.I*(h*density-density*h) != expected:
                raise ArithmeticError("encoded charge/current continuity failed")
        threshold = p["positive_occupation_gaps"][sector]
        for count in range(4):
            support = [mask for mask in range(8) if mask.bit_count() == count]
            block = h.extract(support, support)-3*p["onsite_baselines"][sector]*s.eye(len(support))
            eigenvalues = block.eigenvals()
            if min(eigenvalues) != count*threshold:
                raise ArithmeticError("positive exchange did not retain its stable vacuum and charge thresholds")
        one = h.extract((1, 2, 4), (1, 2, 4))-3*p["onsite_baselines"][sector]*s.eye(3)
        expected = {threshold: 1, threshold+p["exchange_strength"]: 1, threshold+3*p["exchange_strength"]: 1}
        if one.eigenvals() != expected:
            raise ArithmeticError("actual native-derived one-particle band changed")
        report.append({"sector": sector, "native_vacuum_label": p["vacuum_labels"][sector],
            "native_occupied_label": p["occupied_labels"][sector],
            "onsite_baseline": str(p["onsite_baselines"][sector]), "positive_charge_threshold": str(threshold),
            "one_particle_poles_above_vacuum": list(map(str, sorted(expected))),
            "minimum_energy_above_vacuum_by_charge": [str(j*threshold) for j in range(4)],
            "Hamiltonian_all8_columns": _matrix_json(h)})
    # Nonadjacent ordinary carrier exchange is not local CAR hopping.
    x, y = a[0].T*a[0], a[2].T*a[2]
    false_local = x+y-2*x*y-a[0].T*a[2]-a[2].T*a[0]
    discrepancy = s.eye(8)-carrier_swap(0, 2)-false_local
    if discrepancy[6, 3] != -2:
        raise ArithmeticError("the occupied-middle-mode loop countercontrol vanished")
    two = np.asarray(chain_hamiltonian(0).extract((3, 5, 6), (3, 5, 6)), complex)
    # At N2 the interaction yields shifts not equal to sums of distinct
    # one-particle free-fermion energies; the fully occupied state is also
    # unaffected by exchange, instead of accumulating a free-band trace.
    return {"per_chain": report, "independent_chain_count": 3,
        "independent_carriers_per_chain": 3, "total_independent_CAR_modes": 9,
        "full_tensor_Fock_dimension": 512, "inter_chain_mixing": "0",
        "unavoidable_density_interaction_per_bond": str(-2*p["exchange_strength"]),
        "two_particle_is_free_sum_of_one_particle_poles": False,
        "three_particle_exchange_energy": "0",
        "nonadjacent_naive_local_hopping_counterexample": {"before_mask": 3, "after_mask": 6,
            "operator_difference": "-2", "missing_service": "Jordan-Wigner parity transport or an additional local gauge encoding on graphs with loops"},
        "dispersion": "open3cell graph Laplacian: gap+2J*(1-cos(k)), k=0,pi/3,2pi/3",
        "relativistic_Weyl_Dirac_or_3D_local_CAR_dispersion": False}


def global_car_certificate():
    """All512basis states, including cross-chain odd-operator parity signs."""
    def apply(mask, site, dagger):
        occupied = (mask >> site) & 1
        if occupied == int(dagger):
            return None
        return mask ^ (1 << site), (-1)**((mask & ((1 << site)-1)).bit_count())
    def composition(mask, left, right, dagger):
        first = apply(mask, right, dagger)
        if first is None:
            return None
        second = apply(first[0], left, False)
        return None if second is None else (second[0], first[1]*second[1])
    checked = 0
    for mask in range(512):
        for left in range(9):
            for right in range(9):
                for dagger in (False, True):
                    values = {}
                    first = composition(mask, left, right, dagger)
                    if first is not None:
                        values[first[0]] = first[1]
                    # Reverse order is a_right^(dagger) a_left.
                    first_step = apply(mask, left, False)
                    if first_step is not None:
                        second = apply(first_step[0], right, dagger)
                        if second is not None:
                            values[second[0]] = values.get(second[0], 0)+first_step[1]*second[1]
                    values = {key: value for key, value in values.items() if value}
                    expected = {mask: 1} if dagger and left == right else {}
                    if values != expected:
                        raise ArithmeticError("global graded tensor field fails CAR")
                    checked += 1
    return {"ordering": "mode=3*chain+site; a_chain,site includes parity of every earlier chain",
            "independently_occupiable_modes": 9, "basis_states": 512,
            "all_basis_anticommutators_checked": checked,
            "ungraded_tensor_odd_operators_used": False,
            "even_chain_Hamiltonians_commute": True,
            "inter_chain_mixing": "0", "three_chains_identified_with_SM_generations": False}


def exchange_word(duration):
    """Exact positive ordinary carrier-swap bond in its proven code image."""
    duration = _rational(duration)
    theta = float(duration*logical_parameters()["exchange_strength"])/2
    return ([compiler.pulse("I", -theta, [0])]
        + compiler.pauli_rotation("XX", theta)
        + compiler.pauli_rotation("YY", theta)
        + compiler.pauli_rotation("ZZ", theta))


def field_source(sector=0, *, time=FIELD_TIME, steps=STEPS):
    """Actual all8-column Strang word; finite splitting error is retained."""
    sector, time = _sector(sector), _rational(time)
    if type(steps) is not int or not 1 <= steps <= 4096 or time > 32:
        raise ValueError("one through4096 steps and model duration at most32 required")
    dt = time/steps
    first, second = exchange_word(dt/2), exchange_word(dt)
    local_first, local_second = compiler.execute(first, 2), compiler.execute(second, 2)
    target_first = expm(-1j*float(dt/2*logical_parameters()["exchange_strength"])*np.asarray(s.eye(4)-carrier_swap(0, 1, 2), complex))
    target_second = expm(-1j*float(dt*logical_parameters()["exchange_strength"])*np.asarray(s.eye(4)-carrier_swap(0, 1, 2), complex))
    local_error = max(np.linalg.norm(local_first-target_first, 2), np.linalg.norm(local_second-target_second, 2))
    if local_error > SOURCE_GATE:
        raise ArithmeticError("source reduced bond differs from the independently derived native-carrier interaction")
    p = logical_parameters()
    word = [compiler.pulse("I", -float(3*time*p["onsite_baselines"][sector]), [0])]
    for site in range(3):
        mu = float(time*p["positive_occupation_gaps"][sector])/2
        word += [compiler.pulse("I", -mu, [site]), compiler.pulse("Z", mu, [site])]
    for _ in range(steps):
        word += compiler.remap(first, (0, 1))+compiler.remap(second, (1, 2))+compiler.remap(first, (0, 1))
    actual = _freeze(compiler.execute(word, 3))
    h = np.asarray(chain_hamiltonian(sector), complex)
    target = expm(-1j*float(time)*h)
    # Childs et al., PRX11,011020(2021), Proposition16 gives coefficients
    #1/24 and1/12 for these two nested commutators. We use twice that bound.
    J = float(p["exchange_strength"])
    A, B = J*np.asarray(s.eye(8)-carrier_swap(0, 1), complex), J*np.asarray(s.eye(8)-carrier_swap(1, 2), complex)
    comm = A@B-B@A
    nested = np.linalg.norm(A@comm-comm@A, 2)+2*np.linalg.norm(B@comm-comm@B, 2)
    A0, B0 = s.eye(8)-carrier_swap(0, 1), s.eye(8)-carrier_swap(1, 2)
    C0 = A0*B0-B0*A0
    for D in (A0*C0-C0*A0, B0*C0-C0*B0):
        D2 = D*D
        if D != D.T or D2*(D2-12*s.eye(8)) != s.zeros(8) or s.trace(D2) != 48:
            raise ArithmeticError("exact nested-commutator norm certificate changed")
    # Each nested norm is2sqrt3*J³. Directed rational root, not a LAPACK
    # residual, supplies the portable mathematical splitting upper bound.
    split_exact = root_interval(F(3))[1]*F(p["exchange_strength"])**3*F(time)**3/(2*steps**2)
    split_bound = float(split_exact)
    error = float(np.linalg.norm(actual-target, 2))
    if error > split_bound+3*steps*SOURCE_GATE+SOURCE_GATE:
        raise ArithmeticError("source finite splitting exceeds its declared comparison budget")
    return {"unitary": actual, "word": tuple(word), "resources": compiler.resources(word),
            "half_bond_unitary": _freeze(local_first), "full_bond_unitary": _freeze(local_second),
            "splitting_operator_bound_diagnostic": split_bound,
            "splitting_operator_bound_exact_upper": str(split_exact),
            "splitting_primary_source": "https://arxiv.org/abs/1912.08854; Proposition16; doubled two-term bound",
            "dense_target_difference_diagnostic": error, "local_source_difference_diagnostic": float(local_error),
            "source_steps": steps, "model_coordinate_duration": str(time)}


def native_charge_readout():
    """Decode, read the chosen matching charge, re-encode all512 rows."""
    encoded = code.code_matrix()
    decoded = code._apply_permutation(code._encoder_permutation(), encoded, True)
    decoded = compiler.execute(code._mixing_word(6, True), 9, decoded)
    p = logical_parameters()
    charge = [0]*8
    for label in p["occupied_labels"]:
        charge[label] = 1
    read = instruments.diagonal_read(decoded, [1-value for value in charge], targets=(0, 1, 2))
    outputs = tuple(_freeze(code._encode(value)) for value in read.outputs)
    ideal = code.ideal_code_matrix()
    errors = [float(np.linalg.norm(value-ideal*np.array([int(charge[j] == result) for j in range(6)])[None, :], 2))
              for result, value in enumerate(outputs)]
    if max(errors) > SOURCE_GATE:
        raise ArithmeticError("native matching-support charge readout exceeded its full-column source gate")
    return {"outputs": outputs, "charge_labels": charge[:6],
            "pointer_resources": dict(read.resources), "decoder_resources": code.decoding_resources(),
            "encoder_resources": code.encoding_resources(), "operator_error_diagnostic": max(errors),
            "native_one_port_load_is_the_read_charge": False,
            "all_pointer_outcomes_syndromes_and_reference_columns_retained": True}


def two_time_records(unitary):
    """One continuing two-density-read history, with complete state updates."""
    prep_word = gaussian.local_gate("hopping", 1., -math.pi/4).word
    initial = np.zeros((8, 1), complex); initial[1] = 1
    initial = compiler.execute(compiler.remap(prep_word, (0, 1)), 3, initial)
    before = unitary@initial
    first = instruments.pauli_read(before, 0, axis="Z")
    blocks = []
    for bit, branch in enumerate(first.outputs):
        second = instruments.pauli_read(unitary@branch, 1, axis="Z")
        for final, output in enumerate(second.outputs):
            blocks.append((bit, final, output))
    arrays = np.stack([value for _, _, value in blocks])
    probabilities = [float(np.vdot(value, value).real) for _, _, value in blocks]
    mass = sum(probabilities)
    if abs(mass-1) > SOURCE_GATE:
        raise ArithmeticError("complete continuing field instrument lost probability mass")
    joint = np.array(probabilities).reshape(2, 2)
    gram = sum(value@value.conj().T for _, _, value in blocks)
    # Both observers' data are restrictions of this one instrument. Keeping
    # only one observer's reads changes backaction; it is not this marginal.
    return ({"joint_density_records": probabilities,
            "observer_A_first_record_marginal": joint.sum(axis=1).tolist(),
            "observer_B_second_record_marginal": joint.sum(axis=0).tolist(),
            "records_are_common_history_restrictions": True,
            "first_measurement_backaction_retained": True,
            "no_first_read_protocol_claimed_equal": False,
            "initial_state": "source hopping of seed encoded charge1: (|100>+i|010>)/sqrt2",
            "total_encoded_charge_before_and_after": 1,
            "source_preparation_native_pulses": compiler.resources(prep_word)["native_pulses"],
            "source_pointer_native_pulses": first.resources["native_pulses"]+second.resources["native_pulses"]},
            {"initial": _freeze(initial), "complete_two_time_outputs": arrays, "final_unread_density": _freeze(gram)})


def _route_two_mode_word(word, a, b):
    if not 0 <= a < b < 4:
        raise ValueError("ordered modes within the four-mode probe required")
    route = [(j-1, j) for j in range(b, a+1, -1)]
    swap = gaussian.local_gate("fswap").word
    result = [event for sites in route for event in compiler.remap(swap, sites)]
    result += compiler.remap(word, (a, a+1))
    result += [event for sites in reversed(route) for event in compiler.remap(swap, sites)]
    return result


def _mode_word(mode):
    if type(mode) is not int or mode not in (0, 1, 2):
        raise ValueError("one of the three fixed spatial standing-wave probes required")
    vectors = ((1, 1, 1), (1, 0, -1), (1, -2, 1))
    vector = np.asarray(vectors[mode], float)
    vector /= np.linalg.norm(vector)
    a, b, c = vector
    alpha, beta = math.atan2(c, math.hypot(a, b)), math.atan2(b, a)
    first = _route_two_mode_word(gaussian.local_gate("hopping", -1j, alpha).word, 0, 2)
    second = compiler.remap(gaussian.local_gate("hopping", -1j, beta).word, (0, 1))
    word = first+second
    actual = compiler.execute(word, 4)
    one = np.zeros(16, complex); one[1] = 1
    expected = np.zeros(16, complex)
    expected[[1, 2, 4]] = vector
    if np.linalg.norm(actual@one-expected) > SOURCE_GATE:
        raise ArithmeticError("source spatial mode preparation does not match the fixed graph probe")
    return word, _freeze(actual)


def _phase_box(cosine, sine, radius=READOUT_RADIUS):
    """Directed atan2 in the independently declared positive(0,pi) window."""
    c, q = (cosine-radius, cosine+radius), (sine-radius, sine+radius)
    if q[0] <= 0:
        raise ValueError("charged phase does not resolve its declared positive half-circle")
    if c[0] > 0:
        return atan_interval(q[0]/c[1])[0], atan_interval(q[1]/c[0])[1]
    if c[1] < 0:
        low, high = atan_interval(q[0]/(-c[0]))[0], atan_interval(q[1]/(-c[1]))[1]
        pi = pi_interval()
        return pi[0]-high, pi[1]-low
    raise ValueError("this finite probe refuses unresolved quadrant boundaries")


def _binary_probability(value):
    return sum((F(float(z.real))**2+F(float(z.imag))**2 for z in np.asarray(value).ravel()), F())


def charged_phase_records(unitary, *, time=FIELD_TIME, modes=(0, 1, 2)):
    """Actual fixed spatial-mode/source-reference quadrature instruments.

    The state has exactly one encoded charge shared with the reference.
    Preparation/analyzers depend only on the fixed graph mode, never on a
    gap, eigenvalue or measured particle packet. Read means are obtained from
    original unnormalized binary source outputs, not a substituted frequency.
    The interval is conditional on the stated numerical score admission;
    finite shots and physical control accuracy are not supplied here.
    """
    time = _rational(time, positive=True)
    results, retained = [], []
    for mode in modes:
        mode_word, mode_unitary = _mode_word(mode)
        balance = _route_two_mode_word(gaussian.local_gate("hopping", -1j, -math.pi/4).word, 0, 3)
        seed = np.zeros((16, 1), complex); seed[8] = 1
        prep_word = balance+mode_word
        initial = compiler.execute(prep_word, 4, seed)
        evolved = np.kron(np.eye(2), unitary)@initial
        values, output_rows = [], []
        costs = []
        for axis in ("cosine", "sine"):
            coefficient, angle = (-1j, math.pi/4) if axis == "cosine" else (1., -math.pi/4)
            analyzer = _route_two_mode_word(gaussian.local_gate("hopping", coefficient, angle).word, 0, 3)
            word = compiler.inverse(mode_word)+analyzer
            changed = compiler.execute(word, 4, evolved)
            read = instruments.pauli_read(changed, 0, axis="Z")
            probabilities = tuple(_binary_probability(value) for value in read.outputs)
            if abs(sum(probabilities)-1) > F(1, 10**7):
                raise ArithmeticError("charged reference read lost complete probability mass")
            values.append(probabilities[0]-probabilities[1])
            output_rows.append(np.stack(read.outputs))
            costs.append(compiler.resources(word)["native_pulses"]+read.resources["native_pulses"])
        box = _phase_box(*values)
        energy_box = box[0]/F(time), box[1]/F(time)
        results.append({"spatial_mode": mode, "raw_cosine_mean": str(values[0]), "raw_sine_mean": str(values[1]),
            "conditional_record_pole_interval": list(map(str, energy_box)),
            "conditional_record_pole_interval_numeric": list(map(float, energy_box)),
            "source_preparation_native_pulses": compiler.resources(prep_word)["native_pulses"],
            "source_analyzer_pointer_native_pulses": costs,
            "frequency_or_gap_received_by_inverse": False,
            "source_reference_total_encoded_charge": 1})
        retained.append(np.stack(output_rows))
    return results, np.stack(retained)


def uniform_gap_source(sector):
    """Long fixed-time word exact on vacuum and uniform one-particle mode.

    Each positive exchange bond separately annihilates both states, hence
    the ordered word needs no splitting limit for this particular probe.
    Its complete8-column complement is retained but is NOT claimed to equal
    the long-time simultaneous-chain evolution on other initial states.
    """
    sector = _sector(sector)
    p = logical_parameters()
    uniform = s.zeros(8, 1)
    uniform[1] = uniform[2] = uniform[4] = 1
    for bond in range(2):
        if (s.eye(8)-carrier_swap(bond, bond+1))*uniform != s.zeros(8, 1):
            raise ArithmeticError("uniform charged mode is not fixed by every bond")
    word = [compiler.pulse("I", -float(3*GAP_TIME*p["onsite_baselines"][sector]), [0])]
    mu = float(GAP_TIME*p["positive_occupation_gaps"][sector])/2
    for site in range(3):
        word += [compiler.pulse("I", -mu, [site]), compiler.pulse("Z", mu, [site])]
    for bond in range(2):
        word += compiler.remap(exchange_word(GAP_TIME), (bond, bond+1))
    return {"unitary": _freeze(compiler.execute(word, 3)), "word": tuple(word),
            "resources": compiler.resources(word), "splitting_error_on_declared_uniform_probe": "0",
            "whole_complement_long_time_chain_equality_claimed": False}


def uniform_gap_prior():
    """Coarse operator-wide non-alias admission, independent of read records.

    The six native matching means all lie in[18/512,21/512]. Every positive
    oriented difference is therefore at most3/512. Every exchange bond
    annihilates the uniform mode, so at the common programmed duration512
    its phase is in(0,3], strictly below pi. This is one common spectral
    diameter bound, not a separate exact gap supplied to the inverse.
    """
    lo, hi = F(18, 512), F(21, 512)
    if any(not lo <= F(value) <= hi for value in _parameters()[0]) or pi_interval()[0] <= 3:
        raise ArithmeticError("the native architecture no longer admits the common non-alias window")
    return {"all_six_native_mean_energy_enclosure": [str(lo), str(hi)],
            "positive_uniform_charge_pole_prior": ["0", "3/512"],
            "common_programmed_duration": "512", "positive_phase_prior": ["0", "3"],
            "pi_lower_strictly_exceeds_three": True,
            "derived_before_and_without_measured_record_counts": True,
            "separate_exact_gaps_given_to_inverse": False}


def decode_uniform_records(counts):
    """Count-only inverse, with no native energies/operator/packet access."""
    if (type(counts) is not np.ndarray or counts.dtype != np.dtype("int64") or counts.shape != (3, 2, 2)
            or np.any(counts < 0) or np.any(counts > TRIALS) or np.any(counts.sum(axis=2) != TRIALS)):
        raise ValueError("complete original3chain/2quadrature/binary int64 counts required")
    radius = RECORD_RADIUS+RECORD_SOURCE_ALLOWANCE+PHYSICAL_SCORE_ALLOWANCE+SAMPLER_ALLOWANCE
    boxes = []
    for values in counts:
        means = [F(int(row[0])-int(row[1]), TRIALS) for row in values]
        phase = _phase_box(*means, radius=radius)
        boxes.append((phase[0]/GAP_TIME, phase[1]/GAP_TIME))
    return tuple(boxes)


def uniform_gap_records(*, retained_counts=None, seed=20261029):
    if type(seed) is not int or seed < 0:
        raise ValueError("an original nonnegative simulated-record seed required")
    if retained_counts is not None:
        decode_uniform_records(retained_counts)
    prior = uniform_gap_prior()
    generator = np.random.Generator(np.random.PCG64(seed))
    source, outputs, counts, rows, costs = [], [], [], [], []
    for sector in range(3):
        value = uniform_gap_source(sector)
        diagnostic, complete = charged_phase_records(value["unitary"], time=GAP_TIME, modes=(0,))
        branch_outputs = complete[0]
        sector_counts = []
        for setting in range(2):
            raw = [_binary_probability(branch) for branch in branch_outputs[setting]]
            total = sum(raw)
            if abs(total-1) > SAMPLER_ALLOWANCE:
                raise ArithmeticError("uniform-gap source total exceeds sampler normalization budget")
            law = np.asarray([float(q/total) for q in raw])
            law /= law.sum()
            sector_counts.append(generator.multinomial(TRIALS, law) if retained_counts is None else retained_counts[sector, setting])
        counts.append(np.stack(sector_counts)); source.append(value["unitary"]); outputs.append(branch_outputs)
        row = diagnostic[0]
        phase = float(GAP_TIME*logical_parameters()["positive_occupation_gaps"][sector])
        if max(abs(float(F(row["raw_cosine_mean"]))-math.cos(phase)),
               abs(float(F(row["raw_sine_mean"]))-math.sin(phase))) > float(RECORD_SOURCE_ALLOWANCE):
            raise ArithmeticError("uniform native source/reference score exceeds its numerical admission")
        rows.append(_phase_diagnostic(row)); costs.append(value["resources"]["native_pulses"])
    counts = np.stack(counts)
    intervals = decode_uniform_records(counts)
    gaps = logical_parameters()["positive_occupation_gaps"]
    for sector, box in enumerate(intervals):
        if not box[0] <= F(gaps[sector]) <= box[1]:
            raise ArithmeticError("record-derived charged gap missed the native architecture")
    certificate = concentration_certificate(6, TRIALS, math.nextafter(float(RECORD_RADIUS), 0.))
    if F(certificate["familywise_error_upper"]) >= F(1, 100):
        raise ArithmeticError("uniform-gap finite-sampling certificate failed")
    return {"simulated_all_attempt_records": 6*TRIALS, "trials_per_quadrature": TRIALS,
            "fixed_model_coordinate_duration": str(GAP_TIME), "sampler_seed": seed,
            "non_alias_prior_from_native_architecture": prior,
            "clock": "one supplied programmed native-coupling coordinate; no physical clock dictionary derived",
            "probe_parameters_use_empirical_masses": False,
            "all_chains_use_the_same_duration_and_apparatus": True,
            "charge_reference_retained": True, "unknown_positive_common_visibility_cancels_in_phase_ratio": True,
            "inverse_receives_native_gaps_or_Hamiltonian": False,
            "conditional_record_gap_intervals": [list(map(str, box)) for box in intervals],
            "conditional_record_gap_intervals_numeric": [list(map(float, box)) for box in intervals],
            "raw_source_probability_diagnostics": rows,
            "uniform_mode_splitting_error": "0",
            "source_score_allowance_assumed_numerically_checked": str(RECORD_SOURCE_ALLOWANCE),
            "physical_score_allowance_required_not_authenticated": str(PHYSICAL_SCORE_ALLOWANCE),
            "sampler_normalization_allowance": str(SAMPLER_ALLOWANCE),
            "record_radius": str(RECORD_RADIUS), "simultaneous_sampling_certificate": certificate,
            "reduced_source_field_pulses_per_setting_by_chain": costs}, {
            "uniform_gap_source_unitaries": np.stack(source), "uniform_gap_reference_pointer_outputs": np.stack(outputs),
            "uniform_gap_complete_simulated_counts": counts}


def _phase_diagnostic(row):
    result = deepcopy(row)
    for key in ("raw_cosine_mean", "raw_sine_mean"):
        result[key+"_diagnostic"] = float(F(result.pop(key)))
    result.pop("conditional_record_pole_interval")
    return result


def execute(*, retained_counts=None):
    algebra, global_car = exact_field_checks(), global_car_certificate()
    onsite = onsite_source()
    native_bond = native_collective_exchange()
    charge = native_charge_readout()
    arrays = {"native_encoded_six_columns": onsite["encoded_input"],
              "native_orbit_averaged_outputs": onsite["output"],
              "native_charge_pointer_outputs": np.stack(charge["outputs"]),
              "collective_encoded_two_columns": native_bond["encoded_two_columns"],
              "collective_fredkin_macro": native_bond["fredkin_matrix"],
              "collective_ancilla_H": native_bond["ancilla_H"], "collective_ancilla_Z": native_bond["ancilla_Z"],
              "collective_scalar_I": native_bond["scalar_I"], "collective_logical_projection": native_bond["logical_projection"]}
    field_reports, field_outputs, histories, phase_outputs, two_time = [], [], [], [], []
    half_bond, full_bond = None, None
    p = logical_parameters()
    for sector in range(3):
        source = field_source(sector)
        history, retained = two_time_records(source["unitary"])
        phases, phase_data = charged_phase_records(source["unitary"])
        expected = algebra["per_chain"][sector]["one_particle_poles_above_vacuum"]
        for row, pole in zip(phases, expected, strict=True):
            low, high = map(F, row["conditional_record_pole_interval"])
            if not low <= F(pole) <= high:
                raise ArithmeticError("source/reference pole readout missed the native-derived dispersion")
        field_reports.append({key: value for key, value in source.items() if key not in
            ("unitary", "word", "half_bond_unitary", "full_bond_unitary")}
            | {"sector": sector, "source_reference_pole_diagnostics": [_phase_diagnostic(row) for row in phases]})
        field_outputs.append(source["unitary"]); histories.append(retained["complete_two_time_outputs"])
        phase_outputs.append(phase_data); two_time.append(history)
        half_bond, full_bond = source["half_bond_unitary"], source["full_bond_unitary"]
    gaps, gap_arrays = uniform_gap_records(retained_counts=retained_counts)
    arrays.update(field_source_unitaries=np.stack(field_outputs), complete_two_time_density_outputs=np.stack(histories),
                  charged_reference_pointer_outputs=np.stack(phase_outputs),
                  half_bond_source_unitary=half_bond, full_bond_source_unitary=full_bond)
    arrays.update(gap_arrays)
    report = {"schema": "oph.repair_field_attachment.v1",
        "classification": "native-derived finite interacting CAR field under an explicitly chosen correlation-average/collective-exchange architecture",
        "compatibility_domain": {"native_carrier": "same total-two binary66 sector in original512-row code register",
            "legitimate_control_boundary": "tensor product of record_code's six-dimensional matching codes",
            "controls": "completed averaged-potential, whole-carrier exchange, decoded matching-charge pointer",
            "repair_schedule": "ideal complete C/E repair/recovery may interleave only at completed compatible control boundaries",
            "arbitrary_concurrent_repair_inside_source_words_covered": False,
            "finite_source_errors_and_physical_CE_protection_authenticated": False},
        "native_port_local_obstruction": native_obstruction(), "finite_correlation_average": orbit_average(),
        "field_algebra": algebra, "global_graded_CAR": global_car,
        "operational_charge": {"logical_eigenvalues": charge["charge_labels"],
            "definition": "matching-support occupation after stated hole orientations; read by existing inverse encoder, pointer, encoder",
            "native_total_load_per_carrier": "2*I", "native_port_load_compression": "I6/6",
            "electric_color_weak_charge_or_SM_generation_identified": False,
            "complete_native_pointer_outputs_retained": True},
        "collective_native_source": {key: value for key, value in native_bond.items() if key not in
            ("encoded_two_columns", "fredkin_matrix", "ancilla_H", "ancilla_Z", "scalar_I", "logical_projection")},
        "native_onsite_source": {"model_coordinate_duration": str(PROBE_TIME), "resources": onsite["resources"],
            "operator_error_diagnostic": onsite["operator_error_diagnostic"], "complete_unknown_columns": 6, "complete_rows": 512},
        "native_charge_source": {"pointer_resources": charge["pointer_resources"], "decoder_resources": charge["decoder_resources"],
            "encoder_resources": charge["encoder_resources"], "operator_error_diagnostic": charge["operator_error_diagnostic"]},
        "finite_field_source": field_reports, "continuing_density_observers": two_time,
        "uniform_charged_gap_records": gaps,
        "source_contract": {"numerical_source_operator_admission": "1/100000000", "fresh_complete_array_tolerance": "1/500000000",
            "uniform_gap_score_admission": str(RECORD_SOURCE_ALLOWANCE),
            "native_collective_and_reduced_bond_words_identical": False,
            "their_ideal_code_endpoint_actions_identical": True,
            "reduced_field_arrays_are_literal_source_executed": True,
            "full_native_field_numerical_composition_error_certified": False,
            "native_19bit_outputs_retained": "exact factorization from retained512x2 encoder columns, actual8x8 Fredkin and actual ancilla matrices; no padding/failed-return projection",
            "native_field_hardware_tape_or_protection_history_executed": False,
            "physical_time_or_energy_dictionary_derived": False},
        "premises_and_scope": {"chosen_graph": "three open sites per chain, unit coordinate spacing; three disjoint native-carrier chains",
            "chosen_native_control": "matching-support orbit average of H2 plus collective whole-carrier exchange",
            "exchange_coefficient_rule": "J=sum_(a<b)(W**2)_ab=125/1024, fixed before comparisons",
            "six_level_spectator_reinterpreted_as_three_independent_species": False,
            "vacuum_stability_domain": "chosen conserved spectator sector of each chain, not unrestricted native-code ground selection",
            "all_nine_occupation_modes_independent": True,
            "new_empirical_mass_or_calibration_packet_inputs": [],
            "original_CE_cycle_selects_this_architecture_or_coupling": False,
            "derived_unfitted_positive_gaps": list(map(str, p["positive_occupation_gaps"])),
            "gap_ratio_in_ascending_order": [1, 3, 12],
            "measured_SM_hierarchy_or_relativistic_particle_poles_recovered": False,
            "three_dimensional_local_CAR_or_chiral_continuum_derived": False,
            "nearest_result_and_addition": "parameter_selection: H2 has a leaky six-level compression; here a declared finite native average closes it exactly and collective exchange supplies an interacting charge/current field on independent carriers"}}
    return report, arrays


def source_bindings():
    paths = set(selection.source_bindings()) | {"model/common_source/"+name for name in (
        "repair_field_attachment.py", "test_repair_field_attachment.py", "chiral_candidate/gaussian_source.py",
        "mass_response/history.py", "mass_response/model.py", "halfline.py")}
    paths.update({"model/native_modular_reference/phase_inference.py",
        "model/scalar_observations/experimental_readout/stats.py", "model/requirements.txt"})
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths) if not name.endswith(".md")}


def compare_execution(report, arrays, fresh, fresh_arrays):
    retained_metrics = verify_retained_arrays(arrays)
    fresh_metrics = verify_retained_arrays(fresh_arrays)
    _validate_source_diagnostics(report, retained_metrics, "retained")
    _validate_source_diagnostics(fresh, fresh_metrics, "fresh")
    comparison = deepcopy(fresh)
    # These norms describe their own original binary source arrays. A fresh
    # platform may produce a different rounding residual while satisfying
    # exactly the same operator gate. Recompute the archived norms from the
    # archived arrays; retain strict comparison for every other report field.
    for section, metric in retained_metrics.items():
        comparison[section]["operator_error_diagnostic"] = metric
    selection._compare(comparison, report)
    if set(arrays) != set(fresh_arrays):
        raise ValueError("complete retained field artifact keys differ")
    for name, expected in fresh_arrays.items():
        value = arrays[name]
        if type(value) is not np.ndarray or value.dtype != expected.dtype or value.shape != expected.shape:
            raise ValueError("retained field artifact schema differs: "+name)
        if value.dtype.kind in "iu":
            if not np.array_equal(value, expected):
                raise ValueError("retained complete field records differ: "+name)
        elif not np.isfinite(value).all() or not np.allclose(value, expected, atol=REPLAY_TOLERANCE, rtol=0):
            raise ValueError("complete complex field/source/reference outputs differ: "+name)


def verify_retained_arrays(arrays):
    """Independently gate stored native amplitudes, separate from portability."""
    ideal = code.ideal_code_matrix()
    energy = np.asarray(list(map(float, _parameters()[0])))
    charge = [0]*6
    for label in _parameters()[2]:
        charge[label] = 1
    targets = {"native_encoded_six_columns": ideal,
        "native_orbit_averaged_outputs": ideal*np.exp(-1j*float(PROBE_TIME)*energy)[None, :],
        "native_charge_pointer_outputs": np.stack([ideal*np.array([int(n == record) for n in charge])[None, :] for record in (0, 1)]),
        "collective_fredkin_macro": _fredkin_target(),
        "collective_logical_projection": expm(-1j*float(logical_parameters()["exchange_strength"]*PROBE_TIME)*np.asarray(s.eye(4)-carrier_swap(0, 1, 2), complex))}
    errors = {}
    for name, target in targets.items():
        value = arrays.get(name)
        if (type(value) is not np.ndarray or value.dtype != np.dtype("complex128") or value.shape != target.shape
                or not np.isfinite(value).all()):
            raise ValueError("retained native target schema differs: "+name)
        difference = value-target
        error = (max(float(np.linalg.norm(branch, 2)) for branch in difference)
                 if difference.ndim == 3 else float(np.linalg.norm(difference, 2)))
        if error > SOURCE_GATE:
            raise ValueError("retained native amplitudes exceed source operator gate: "+name)
        errors[name] = error
    encoded = arrays.get("collective_encoded_two_columns")
    labels = logical_parameters()["vacuum_labels"][0], logical_parameters()["occupied_labels"][0]
    if (type(encoded) is not np.ndarray or not np.array_equal(encoded, arrays["native_encoded_six_columns"][:, labels])):
        raise ValueError("retained collective encoder factor lost its exact in-archive ancestry")
    macro, H, phase, scalar = (arrays.get(name) for name in ("collective_fredkin_macro", "collective_ancilla_H", "collective_ancilla_Z", "collective_scalar_I"))
    if any(type(value) is not np.ndarray or value.dtype != np.dtype("complex128") for value in (macro, H, phase, scalar)):
        raise ValueError("complete collective source macro arrays required")
    _collective_local_gates(macro, H, phase, scalar)
    incoming, returned, failed = _collective_outputs(encoded, macro, H, phase, scalar)
    theta = float(PROBE_TIME*logical_parameters()["exchange_strength"])
    swapped = incoming.reshape(512, 512, 4).transpose(1, 0, 2).reshape(1 << 18, 4)
    target = np.exp(-1j*theta)*(math.cos(theta)*incoming+1j*math.sin(theta)*swapped)
    if max(np.linalg.norm(returned-target, 2), np.linalg.norm(failed, 2)) > SOURCE_GATE:
        raise ValueError("retained collective factor output exceeds full-column source gate")
    if np.linalg.norm(incoming.conj().T@returned-arrays["collective_logical_projection"], 2) > REPLAY_TOLERANCE:
        raise ValueError("retained collective projection differs from its complete factor outputs")
    q = np.zeros((8, 2), complex); q[0, 0] = 1; q[[1, 2, 4], 1] = 1/math.sqrt(3)
    number = np.asarray(occupation(), complex)
    for name in ("field_source_unitaries", "uniform_gap_source_unitaries"):
        if (type(arrays.get(name)) is not np.ndarray or arrays[name].shape != (3, 8, 8)
                or arrays[name].dtype != np.dtype("complex128") or not np.isfinite(arrays[name]).all()):
            raise ValueError("complete all-Fock field source schema differs: "+name)
        for sector, U in enumerate(arrays[name]):
            if max(np.linalg.norm(U.conj().T@U-np.eye(8), 2), np.linalg.norm(U@number-number@U, 2)) > SOURCE_GATE:
                raise ValueError("retained all-Fock phase/charge source operator gate failed")
            if name == "uniform_gap_source_unitaries":
                p = logical_parameters()
                phases = np.exp(-1j*float(GAP_TIME)*np.array([float(3*p["onsite_baselines"][sector]),
                    float(3*p["onsite_baselines"][sector]+p["positive_occupation_gaps"][sector])]))
                if np.linalg.norm(U@q-q*phases[None, :], 2) > SOURCE_GATE:
                    raise ValueError("retained uniform charged phase does not meet its source gate")
    return {"native_onsite_source": errors["native_orbit_averaged_outputs"],
            "native_charge_source": errors["native_charge_pointer_outputs"]}


def _validate_source_diagnostics(report, metrics, origin):
    """Bind each rounding residual to its own arrays, not another backend."""
    for section, expected in metrics.items():
        value = report[section]["operator_error_diagnostic"]
        if (type(value) is not float or not math.isfinite(value) or not 0 <= value <= SOURCE_GATE
                or abs(value-expected) > RETAINED_DIAGNOSTIC_TOLERANCE):
            raise ValueError(origin+" "+section+" diagnostic differs from its original source arrays")


def _finite_json(value):
    if type(value) is dict:
        for item in value.values():
            _finite_json(item)
    elif type(value) is list:
        for item in value:
            _finite_json(item)
    elif type(value) is float and not math.isfinite(value):
        raise ValueError("nonfinite native-field JSON is inadmissible")


def check(path):
    path = Path(path)
    retained = strict_json.load(path)
    _finite_json(retained)
    pins = source_bindings()
    if type(retained) is not dict or retained.get("source_bindings") != pins:
        raise ValueError("native field closed source bindings differ")
    artifact = path.with_suffix(".npz")
    if retained.get("artifact") != {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError("native field complete artifact custody differs")
    with np.load(artifact, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError("duplicate native field artifact array names")
        arrays = {name: archive[name] for name in archive.files}
    verify_retained_arrays(arrays)
    fresh, fresh_arrays = execute(retained_counts=arrays.get("uniform_gap_complete_simulated_counts"))
    scientific = {key: value for key, value in retained.items() if key not in ("source_bindings", "artifact")}
    compare_execution(scientific, arrays, fresh, fresh_arrays)
    if source_bindings() != pins:
        raise ValueError("native field sources changed during replay")
    return scientific


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--output", type=Path)
    choice.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print("native correlation/exchange finite CAR field complete replay: PASS")
        return
    artifact = args.output.with_suffix(".npz")
    if args.output.exists() or artifact.exists():
        parser.error("fresh JSON and NPZ output paths required")
    pins = source_bindings()
    report, arrays = execute()
    if source_bindings() != pins:
        raise RuntimeError("native field source changed during execution")
    report["source_bindings"] = pins
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with artifact.open("xb") as stream:
        np.savez_compressed(stream, **arrays)
    report["artifact"] = {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False)+"\n")
    print(json.dumps(report["uniform_charged_gap_records"]["conditional_record_gap_intervals_numeric"]))


if __name__ == "__main__":
    main()
