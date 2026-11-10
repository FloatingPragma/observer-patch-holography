"""Complete finite phase readout of the admitted charged source.

The large source is replayed by its residual verifier. This module certifies
the explicit apparatus and all-outcome bounds; it does not diagonalize H.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
from model.affine_gauss_full_spectrum import certificate as full

PARENT = "model/affine_gauss_full_spectrum/results.json"
PARENT_SHA = "0121f60ec22ddc56bf540dab5195892918efe8e80048f4331518224c588d0250"
FINITE_CLOCK_PROOF = "model/coherent_defect_wave/chiral_joint/finite_clock/PROOF.md"
M = 2**20
SCALE = 16
HALF_WINDOW = F(3, 100)
SUBSTEPS = 4096


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def floor(value):
    value = F(value)
    return value.numerator // value.denominator


def ceil(value):
    return -floor(-F(value))


def outward(value, upper=False, places=24):
    value = F(value)
    scale = 10**places
    number = ceil(value * scale) if upper else floor(value * scale)
    return str(F(number, scale))


def interval(lo, hi):
    return [outward(lo), outward(hi, True)]


def tail_bound(size, distance):
    """Circle-grid tail outside distance d, derived in PROOF.md."""
    distance = F(distance)
    if not isinstance(size, int) or size < 2 or not 0 < distance <= F(1, 2):
        raise ValueError("finite pointer size and circle distance in (0,1/2] required")
    return min(F(1), 1 / (2 * size**2 * distance**2) + 1 / (2 * size * distance))


def acceptance_indices(center, half_window=HALF_WINDOW, size=M):
    center, half_window = F(center), F(half_window)
    if not 0 < center - half_window < center + half_window < SCALE:
        raise ValueError("the declared nonwrapping positive energy window is required")
    return ceil(size * (center - half_window) / SCALE), floor(size * (center + half_window) / SCALE)


def phase_kernel(energy, size=M):
    """Numerical diagnostic only, all finite pointer outcomes, positive Fourier sign."""
    x = np.arange(size, dtype=float) / size - float(energy) / SCALE
    # Center on the circle. The sinc ratio extends through every removable zero.
    x -= np.rint(x)
    return np.exp(1j * np.pi * (size - 1) * x) * np.sinc(size * x) / np.sinc(x)


def rational_clock_identities():
    """Exact denominator and integral remainder identities, independently replayed."""
    x, z, u = sp.symbols("x z u", real=True)
    a, b = 1 - x**2 / 10, x / 2 - x**3 / 120
    denominator = 1 + x**2 / 20 + x**4 / 600 + x**6 / 14400
    assert sp.expand(a*a + b*b - denominator) == 0
    p = lambda v: 1 + v / 2 + v**2 / 10 + v**3 / 120
    # Integration by parts: I_n = exp(z)/z - n I_(n-1)/z.
    integrals = [(sp.exp(z) - 1) / z]
    for degree in range(1, 7):
        integrals.append(sp.exp(z) / z - degree * integrals[-1] / z)
    weight = sp.Poly(u**3 * (1-u)**3, u)
    integral = sum(weight.nth(j) * integrals[j] for j in range(7))
    assert sp.simplify(sp.exp(z)*p(-z) - p(z) + z**7*integral/720) == 0
    assert sp.integrate(u**3*(1-u)**3, (u, 0, 1)) == sp.Rational(1, 140)
    return {"unitarity_denominator": str(denominator),
            "single_step_norm_error": "abs(t B)^7 / 100800",
            "n_step_norm_error": "abs(t B)^7 / (100800 n^6)"}


def clock_error(size=M, substeps=SUBSTEPS):
    # Raw vacuum energy <45, raw charged-sector H norm <57. For t0=pi/8,
    # pi<22/7 gives t0*57<627/28 and t0*45<495/28.
    return (size - 1) * (F(627, 28)**7 + F(495, 28)**7) / (100800 * substeps**6)


def car_step(word, mode, creation, modes=13):
    if not 0 <= mode < modes or not 0 <= word < 2**modes:
        raise ValueError("invalid finite CAR word")
    occupied = bool(word & (1 << mode))
    if occupied == creation:
        return None
    return word ^ (1 << mode), (-1)**((word & ((1 << mode)-1)).bit_count())


def verify_full_car(modes=13):
    """All occupation words and every creator/annihilator pair, integer signs."""
    def compose(word, left, right, left_creation=False):
        first = car_step(word, right, True, modes)
        if first is None:
            return {}
        second = car_step(first[0], left, left_creation, modes)
        return {} if second is None else {second[0]: first[1]*second[1]}

    for word in range(2**modes):
        for i in range(modes):
            for j in range(modes):
                total = compose(word, i, j)
                first = car_step(word, i, False, modes)
                if first is not None:
                    second = car_step(first[0], j, True, modes)
                    if second is not None:
                        total[second[0]] = total.get(second[0], 0) + first[1]*second[1]
                total = {k: v for k, v in total.items() if v}
                assert total == ({word: 1} if i == j else {})
                creators = compose(word, i, j, True)
                for target, sign in compose(word, j, i, True).items():
                    creators[target] = creators.get(target, 0) + sign
                assert not any(creators.values())
    return {"modes": modes, "occupation_words": 2**modes,
            "checked_mixed_CAR_columns": 2**modes * modes**2,
            "checked_creator_anticommutator_columns": 2**modes * modes**2,
            "signs_exact_integers": True}


@lru_cache(None)
def source():
    if digest(ROOT / PARENT) != PARENT_SHA:
        raise ValueError("the complete charged source receipt changed")
    saved = json.loads((ROOT / PARENT).read_text())
    full.coarse.validate_pins(saved["input_sha256"])
    if saved != full.certificate():
        raise ValueError("complete actual-H source replay differs")
    return saved


def source_bounds(saved):
    center = F(saved["energy_center_exact_binary_value"])
    radius = F(saved["arithmetic_admission"]["maximum_admitted_actual_H_residual"])
    coarse = json.loads((ROOT / saved["parent_receipt"]).read_text())
    lower2 = F(coarse["spectral"]["low_charged_levels"][1]["actual_H_energy_interval"][0])
    flo, fhi = map(F, saved["complete_ground_preparation"]["conditional_ground_window_probability_interval"])
    assert 0 < flo < fhi < 1 and center + radius < lower2 < 12
    return center, radius, lower2, flo, fhi


def record_bounds(saved, finite_clock=True):
    center, radius, lower2, flo, fhi = source_bounds(saved)
    lo, hi = acceptance_indices(center)
    ground_distance = (HALF_WINDOW - radius) / SCALE
    # Both possible directions around the phase circle are controlled.
    excited_distance = min((lower2-center-HALF_WINDOW)/SCALE,
                           1-(F(12)-center+HALF_WINDOW)/SCALE)
    alpha, beta = tail_bound(M, ground_distance), tail_bound(M, excited_distance)
    gate = clock_error() if finite_clock else F(0)
    assert 0 < alpha+gate < 1 and 0 < beta+gate < 1
    good = flo * (1-alpha-gate)
    bad = (1-flo) * (beta+gate)
    accepted_hi = fhi + (1-fhi)*(beta+gate)
    infidelity = bad/(good+bad)
    return {
        "apparatus": "finite R33 source gates and finite pointer phase" if finite_clock else "exact controlled source exponential",
        "accepted_index_interval_inclusive": [lo, hi],
        "accepted_pointer_outcomes": hi-lo+1,
        "rejected_pointer_outcomes": M-hi+lo-1,
        "ground_phase_distance_to_rejection_lower": str(ground_distance),
        "excited_phase_distance_to_acceptance_lower": str(excited_distance),
        "ideal_ground_rejection_upper": outward(alpha, True),
        "ideal_excited_acceptance_upper": outward(beta, True),
        "complete_controlled_clock_operator_error_upper": outward(gate, True),
        "ground_rejection_including_clock_upper": outward(alpha+gate, True),
        "excited_acceptance_including_clock_upper": outward(beta+gate, True),
        "pump_success_probability": "1/4",
        "pump_failure_probability": "3/4",
        "conditional_detector_acceptance_interval": interval(good, accepted_hi),
        "unconditional_accepted_record_interval": interval(good/4, accepted_hi/4),
        "unconditional_rejected_record_interval": interval((1-accepted_hi)/4, (1-good)/4),
        "unconditional_excited_false_acceptance_upper": outward(bad/4, True),
        "accepted_source_ground_fidelity_lower": outward(1-infidelity),
        "accepted_source_ground_infidelity_upper": outward(infidelity, True),
        "exact_ground_projector_implemented": False,
        "exact_three_aggregate_records": ["3/4", "p_accept/4", "(1-p_accept)/4"],
    }


def calculate():
    saved = source()
    center, radius, lower2, flo, fhi = source_bounds(saved)
    vector = np.load(ROOT / "model/affine_gauss_full_spectrum/trial_state.npy", allow_pickle=False)
    vacuum_coefficients = list(map(lambda x: F.from_float(float(x)), vector[::5456]))
    normalization_squared = sum((x*x for x in vacuum_coefficients), F(0))
    assert len(vacuum_coefficients) == 13 and normalization_squared > 0
    # The square-root normalization is an exact algebraic coefficient service.
    assert sum((x*x/normalization_squared for x in vacuum_coefficients), F(0)) == 1
    mp.iv.dps = 80
    capacity = mp.iv.log(2**13*14*2**12*5456*14*M)
    capacity_interval = interval(full.coarse.endpoint(capacity), full.coarse.endpoint(capacity, True))
    return {
        "schema": "oph.charged_finite_phase_readout.v1",
        "parent_receipt": PARENT, "parent_receipt_sha256": PARENT_SHA,
        "source_replay": "full70928-state actual-H residual verifier; all parent pins; no eigensolver",
        "same_source_raw_H_neutral_initial_state_reference_and_clock": True,
        "source": {"physical_dimension": 44695552, "charged_sector_dimension": 70928,
            "source_modes": 13, "reference_charge_states": 14, "longitudinal_bits": 12,
            "photon_dimension": 5456, "total_photon_cutoff": 3,
            "relative_charged_H_norm_upper": "12", "raw_vacuum_energy_upper": "45",
            "raw_charged_H_norm_upper": "57", "reference_energy_in_charged_sector": "1",
            "charged_ground_energy_interval": saved["actual_charged_ground_energy_interval"],
            "second_charged_level_lower": str(lower2)},
        "pump": {"initial_state": "exact actual neutral vacuum with matched source registers",
            "mode_exact_numerators": list(map(str, vacuum_coefficients)),
            "mode_normalization_squared": str(normalization_squared),
            "unitary": "Pphys [sqrt(3)/2 I + (D_z-D_z†)/2] Pphys + (I-Pphys)",
            "all_number_outcomes": list(range(14)),
            "neutral_record_probability": "3/4", "charged_record_probability": "1/4",
            "other_number_record_probability": "0",
            "source_and_reference_number_match_on_every_nonzero_branch": True,
            "charged_success_ground_fidelity_interval": [str(flo), str(fhi)]},
        "phase_apparatus": {"pointer_dimension": M, "pointer_qubits": 20,
            "number_pointer_dimension": 14, "all_joint_pointer_records": 14*M,
            "record_index_domain": "(n,k) in {0,...,13} x {0,...,M-1}",
            "phase_pointer_on_number_failure": "k=0; all other k retained with zero weight",
            "unitary": "sum_j |j><j| tensor exp(-i j t0 Hraw) exp(+i j t0 Evac)",
            "base_time": "t0=2 pi/16 in the original source clock",
            "positive_Fourier_transform": "Fplus[k,j]=exp(+2 pi i k j/M)/sqrt(M)",
            "displayed_energy": "16 k/M above the source neutral vacuum",
            "accepted_energy_window": [str(center-HALF_WINDOW), str(center+HALF_WINDOW)],
            "source_vacuum_pointer_phase_explicit": True,
            "kraus": "K_k=(1/M) sum_(j=0)^(M-1) exp(2 pi i j k/M) U^j",
            "complete_effect_sum": "sum_k K_k† K_k=I on the whole charged sector",
            "exact_joint_law": "p(0,0)=3/4; p(1,k)=||K_k z||²/4; all other p(n,k)=0",
            "conditional_energy_range_for_alias_control": ["0", "12"],
            "higher_source_charge_sectors_discarded": False,
            "higher_and_nonphysical_sectors": "retained; extra phase interaction is identity off Pphys P_(N=1); autonomous source evolution continues",
            "common_source_evolution_retained": "all complete Kraus maps include exp(-i T Hraw); no branch is silently paused",
            "controlled_extra_interaction": "Pphys P1 (Hraw-Evac) tensor |1><1|, in addition to autonomous Hraw"},
        "exact_exponential_record_bounds": record_bounds(saved, False),
        "finite_rational_gate_record_bounds": record_bounds(saved, True),
        "finite_clock": {"exact_identities": rational_clock_identities(),
            "base_gate": "R33(t0 Hraw/n)^n R33(-t0 Evac/n)^n on Pphys P1; identity elsewhere",
            "substeps_per_raw_or_phase_factor": SUBSTEPS,
            "controlled_binary_power_factors": 20,
            "total_R33_factors_in_serial_binary_power_realization": 2*SUBSTEPS*(M-1),
            "largest_coherent_source_time": "(M-1) pi/8",
            "operator_error_upper_exact": str(clock_error()),
            "commutes_with_actual_H_all_Gauss_and_N": True,
            "eigenvectors_changed": False,
            "approximate_gate_eigenphase_identified_as_actual_energy": False,
            "declared_exact_coefficients_and_finite_linear_solve_service": True,
            "efficient_local_hardware_circuit_proved": False,
            "large_source_QPE_matrix_numerically_executed": False},
        "capacity": {"source_full_kinematic_dimension": 2**13*14*2**12*5456,
            "source_physical_dimension": 44695552,
            "retained_apparatus_dimension": 14*M,
            "source_plus_apparatus_log_dimension_formula": "ln(2^13*14*2^12*5456)+ln(14*2^20)",
            "source_plus_apparatus_log_dimension_nats_interval": capacity_interval,
            "log_dimension_nats_diagnostic": math.log(2**13*14*2**12*5456)+math.log(14*M),
            "extra_gate_synthesis_workspace_counted": False,
            "pointer_size_derived_from_memory_budget": False},
        "exact_CAR_replay": verify_full_car(),
        "scope": {"finite_energy_detector": True, "exact_ground_projection": False,
            "translation_invariant_mass_shell": False, "Lorentz_spin": False,
            "Standard_Model_species": None, "measured_inputs": [], "GeV_mass": None},
    }
