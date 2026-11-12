"""The declared 432-state register and its retrospective lepton comparison.

Q_o=I-K^2+(-1)^o K, K=|p,1><m,0|+|m,0><p,1|;
U0=T (VW tensor VW tensor I3), W=Q R_orientation.
The coherent analyzer compresses U(theta) to
M=(sqrt(2)/3)[I/sqrt(6)+(exp(i theta)R+exp(-i theta)R^-1)/sqrt(12)].
The mass-read convention squares its three Fourier amplitudes. It is not
a theorem identifying interacting physical masses or a prospective test.
"""
from fractions import Fraction as F
import json

import mpmath as mp
import numpy as np


SOURCE_REFS = [
    "model/lepton_register_round/register.py",
    "model/lepton_register_round/kick.py",
    "model/lepton_register_round/certificate.py",
    "model/lepton_cycle/experimental_inputs.json",
    "model/lepton_cycle/compare.py",
    "model/lepton_mass_compiler/PROOF.md",
]


def ix(load, record, orientation):
    return 4 * load + 2 * record + orientation


def object_tick(oriented=True):
    """Return 2VW, derived from the six-state seam toggle K."""
    minus = np.array([1, 0, -1, 0, 0, 0], dtype=np.int64)
    plus = np.array([0, 1, 0, 1, 0, 0], dtype=np.int64)
    twice_k = np.outer(minus, plus) + np.outer(plus, minus)
    square = twice_k @ twice_k
    if np.any(square % 2):
        raise ArithmeticError("toggle square must have half-integer coefficients")
    event = np.zeros((12, 12), dtype=np.int64)
    for o in (0, 1):
        event[o::2, o::2] = 2 * np.eye(6, dtype=np.int64) - square // 2 + (1 - 2 * o) * twice_k
    tick = event[:, np.arange(12) ^ 1].copy()
    if oriented:
        for load in range(3):
            tick[ix(load, 1, 1)] *= -1
    return tick


def register_tick(oriented=True):
    """Return the literal integer matrix 4U0; retain all 432 coordinates."""
    obj = object_tick(oriented)
    pair = np.kron(obj, obj)
    out = np.zeros((432, 432), dtype=np.int64)
    for pair_row in range(144):
        c = pair_row // 12
        shift = (1 - 2 * (c % 2)) if (c // 2) % 2 else 0
        for register in range(3):
            out[3 * pair_row + (register + shift) % 3, register::3] = pair[pair_row]
    return out


def verify_tick(tick):
    if tick.shape != (432, 432) or not np.array_equal(tick.T @ tick, 16 * np.eye(432, dtype=np.int64)):
        raise ArithmeticError("complete tick is not unitary")


def source_and_analyzer():
    """Integer coefficients of 6|in> and 4sqrt(3)|Phi>."""
    tail, idle = np.zeros(12, dtype=np.int64), np.zeros(12, dtype=np.int64)
    tail[0:2] = 1
    idle[8:10] = 1
    source = (np.kron(tail, tail) + 2 * np.kron(idle, idle),
              np.kron(tail, idle) + np.kron(idle, tail))
    flags = [np.zeros(144, dtype=np.int64) for _ in range(3)]
    for oc in (0, 1):
        for le in (0, 1):
            for oe in (0, 1):
                flags[0][12 * ix(2, 0, oc) + ix(le, 0, oe)] = 1
    for lc in (0, 1):
        for oe in (0, 1):
            for oc in (0, 1):
                flags[1 + oc][12 * ix(lc, 1, oc) + ix(2, 0, oe)] = 1
    return source, flags


def compression_coefficients(tick):
    """Exact M_j=(A_j+sqrt(2) B_j)/(96sqrt(3)), j=hold,+,-."""
    source, flags = source_and_analyzer()
    a = np.zeros((3, 3, 3), dtype=np.int64)
    b = np.zeros_like(a)
    for col in range(3):
        basis = np.eye(3, dtype=np.int64)[:, col]
        out = [tick @ np.kron(vec, basis) for vec in source]
        for row in range(3):
            h, plus, minus = flags
            a[0, row, col], b[0, row, col] = 2 * (h @ out[1][row::3]), h @ out[0][row::3]
            for j, flag in enumerate((plus, minus), 1):
                a[j, row, col], b[j, row, col] = 2 * (flag @ out[0][row::3]), 2 * (flag @ out[1][row::3])
    return a, b


def verify_compression(tick):
    a, b = compression_coefficients(tick)
    shift = np.roll(np.eye(3, dtype=np.int64), 1, axis=0)
    expected_a, expected_b = np.zeros_like(a), np.zeros_like(b)
    expected_a[0] = 32 * np.eye(3, dtype=np.int64)
    expected_b[1], expected_b[2] = 16 * shift, 16 * shift.T
    if not np.array_equal(a, expected_a) or not np.array_equal(b, expected_b):
        raise ArithmeticError("coherent register compression differs from the balanced circulant")
    return a, b


def exact_read_counts(tick):
    source, _ = source_and_analyzer()
    x0 = np.array([1, 0, 0], dtype=np.int64)
    a, b = (tick @ np.kron(vec, x0) for vec in source)
    if np.any(a * b):
        raise ArithmeticError("read probabilities require the full radical cross term")
    weights = [F(int(x*x + 2*y*y), 576) for x, y in zip(a, b)]
    if sum(weights) != 1:
        raise ArithmeticError("full source weight was lost")
    counts = [F(0), F(0), F(0)]
    for j, weight in enumerate(weights):
        c, e = j // 36, (j // 3) % 12
        repair_plus = (c // 2) % 2 == 1 and c % 2 == 0
        agree = e // 4 < 2 and (e // 2) % 2 == 0
        for k, selected in enumerate((repair_plus, agree, repair_plus and agree)):
            if selected:
                counts[k] += weight
    return counts


def spectrum(theta=F(2, 9), balance_squared=F(1, 2)):
    """Independent high-precision three-character evaluation; no targets."""
    theta = mp.mpf(theta.numerator) / theta.denominator
    balance = 2 * mp.sqrt(mp.mpf(balance_squared.numerator) / balance_squared.denominator)
    return [1 + balance * mp.cos(theta + 2 * mp.pi * k / 3) for k in range(3)]


def comparison(theta=F(2, 9), balance_squared=F(1, 2)):
    with mp.workdps(70):
        lam = spectrum(theta, balance_squared)
        if min(lam) <= 0:
            raise ArithmeticError("Koide square-root identity requires the positive chamber")
        masses = [x*x for x in lam]
        mu, tau = masses[2]/masses[1], masses[0]/masses[1]
        electron = mp.mpf("0.51099895069")
        electron_u = mp.mpf("0.00000000016")
        mu_ref, mu_u = mp.mpf("206.7682827"), mp.mpf("0.0000046")
        tau_ref, tau_u = mp.mpf("1776.93"), mp.mpf("0.09")
        target_tau = tau_ref/electron
        a, b = tau_u/electron, tau_ref*electron_u/electron**2
        uncertainty_range = [abs(a-b), a+b]
        out = {
            "muon_electron_ratio": mu, "tau_electron_ratio": tau,
            "koide_ratio": sum(masses)/sum(lam)**2,
            "muon_relative_residual_ppm": (mu/mu_ref-1)*10**6,
            "muon_reference_only_standardized_discrepancy": (mu-mu_ref)/mu_u,
            "tau_relative_residual_ppm": (tau/target_tau-1)*10**6,
            "tau_reference_only_discrepancy_low": (tau-target_tau)/uncertainty_range[1],
            "tau_reference_only_discrepancy_high": (tau-target_tau)/uncertainty_range[0],
            "electron_calibrated_electron_GeV": electron/1000,
            "electron_calibrated_muon_GeV": mu*electron/1000,
            "electron_calibrated_tau_GeV": tau*electron/1000,
        }
        return {key: mp.nstr(value, 30) for key, value in out.items()}


def run():
    tick = register_tick()
    verify_tick(tick)
    a, b = verify_compression(tick)
    counts = exact_read_counts(tick)
    exclusive = counts[0] + counts[1] - 2 * counts[2]
    if exclusive != F(2, 9):
        raise ArithmeticError("exclusive read count changed")
    try:
        verify_compression(register_tick(oriented=False))
    except ArithmeticError:
        broken_orientation_rejected = True
    else:
        raise ArithmeticError("missing-orientation mutation was not detected")
    balance_squared = F(2*int(b[1, 1, 0])**2, int(a[0, 0, 0])**2)
    values = comparison(exclusive, balance_squared)
    return {
        "schema": "oph.openwave.leptons.v1", "computational_pass": True,
        "full_criterion_met": False, "criterion_rows": [2, 11], "source_refs": SOURCE_REFS,
        "checks": {"unitarity_exact_integer": True, "circulant_coefficients_exact": True,
                   "complete_source_probability_retained": True, "exclusive_read_count_exact": True,
                   "three_characters_of_declared_C3": True},
        "evidence": {"complete_tick_dimension": 432, "repair_plus_agree_joint": list(map(str, counts)),
                     "exclusive_count": str(exclusive), "balance_squared": str(balance_squared),
                     "coefficient_A": a.tolist(), "coefficient_B": b.tolist()},
        "results": values,
        "mutations": {"oriented_record_removed_rejected": broken_orientation_rejected,
                      "exact_muon_agreement_rejected": abs(float(values["muon_reference_only_standardized_discrepancy"])) > 5},
        "inputs": {"declared_phase": "2/9", "off_diagonal_to_diagonal_squared": "1/2",
                   "squared_eigenvalue_interpretation": "declared effect-rate to free-mass compiler",
                   "mode_assignment": {"electron": 1, "muon": 2, "tau": 0},
                   "electron_anchor_MeV": "0.51099895069", "electron_anchor_standard_uncertainty_MeV": "0.00000000016",
                   "muon_electron_reference": "206.7682827", "muon_electron_reference_standard_uncertainty": "0.0000046",
                   "tau_reference_MeV": "1776.93", "tau_reference_standard_uncertainty_MeV": "0.09",
                   "reference_editions": ["CODATA2022", "PDG2024"],
                   "empirical_exposure": "Particle masses and Koide relations were known during architecture development.",
                   "prospective_or_blind_prediction": False},
        "scope": {"physical_mass_poles_derived": False, "three_physical_mass_minima_proved": False,
                  "absolute_GeV_scale_derived": False, "muon_residual_cause_derived": False,
                  "physical_theory_uncertainty_available": False, "joint_covariance_or_chi_square_available": False,
                  "full_native_432_state_tick_rebuilt": True, "complete_phase_and_mass_cycle_reexecuted": False,
                  "conditional_source": "Two independent seam objects, declared preparation/analyzer/orientation and C3 register."},
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
