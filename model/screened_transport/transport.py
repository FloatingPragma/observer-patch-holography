"""Current response and corrected channel of the canonical screened q=1 code.

The current form describes one native attempted repair. The logical channel
describes the distinct cycle with record-selected recovery and core reset.
Neither is a physical Hamiltonian, proper-time action, or momentum spectrum.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import sympy as sp

from tools import strict_json
from ..exchange_sector import sector
from ..observer import record_recovery as records
from ..observer.source_discrimination import circuit
from ..particle_screening import screening
from ..screened_recovery import protocol

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "transport_results.json"
PROFILE = (1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
BASE_SOURCES = (
    "model/particle_screening/screening.py",
    "model/particle_screening/screening_results.json",
    "model/screened_recovery/protocol.py",
    "model/screened_recovery/experiment.py",
    "model/screened_recovery/screened_recovery_results.json",
    "model/screened_recovery/campaign.json",
    "model/screened_recovery/README.md",
    "model/screened_recovery/test_protocol.py",
    "model/observer/record_recovery.py",
    "model/observer/source_discrimination/circuit.py",
    "model/clock_defect_energy/mode_energy_results.json",
    "model/clock_defect_energy/reference_balance_results.json",
    "Lean/Geometry/SourceActionMeasure.lean",
    "Lean/Screen/RegionalContinuity.lean",
    "Lean/Geometry/InternalEnergyInertia.lean",
    "Lean/Geometry/ProperTimeInternalAction.lean",
    "Lean/Dynamics/WardLimitManifest.lean",
)
PACKET_SOURCES = (
    "model/screened_transport/__init__.py",
    "model/screened_transport/transport.py",
    "model/screened_transport/test_transport.py",
    "model/screened_transport/README.md",
    "model/screened_transport/campaign.json",
)


def digest(path: str) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def campaign_binding() -> dict:
    contract = strict_json.load(HERE / "campaign.json")
    expected = {
        "schema": "oph.screened_transport.campaign.v1",
        "campaign_id": "screened-transport-20261028-v1",
        "predecessor_git_commit": "1cfa6106e23d8d1279083de4ac38eec77f16a1d4",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        "boundary_status": "FROZEN_DIAGNOSTIC_CONTRACT",
        "whole_program_frozen": False,
        "native_family_changed": False,
        "integrated_gate_accepted": False,
        "physical_comparison_count_added": 0,
    }
    if any(type(contract.get(k)) is not type(value) or contract[k] != value
           for k, value in expected.items()):
        raise ValueError("screened transport campaign scope changed")
    sources = contract.get("source_sha256")
    if (type(sources) is not dict or set(sources) != set(BASE_SOURCES)
            or any(type(sources[name]) is not str or sources[name] != digest(name)
                   for name in BASE_SOURCES)):
        raise ValueError("screened transport predecessor custody changed")
    return {"campaign_id": contract["campaign_id"],
            "contract_sha256": digest("model/screened_transport/campaign.json")}


def _column_state(pair: tuple[int, int]) -> tuple[int, ...]:
    return tuple(int(port in pair) for port in range(12))


def current_response() -> dict:
    """Compute the same-code current form from formula and complete records."""
    columns = screening.sharp_code(screening.CORE)
    if columns is None or columns != protocol.SUPPORTS or len(columns) != 2:
        raise ArithmeticError("screened particle/recovery code identities differ")
    seams = records.seams()
    if len(seams) != 30 or any(tuple(sorted(pair)) in screening.SEAMS
                                for col in columns for pair in col):
        raise ArithmeticError("screened code lost native nonedge support")
    occupations = []
    for col in columns:
        occupations.append(tuple(F(sum(port in pair for pair in col), len(col))
                                 for port in range(12)))
    if occupations[0] != occupations[1]:
        raise ArithmeticError("one-body occupation reveals logical value")
    w = occupations[0]
    expected_w = (F(1, 2), F(1, 2), F(1, 4), F(1, 4), F(1, 4),
                  F(), F(1, 4), F(), F(), F(), F(), F())
    if w != expected_w:
        raise ArithmeticError("canonical screened occupations changed")
    matrix = sp.zeros(12)
    for a, b in seams:
        incidence = sp.zeros(12, 1)
        incidence[a], incidence[b] = -1, 1
        matrix += sp.Rational(w[a]+w[b], 60) * incidence * incidence.T

    complete_columns = {}
    for law in ("C", "E"):
        for logical, col in enumerate(columns):
            direct = sp.zeros(12)
            checked = 0
            for pair in col:
                state = _column_state(pair)
                for seam in seams:
                    for coin in (0, 1):
                        output, record = records.recorded_step(law, state, seam, coin)
                        signed = records.record_current(law, record)
                        delta = sp.Matrix([after-before for before, after
                                           in zip(state, output, strict=True)])
                        expected = sp.zeros(12, 1)
                        expected[seam[0]], expected[seam[1]] = -signed, signed
                        if delta != expected:
                            raise ArithmeticError("native current lost finite continuity")
                        direct += delta * delta.T / (len(col)*len(seams)*2)
                        checked += 1
            if direct != matrix:
                raise ArithmeticError("native complete records disagree with code current form")
            complete_columns[f"{law}{logical}"] = checked
    if set(complete_columns.values()) != {240} or matrix.rank() != 11:
        raise ArithmeticError("same-code current response rank changed")
    if sp.trace(matrix) != sp.Rational(1, 3):
        raise ArithmeticError("same-code movement rate changed")
    lam = sp.Symbol("lambda")
    characteristic = sp.factor((240*matrix).charpoly(lam).as_expr())
    expected_characteristic = (lam*(lam-4)
        *(lam**3-30*lam**2+254*lam-528)
        *(lam**3-16*lam**2+80*lam-126)
        *(lam**4-30*lam**3+310*lam**2-1270*lam+1728))
    if sp.expand(characteristic-expected_characteristic) != 0:
        raise ArithmeticError("current-form characteristic polynomial changed")
    v = sp.Matrix(PROFILE)
    form = (v.T*matrix*v)[0]/2
    if form != sp.Rational(13, 120):
        raise ArithmeticError("screened code mode current form changed")
    mode = []
    for col in columns:
        values = tuple(sum(PROFILE[port] for port in pair) for pair in col)
        mean = F(sum(values), len(values))
        variance = sum(((F(value)-mean)**2 for value in values), F())/len(values)
        if mean:
            raise ArithmeticError("screened mode has a logical first moment")
        mode.append({"support_values": list(values), "mean": str(mean),
                     "variance": str(variance), "current_form_over_variance": str(F(form)/variance)})
    if [(item["variance"], item["current_form_over_variance"]) for item in mode] != [
            ("1/2", "13/60"), ("5/2", "13/300")]:
        raise ArithmeticError("one common wave ratio was incorrectly inferred")
    return {
        "code_core": list(screening.CORE),
        "logical_supports": [[list(pair) for pair in col] for col in columns],
        "logical_independent_occupations": [str(value) for value in w],
        "laws_and_logical_columns_complete_attempts": complete_columns,
        "law_C_E_second_moments_equal": True,
        "source_current_second_moment_matrix": [[str(F(value)) for value in row] for row in matrix.tolist()],
        "current_form_formula": "C=(1/60) sum_seams (w_a+w_b)(e_b-e_a)(e_b-e_a)^T; D(v)=v^T C v/2",
        "current_form_rank": matrix.rank(),
        "current_form_characteristic_polynomial_of_240C": str(characteristic),
        "nonzero_current_probability_per_uniform_native_attempt": str(F(sp.trace(matrix))/2),
        "mode_profile": list(PROFILE),
        "mode_current_form": str(F(form)),
        "mode_logical_columns": mode,
        "mode_variance_any_normalized_logical_input": "|alpha|^2/2+5|beta|^2/2",
        "prior_profile_has_logical_input_independent_wave_ratio": False,
        "stationary_uniform_source_mode_ratio_from_separate_process": "1/10",
        "physical_energy_or_inertial_curvature_inferred": False,
    }


def corrected_cycle() -> dict:
    """Complete flagged logical channel and one-attempt core excursions."""
    service = protocol.exact_service_certificate()
    if (service["exact_branch_identity"] !=
            "R_h K_h V = sqrt(p_h) V; zero and abort slots have p_h=0"
            or not service["arbitrary_unknown_reference_preserved"]
            or service["q1_probability_after_label_reset"] != "1"):
        raise ArithmeticError("source-corrected branch identity changed")
    retained = strict_json.load(ROOT / "model/screened_recovery/screened_recovery_results.json")
    if (retained["service"]["positive_slots_by_law"] != {"C": 90, "E": 118}
            or not retained["source_replay"]["all_thirty_seam_recovery_source_replay"]
            or retained["source_replay"]["physical_service_or_roundoff_certificate"] is not False):
        raise ArithmeticError("all positive source recoveries were not retained")
    result = {}
    for law in ("C", "E"):
        distributions = Counter()
        positive = 0
        for seam in protocol.SEAMS:
            total = F()
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    p = protocol.branch_weight(law, seam, coin, raw)
                    total += p
                    if not p:
                        continue
                    positive += 1
                    moved = protocol.core_after_branch(law, seam, coin, raw)
                    if moved is None:
                        raise ArithmeticError("positive complete flag has no core")
                    distributions[moved] += p / len(protocol.SEAMS)
            if total != 1:
                raise ArithmeticError("complete screened flags lost probability")
        expected = {protocol.CORE: F(14, 15),
                    (0, 5): F(1, 120), (0, 7): F(1, 120),
                    (0, 8): F(1, 120), (0, 9): F(1, 120),
                    (1, 5): F(1, 120), (1, 7): F(1, 120),
                    (1, 10): F(1, 120), (1, 11): F(1, 120)}
        if dict(distributions) != expected or positive != retained["service"]["positive_slots_by_law"][law]:
            raise ArithmeticError("complete flagged core excursion law changed")
        resource = retained["resources"]["laws"][law]
        expected_pulses = "37499681/60" if law == "C" else "36870437/60"
        if (resource["uniform_seam_expected_native_plus_recovery_pulses"] != expected_pulses
                or resource["active_core_reset_probability_uniform_seam"] != "1/15"):
            raise ArithmeticError("same-code uniform-seam recovery resource control changed")
        result[law] = {
            "positive_complete_raw_flags": positive,
            "pre_reset_core_distribution": {str(k): str(v) for k, v in sorted(distributions.items())},
            "pre_reset_core_change_probability": str(1-distributions[protocol.CORE]),
            "post_reset_core_distribution": {str(protocol.CORE): "1"},
            "ideal_corrected_logical_superoperator_eigenvalues": ["1"]*4,
            "ideal_corrected_logical_generator_eigenvalues_per_attempt": ["0"]*4,
            "logical_channel_preserves_any_untouched_reference": True,
            "same_code_uniform_seam_expected_source_pulses_from_recovery_receipt": expected_pulses,
            "retained_classical_bits_per_attempt_including_core_archive":
                resource["fixed_format_total_classical_bits_per_uniform_attempt"],
        }
    return {
        "laws": result,
        "effective_attempt_boundary_core_transition": "delta_(0,1) to (0,1)",
        "uniform_independent_attempt_core_excursions": "probability any reset in t attempts = 1-(14/15)^t; expected reset count = t/15",
        "time_parameter": "native attempted repair count, not a calibrated observer clock",
        "unreset_fixed_code_q1_probability": "29/30",
        "source_compiler_resource_counts_owned_by_screened_recovery": True,
        "source_execution_error_is_exactly_zero_claimed": False,
        "common_logical_phase_or_rest_gap_fixed_by_channel": False,
        "continuous_momentum_or_coherent_dispersion_constructed": False,
        "gravitational_stress_constructed": False,
    }


def report() -> dict:
    binding = campaign_binding()
    source_mode = strict_json.load(ROOT / "model/clock_defect_energy/mode_energy_results.json")
    if (source_mode["native_record_forms"]["normalized_dirichlet_rate"] != "1/10"
            or source_mode["native_record_forms"]["mode_profile"] != list(map(str, PROFILE))
            or source_mode["native_record_forms"]["physical_energy_identified"] is not False):
        raise ArithmeticError("separate uniform-source mode control changed")
    return {
        "schema": "oph.screened_transport.v1",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        **binding,
        "native_family_changed": False,
        "cp04_accepted": False,
        "physical_comparisons_added": 0,
        "current_response": current_response(),
        "corrected_cycle": corrected_cycle(),
        "source_sha256": {name: digest(name) for name in (*BASE_SOURCES, *PACKET_SOURCES)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    choice.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = report()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8", newline="\n")
    elif args.check:
        actual = strict_json.load(args.check)
        if json.dumps(actual, sort_keys=True, allow_nan=False) != json.dumps(result, sort_keys=True, allow_nan=False):
            raise ValueError("screened transport receipt differs from exact source replay")
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
