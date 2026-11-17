"""Exact action ambiguity and counted work for the active sharp-q=1 code.

Only the already committed uniform-seam C/E law is evaluated. Sector
reweightings are stationary mathematical controls, not admitted physical
energies. Pulse counts and currents are per attempted repair, not rest energy.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..clock_defect_energy import defect_history, reference_balance
from ..common_source import encoded_memory as memory, native_routing as routing
from ..exchange_sector import sector
from ..particle_screening import screening
from ..observer import record_recovery as records
from ..observer.source_discrimination import circuit
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "energy_bridge_results.json"
BASE_SOURCES = (
    "model/exchange_sector/sector.py",
    "model/exchange_sector/sector_results.json",
    "model/clock_defect_energy/reference_balance.py",
    "model/clock_defect_energy/reference_balance_results.json",
    "model/clock_defect_energy/defect_history.py",
    "model/clock_defect_energy/mode_energy_results.json",
    "model/common_source/encoded_memory.py",
    "model/common_source/native_routing.py",
    "model/observer/record_recovery.py",
    "model/observer/source_discrimination/circuit.py",
    "model/particle_screening/campaign.json",
    "model/particle_screening/README.md",
    "model/particle_screening/screening.py",
    "model/particle_screening/screening_results.json",
    "model/particle_screening/test_screening.py",
    "Lean/Geometry/SourceActionMeasure.lean",
    "Lean/Geometry/InternalEnergyInertia.lean",
    "Lean/Geometry/ProperTimeInternalAction.lean",
    "Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean",
    "Lean/Screen/RegionalContinuity.lean",
    "Lean/Dynamics/WardLimitManifest.lean",
)
PACKET_SOURCES = (
    "model/active_code_energy/__init__.py",
    "model/active_code_energy/energy_bridge.py",
    "model/active_code_energy/test_energy_bridge.py",
    "model/active_code_energy/README.md",
    "model/active_code_energy/campaign.json",
)
# Import the particle packet's canonical code geometry. Its all-seam KL
# certificate is source-bound below; only its current is tested here.
SCREENED_SUPPORTS = screening.sharp_code(sector.CORE)


def digest(path: str) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def campaign_binding() -> dict:
    campaign = strict_json.load(HERE / "campaign.json")
    required = {
        "schema": "oph.active_code_energy.campaign.v1",
        "campaign_id": "active-code-energy-20261028-v1",
        "predecessor_git_commit": "7ab50c9fa326b86344ba0dda6ca56804eb094fb7",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        "boundary_status": "FROZEN_DIAGNOSTIC_CONTRACT",
        "whole_program_frozen": False,
        "native_family_changed": False,
        "physical_comparison_count_added": 0,
        "integrated_gate_accepted": False,
    }
    if any(type(campaign.get(k)) is not type(v) or campaign[k] != v
           for k, v in required.items()):
        raise ValueError("active-code energy campaign scope changed")
    if (type(bound_pins(campaign.get("source_sha256"))) is not dict
            or set(bound_pins(campaign["source_sha256"])) != set(BASE_SOURCES)
            or any(campaign["source_sha256"][name] != digest(name)
                   for name in BASE_SOURCES)):
        raise ValueError("predecessor source-bound input changed")
    return {"campaign_id": campaign["campaign_id"],
            "contract_sha256": digest("model/active_code_energy/campaign.json")}


def _swap(pair: tuple[int, int], seam: tuple[int, int]) -> tuple[int, int]:
    a, b = seam
    return tuple(sorted(b if p == a else a if p == b else p for p in pair))


def sector_reweighting() -> dict:
    """Test two distinct positive stationary references on the *same* law."""
    pairs, edges, *_ = defect_history.geometry()
    cores = sector.CORES
    if len(pairs) != 66 or len(cores) != 66 or len(edges) != 300:
        raise ArithmeticError("source/core census changed")
    kernels = {law: reference_balance.kernel(law) for law in ("C", "E")}
    if kernels["C"] != kernels["E"]:
        raise ArithmeticError("C/E population kernels differ")
    kernel = kernels["C"]
    edge_set = {tuple(sorted(edge)) for edge in edges}
    for i in range(66):
        for j in range(66):
            expected = F(1, 60) if (i, j) in edge_set or (j, i) in edge_set else F()
            if i != j and kernel[i][j] != expected:
                raise ArithmeticError("complete source kernel lost an edge/current")
        if sum(kernel[i]) != 1:
            raise ArithmeticError("source kernel lost a zero-current hold")

    charges = {(i, c): len(set(pair) & set(core))
               for i, pair in enumerate(pairs) for c, core in enumerate(cores)}
    sizes = Counter(charges.values())
    if sizes != {0: 2970, 1: 1320, 2: 66}:
        raise ArithmeticError("co-moving charge sector size changed")
    original_mass = {q: F(sizes[q] * sector.WEIGHTS[q], 66) for q in range(3)}
    if original_mass != {0: F(2, 5), 1: F(5, 9), 2: F(2, 45)}:
        raise ArithmeticError("transported reference sector masses changed")
    rows = {}
    for factor in (F(1), F(2)):
        z = sum((original_mass[q] * (factor if q == 1 else 1)
                 for q in range(3)), F())
        weight = {(i, c): sector.WEIGHTS[q] * (factor if q == 1 else 1) / (66*z)
                  for (i, c), q in charges.items()}
        if sum(weight.values(), F()) != 1 or min(weight.values()) <= 0:
            raise ArithmeticError("sector-reweighted reference is not faithful")
        checks = 0
        gradient = F()
        for i, j in edges:
            removed = set(pairs[i]) - set(pairs[j])
            added = set(pairs[j]) - set(pairs[i])
            if len(removed) != 1 or len(added) != 1:
                raise ArithmeticError("native edge is not a one-particle current")
            seam = tuple(sorted((*removed, *added)))
            for c, core in enumerate(cores):
                moved = _swap(core, seam)
                d = cores.index(moved)
                if (charges[i, c] != charges[j, d]
                        or weight[i, c] * kernel[i][j]
                        != weight[j, d] * kernel[j][i]):
                    raise ArithmeticError("reweighted reference lost native balance")
                # For every U(q), each term has U(q')-U(q)=0. Check a
                # nonconstant representative rather than declaring a value.
                u = (F(0), F(1), F(3))
                gradient += weight[i, c] * kernel[i][j] * (u[charges[j, d]] - u[charges[i, c]])**2
                checks += 1
        if gradient:
            raise ArithmeticError("a q-only potential acquired native gradient cost")
        mass = {q: sum((w for key, w in weight.items() if charges[key] == q), F())
                for q in range(3)}
        if any(sum((weight[i, c] for c in range(66) if charges[i, c] == 1), F())
               / mass[1] != F(1, 66) for i in range(66)):
            raise ArithmeticError("stationary q1 source marginal ceased to be uniform")
        rows[str(factor)] = {
            "stationary_sector_masses": {str(q): str(mass[q]) for q in range(3)},
            "faithful_reversible_edge_core_checks": checks,
            "nonconstant_q_potential_native_gradient_form": str(gradient),
            "fixed_q1_conditional_state_weight": str(next(
                w / mass[1] for key, w in weight.items() if charges[key] == 1)),
            "conditional_q1_source_marginal_each_state": "1/66",
        }
    if rows["1"]["fixed_q1_conditional_state_weight"] != rows["2"]["fixed_q1_conditional_state_weight"]:
        raise ArithmeticError("sector reweighting changed conditional q1 records")
    if rows["2"]["stationary_sector_masses"] != {"0": "9/35", "1": "5/7", "2": "1/35"}:
        raise ArithmeticError("factor-two sector contrast changed")
    predecessor = reference_balance.co_moving()
    if (predecessor["joint_charge_native_gradient_form"] != "0"
            or not predecessor["charge_projectors_preserved_by_each_complete_record_block"]):
        raise ArithmeticError("complete raw record charge control changed")
    return {
        "joint_source_core_states": len(charges),
        "connected_charge_sector_sizes": {str(k): sizes[k] for k in range(3)},
        "C_E_complete_population_kernels_equal": True,
        "complete_raw_C_E_record_blocks_preserve_q": True,
        "reweight_q1_relative_to_q0_by_factor_two": rows,
        "relative_modular_gap_shift_q1_minus_q0": "-log(2)",
        "at_one_declared_common_beta_relative_energy_gap_shift": "-log(2)/beta",
        "one_physical_beta_or_energy_action_selected": False,
        "stationarity_scope": "reduced diagonal source/core population law under uniform seams; no stationary retained archive",
        "statement": "arbitrary positive sector reweightings are stationary; q-only action offsets have zero native gradient form",
    }


def active_code_work() -> dict:
    """Count currents for both codes; recovery resources for the old code."""
    particle = strict_json.load(ROOT / "model/particle_screening/screening_results.json")
    if (SCREENED_SUPPORTS is None
            or [[list(pair) for pair in col] for col in SCREENED_SUPPORTS]
            != particle["all_seam_code_geometry"]["edge_core_code_columns"]
            or not particle["all_seam_code_geometry"]["arbitrary_reference_exact_recovery_exists_per_retained_flag"]
            or particle["all_seam_code_geometry"]["source_encoder_and_recovery_compiled_for_this_new_code"] is not False):
        raise ArithmeticError("screened candidate differs from the source-bound particle packet")
    pairs, *_ = defect_history.geometry()
    kernel = reference_balance.kernel("C")
    activity = {pair: 1-kernel[i][i] for i, pair in enumerate(pairs)}
    supports = tuple(tuple(sorted(pair)) for logical in memory.LOGICAL_PAIRS for pair in logical)
    if supports != ((0, 2), (1, 3), (0, 3), (1, 2)):
        raise ArithmeticError("same q1 source code changed")
    code_activity = sum((activity[pair] for pair in supports), F()) / 4
    uniform = sum(activity.values(), F()) / 66
    q1_fixed = tuple(pair for pair in pairs if len(set(pair) & set(sector.CORE)) == 1)
    q1_fixed_background = sum((activity[pair] for pair in q1_fixed), F()) / len(q1_fixed)
    fixed_reference_background = sum((sector.WEIGHTS[sector.charge(x)] * activity[pair]
                                      for x, pair in zip(sector.BASIS, pairs, strict=True)), F())
    if (code_activity, uniform, q1_fixed_background, fixed_reference_background) != (
            F(1, 6), F(5, 33), F(23, 150), F(341, 2250)):
        raise ArithmeticError("same-code native current activities changed")
    screened_rows = []
    for logical in SCREENED_SUPPORTS:
        if any(len(set(pair) & set(sector.CORE)) != 1 for pair in logical):
            raise ArithmeticError("screened successor left sharp q1")
        rates = tuple(activity[pair] for pair in logical)
        if rates != (F(1, 6),) * 4:
            raise ArithmeticError("screened successor current rate changed")
        screened_rows.append({"pair_supports": [list(pair) for pair in logical],
                              "basis_current_rates": [str(rate) for rate in rates],
                              "logical_current_rate": str(sum(rates, F()) / 4)})
    screened_core_motion = {}
    native_seams = set(records.seams())
    for law in ("C", "E"):
        columns = []
        for logical in SCREENED_SUPPORTS:
            changes = exits = 0
            for pair in logical:
                state = tuple(int(port in pair) for port in range(12))
                for seam in records.seams():
                    for coin in (0, 1):
                        _, record = records.recorded_step(law, state, seam, coin)
                        moved = sector.moved_core(sector.CORE, law, record)
                        changes += moved != sector.CORE
                        exits += moved not in native_seams
            denominator = len(logical) * len(native_seams) * 2
            columns.append({"core_value_change_per_uniform_seam_attempt": str(F(changes, denominator)),
                            "transported_core_exits_native_seam_set_per_attempt": str(F(exits, denominator))})
        if columns != [{"core_value_change_per_uniform_seam_attempt": "1/15",
                        "transported_core_exits_native_seam_set_per_attempt": "1/30"}] * 2:
            raise ArithmeticError("screened successor native core motion changed")
        screened_core_motion[law] = columns

    laws = {}
    for law in ("C", "E"):
        weighted_recovery = F()
        reset_probability = F()
        abort_probability = F()
        positive_branches = []
        for coin in (0, 1):
            for raw in range(1 << circuit.RECORD_BITS[law]):
                p = memory.branch_weight(law, sector.SEAM, coin, raw)
                resources = memory.recovery_resources(law, sector.SEAM, coin, raw)
                weighted_recovery += p * resources["native_pulses"]
                abort_probability += p * resources["raw_abort"]
                if p:
                    value = circuit.record_values(law, coin)[raw]
                    moved = sector.moved_core(sector.CORE, law, (*sector.SEAM, coin, value))
                    if moved != sector.CORE:
                        reset_probability += p
                    positive_branches.append({"coin": coin, "raw": raw,
                                              "probability": str(p),
                                              "recovery_native_pulses": resources["native_pulses"]})
        repair = routing.routed_repair_resources(law, sector.SEAM)["native_pulses"]
        if sum((F(b["probability"]) for b in positive_branches), F()) != 1:
            raise ArithmeticError("code branch probability not normalized")
        if abort_probability or reset_probability != F(1, 4):
            raise ArithmeticError("code abort or core reset probability changed")
        laws[law] = {
            "positive_raw_branches": positive_branches,
            "raw_abort_probability_on_code": str(abort_probability),
            "mean_recovery_native_pulses_per_attempt": str(weighted_recovery),
            "fixed_routed_repair_native_pulses_per_attempt": repair,
            "mean_repair_plus_recovery_native_pulses_per_attempt": str(repair+weighted_recovery),
            "initial_encoder_native_pulses_per_preparation": memory.encoding_resources()["native_pulses"],
            "core_reset_value_change_probability": str(reset_probability),
            "archive_bits_per_nonabort_reset_entry": 15,
            "archive_entries_per_positive_code_attempt": 1,
        }
    if (laws["C"]["mean_recovery_native_pulses_per_attempt"] != "52670"
            or laws["E"]["mean_recovery_native_pulses_per_attempt"] != "89999/2"):
        raise ArithmeticError("source-executed recovery work changed")
    return {
        "code_core": list(sector.CORE),
        "source_executed_repair_seam": list(sector.SEAM),
        "prior_source_executed_code_nonzero_current_per_uniform_seam_attempt": str(code_activity),
        "screened_successor_code": {
            "particle_campaign_id": particle["campaign_id"],
            "logical_columns": screened_rows,
            "any_logical_input_nonzero_current_per_uniform_seam_attempt": "1/6",
            "basis_supports_all_have_same_degree": 10,
            "native_core_motion_by_law_and_logical_column": screened_core_motion,
            "all_seam_recovery_exists_algebraically_in_particle_packet": True,
            "source_pulse_recovery_executed": False,
        },
        "uniform_source_or_stationary_q1_background": str(uniform),
        "excess_over_stationary_q1_background_per_attempt": str(code_activity-uniform),
        "fixed_core_q1_conditional_background_not_stationary": str(q1_fixed_background),
        "excess_over_fixed_core_q1_conditional_background_per_attempt": str(code_activity-q1_fixed_background),
        "fixed_core_unconditional_reference_background_not_stationary": str(fixed_reference_background),
        "excess_over_fixed_core_unconditional_reference_per_attempt": str(code_activity-fixed_reference_background),
        "mean_native_current_excess_times_declared_cost_epsilon": "epsilon/66 per attempted uniform-seam repair",
        "power_with_supplied_attempt_interval_delta_tau": "epsilon/(66*delta_tau)",
        "uniform_seam_activity_and_one_seam_executed_recovery_form_one_protocol": False,
        "prior_source_executed_code_recovery_laws": laws,
        "physical_cost_per_pulse_or_current_supplied": False,
        "kinetic_curvature_for_same_code_measured": False,
        "rest_energy_or_gravitational_stress_inferred": False,
    }


def report() -> dict:
    campaign = campaign_binding()
    other = strict_json.load(ROOT / "model/clock_defect_energy/mode_energy_results.json")
    if (other["response"]["rest_gap"] != "E_n(0)=sqrt(Q)*theta=a"
            or other["response"]["ideal_independent_curvature"] != "1/a"
            or other["scope"]["localized_defect_identified"] is not False):
        raise ArithmeticError("independent product-mode control changed")
    return {
        "schema": "oph.active_code_energy.v1",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        **campaign,
        "predecessor_git_commit": "7ab50c9fa326b86344ba0dda6ca56804eb094fb7",
        "native_family_changed": False,
        "cp04_accepted": False,
        "physical_comparisons_added": 0,
        "sector_action_ambiguity": sector_reweighting(),
        "active_code_work": active_code_work(),
        "independent_product_mode_control": {
            "ideal_rest_gap": other["response"]["rest_gap"],
            "ideal_kinetic_curvature": other["response"]["ideal_independent_curvature"],
            "same_active_q1_code": False,
            "same_source_or_momentum_dictionary": False,
            "physical_mass_transfer_allowed": False,
        },
        "source_sha256": bound_pins({name: digest(name) for name in (*BASE_SOURCES, *PACKET_SOURCES)}),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    choice.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = report()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif args.check:
        actual = strict_json.load(args.check)
        if json.dumps(actual, sort_keys=True, allow_nan=False) != json.dumps(result, sort_keys=True, allow_nan=False):
            raise ValueError("active-code energy receipt differs from source replay")
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
