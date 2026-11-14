"""Complete ordered charge/energy/transport records on the frozen C10 clock."""
from pathlib import Path
from fractions import Fraction
import hashlib
import importlib.util
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
NEXT = ROOT / "openwave-tuning/next"
CLOCK = HERE.parent / "clock"
OP = NEXT / "ROW1_NATIVE/routes/quantum_operational"
sys.path[:0] = [str(CLOCK), str(OP), str(NEXT / "C10")]
from common_clock import CommonClock
from spectral_readout import GroupedEnergyInstrument
import operational_readout as op


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trace(state):
    return float(np.trace(state).real)


def validate_pins(document, root=ROOT):
    for name, expected in document["source_sha256"].items():
        if sha(root / name) != expected:
            raise ValueError("source pin changed: " + name)


def rectangular_apply(matrix, target, state, bits=5):
    """Apply a rectangular export branch; retain every spectator coordinate."""
    axes = [bits - 1 - target] + [i for i in range(bits) if i != bits - 1 - target] + [bits]
    arranged = state.reshape([2] * bits + [state.shape[1]]).transpose(axes)
    return (matrix @ arranged.reshape(2, -1)).reshape(-1, state.shape[1])


def export_families(oracle, pointer_run, arrays):
    """Explicit first-failure Kraus maps into flagged non-C10 output spaces."""
    path = OP.parent / "auditor_bridge/native_pointer_bridge.py"
    spec = importlib.util.spec_from_file_location("particle_pointer", path)
    provider = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(provider)
    records, _ = provider.actual_pointer()
    families = {name: [] for name in ("oracle", "pointer_first", "pointer_second")}
    for event in range(len(oracle.word)):
        for complement, k in enumerate(oracle.failure(event)):
            families["oracle"].append(({"event": event, "complement": complement}, k))
    for event in range(len(pointer_run.word)):
        for complement, k in enumerate(pointer_run.failure(event)):
            meta = {"event": event, "complement": complement}
            families["pointer_first"].append((meta, rectangular_apply(k, 3, oracle.success)))
            for native_record, (record, operators) in enumerate(records["C"]):
                for hidden, first in enumerate(operators):
                    state = op.compiler._apply(oracle.success, first, [3], 5)
                    families["pointer_second"].append((
                        {**meta, "first_native_record": native_record,
                         "first_record": record, "first_hidden_output": hidden},
                        rectangular_apply(k, 4, state)))
    rows = {}
    effects = {}
    for name, family in families.items():
        rows[name] = []
        effects[name] = np.zeros((8, 8), complex)
        for i, (meta, k) in enumerate(family):
            key = f"abort_{name}_{i}"
            arrays[key] = k
            effects[name] += k.conj().T @ k
            rows[name].append({**meta, "local_kraus_key": key,
                               "retained_local_output_dimension": len(k)})
    return rows, effects


class ParticleChain:
    def __init__(self):
        self.clock_contract = json.loads((CLOCK / "CLOCK_CONTRACT.json").read_text())
        self.operational = json.loads((OP / "OPERATIONAL_QUANTUM_REVIEW.json").read_text())
        validate_pins(self.clock_contract)
        validate_pins(self.operational)
        for key, field in (("receipt", "receipt_sha256"), ("state_archive", "state_archive_sha256")):
            if sha(CLOCK / self.clock_contract[key]) != self.clock_contract[field]:
                raise ValueError("clock output custody failed")
        if sha(OP / self.operational["archive"]) != self.operational["archive_sha256"]:
            raise ValueError("operational archive custody failed")
        self.clock = CommonClock()
        self.model = self.clock.model
        self.energy = GroupedEnergyInstrument(self.model)
        self.arrays = dict(self.energy.arrays)
        self.indices = op.encoded_indices(self.model.states)
        self.arrays["coherent_code_indices"] = self.indices
        decoded = op.decoded_occupations(
            np.array([f for _, f in self.model.states]), 3 - self.model.numbers)
        self.regional = decoded[:, (0, 1, 3)].sum(axis=1)
        self.arrays["regional_occupation"] = self.regional
        oracle, _, _, _ = op.compiled_oracle(3)
        instruments, errors, pointer_run, _ = op.native_pointer_instrument(oracle)
        self.maps = instruments["C"]["branches"]
        self.failures, self.failure_effects = export_families(oracle, pointer_run, self.arrays)
        self.failure_effect_error = max(float(np.max(abs(
            self.failure_effects[name] - instruments["C"]["failure_effects"][name])))
            for name in self.failure_effects)
        self.map_error = errors["C"]["coarse_instrument_matrix_unit_error"]
        for q, operators in self.maps.items():
            self.arrays[f"C_Q{q}_local_kraus"] = np.stack(operators)
        self.reset_exact = Fraction(self.operational["unconditional_apparatus_laws"]["C"]["reset_success_exact"])
        self.reset = float(self.reset_exact)
        self.reset_fail = float(1 - self.reset_exact)
        self.sector_offsets = np.cumsum([0] + [len(ids)**2 for ids, *_ in self.model.spectral.values()])
        self.arrays["number_block_offsets"] = self.sector_offsets
        self.block_states = []
        self.number_block_roundoff = 0.

    def pack(self, rho):
        self.clock.require_block_diagonal(rho)
        packed = np.concatenate([rho[np.ix_(ids, ids)].ravel() for ids, *_ in self.model.spectral.values()])
        self.block_states.append(packed)
        return len(self.block_states) - 1

    def reduced(self, rho):
        flat = self.indices.ravel()
        return np.einsum("arcr->ac", rho[np.ix_(flat, flat)].reshape(8, 48, 8, 48))

    def measure(self, rho):
        reduced = self.reduced(rho)
        stops = {"reset_timeout": self.reset_fail * trace(rho)}
        stops.update({name: self.reset * float(np.trace(effect @ reduced).real)
                      for name, effect in self.failure_effects.items()})
        states = [self.reset * op.apply_local_channel(rho, self.maps[q], self.indices) for q in range(4)]
        for state in states:
            off = state.copy()
            for ids, *_ in self.model.spectral.values():
                off[np.ix_(ids, ids)] = 0
            self.number_block_roundoff = max(self.number_block_roundoff, float(np.max(abs(off))))
        residual = abs(sum(stops.values()) + sum(map(trace, states)) - trace(rho))
        if residual > 2e-12:
            raise ValueError("detector loses complete probability")
        return states, stops, residual

    def energy_coefficients(self, rho):
        """Every grouped H/Q outcome, including zero and degenerate groups."""
        coefficients = np.zeros((len(self.energy.groups), 2, 2), complex)
        # Maximal degeneracy is checked, so padded zeros have no state meaning.
        transformed = {q: vectors.conj().T @ rho[np.ix_(ids, ids)] @ vectors
                       for q, (ids, _, vectors) in self.model.spectral.items()}
        cursor = {q: 0 for q in self.model.spectral}
        for i, group in enumerate(self.energy.groups):
            q, d = group["source_number"], group["degeneracy"]
            if d > 2:
                raise ValueError("changed spectral degeneracy")
            start = cursor[q]
            coefficients[i, :d, :d] = transformed[q][start:start+d, start:start+d]
            cursor[q] += d
        return coefficients

    def reconstruct_energy(self, coefficient, group_index):
        group = self.energy.groups[group_index]
        ids = self.arrays[group["position_indices_archive_key"]]
        basis = self.arrays[group["basis_archive_key"]]
        d = group["degeneracy"]
        out = np.zeros((384, 384), complex)
        out[np.ix_(ids, ids)] = basis @ coefficient[:d, :d] @ basis.conj().T
        return out


def run():
    chain = ParticleChain()
    clock, model = chain.clock, chain.model
    clock_archive = np.load(CLOCK / "common_clock_states.npz", allow_pickle=False)
    calibration = json.loads((NEXT / "C11/CHARGED_ENERGY_UNITS.json").read_text())
    validate_pins(calibration, NEXT.parent)
    preparation, reservoir = clock.preparation()
    input_error = float(np.max(abs(preparation - clock_archive["input_state"])))
    groups = [{**group, "global_charge": group["source_number"] - 3,
               "excitation_energy": group["energy"] - model.ground,
               "calibrated_excitation_GeV": (group["energy"] - model.ground) * calibration["calibration"]["GeV_per_source_energy"]}
              for group in chain.energy.groups]
    times = chain.clock_contract["numerical_member"]
    duration, rate = times["duration"], times["native_attempt_rate"]
    labels, first_rows, second_rows, terminals, coefficients = [], [], [], [], []
    transition = np.zeros((12, 4, 3, 4))
    failure_total = 0.
    complete_error = energy_error = 0.
    work_first = work_second = work_clock = 0.
    continued_sum = np.zeros((384, 384), complex)
    selected = None
    for reservoir_row in chain.clock_contract["preparation_labels"]:
        for first_count, key in enumerate(reservoir_row["complete_final_count_state_archive_keys"]):
            ci = len(labels)
            rho = clock_archive[key]
            labels.append({"clock_input_index": ci, "reservoir_initial": reservoir_row["reservoir_initial"],
                           "reservoir_final": reservoir_row["reservoir_final"], "first_count_class": first_count,
                           "clock_archive_key": key, "probability": trace(rho)})
            first, stops, error = chain.measure(rho)
            complete_error = max(complete_error, error)
            failure_total += sum(stops.values())
            first_rows.append({"clock_input_index": ci, "input_state": {"clock_archive_key": key}, "stops": stops,
                               "conclusive_probabilities": list(map(trace, first))})
            work_first += float(np.einsum("ij,ji->", clock.h, sum(first) - chain.reset * rho).real)
            for q1, state in enumerate(first):
                continued_counts = clock.count_split(state, duration, rate)
                work_clock += float(np.einsum("ij,ji->", clock.h, sum(continued_counts) - state).real)
                for second_count, continued in enumerate(continued_counts):
                    block_index = chain.pack(continued)
                    continued_sum += continued
                    second, stops, error = chain.measure(continued)
                    complete_error = max(complete_error, error)
                    failure_total += sum(stops.values())
                    second_rows.append({"clock_input_index": ci, "first_regional_charge": q1,
                                        "second_count_class": second_count, "number_block_state_index": block_index,
                                        "stops": stops, "conclusive_probabilities": list(map(trace, second))})
                    work_second += float(np.einsum("ij,ji->", clock.h, sum(second) - chain.reset * continued).real)
                    for q2, terminal in enumerate(second):
                        coefficient = chain.energy_coefficients(terminal)
                        probabilities = np.trace(coefficient, axis1=1, axis2=2).real
                        energy_error = max(energy_error, abs(float(probabilities.sum()) - trace(terminal)))
                        ti = len(terminals)
                        terminals.append({"terminal_index": ti, "clock_input_index": ci,
                                          "first_regional_charge": q1, "second_count_class": second_count,
                                          "second_regional_charge": q2, "number_block_state_index": block_index,
                                          "probability": trace(terminal)})
                        coefficients.append(coefficient)
                        transition[ci, q1, second_count, q2] = trace(terminal)
                        for gi, group in enumerate(groups):
                            if group["global_charge"] == 1 and group["excitation_energy"] > 1e-10:
                                candidate = (float(probabilities[gi]), ti, gi)
                                if selected is None or candidate[0] > selected[0]:
                                    selected = candidate
    chain.arrays["continued_number_blocks"] = np.array(chain.block_states)
    chain.arrays["terminal_energy_coefficients"] = np.array(coefficients)
    chain.arrays["regional_transition_probabilities"] = transition
    chain.arrays["continued_ensemble_state"] = continued_sum
    # Continue the most probable positive Q=+1 complete terminal outcome.
    probability, ti, gi = selected
    excitation = chain.reconstruct_energy(coefficients[ti][gi], gi)
    chain.arrays["selected_excitation_state"] = excitation
    witness = []
    for count, state in enumerate(clock.count_split(excitation, duration, rate)):
        key = f"selected_excitation_continued_count{count}"
        chain.arrays[key] = state
        outputs, stops, error = chain.measure(state)
        for q, output in enumerate(outputs):
            chain.arrays[f"{key}_Q{q}"] = output
        witness.append({"count_class": count, "pre_detector_state_key": key, "stops": stops,
                        "conclusive_regional_charge_probabilities": list(map(trace, outputs))})
    # Measurement-order witness on an actual clock output, with full poststates.
    actual = clock_archive["final_state"]
    regional = sum(op.coarse_branches(actual, chain.regional).values())
    ecoeff = chain.energy_coefficients(actual)
    energy_first = sum(chain.reconstruct_energy(ecoeff[i], i) for i in range(len(groups)))
    energy_then_regional = sum(op.coarse_branches(energy_first, chain.regional).values())
    rcoeff = chain.energy_coefficients(regional)
    regional_then_energy = sum(chain.reconstruct_energy(rcoeff[i], i) for i in range(len(groups)))
    order_distance = float(np.linalg.norm(energy_then_regional - regional_then_energy, ord="nuc") / 2)
    chain.arrays["energy_then_regional_poststate"] = energy_then_regional
    chain.arrays["regional_then_energy_poststate"] = regional_then_energy
    global_commutator = float(np.max(abs(clock.h * (model.numbers[None, :] - model.numbers[:, None]))))
    regional_commutator = float(np.linalg.norm(clock.h * (chain.regional[None, :] - chain.regional[:, None]), 2))
    # Actual charged excitation transport is not reduced to probability conservation.
    q_initial = terminals[ti]["second_regional_charge"]
    selected_after = sum(chain.arrays[row["pre_detector_state_key"]] for row in witness)
    conditional_state_change = float(np.linalg.norm(selected_after - excitation, ord="nuc") / (2 * probability))
    native_intervention = clock.semigroup(excitation, duration, rate) - clock.semigroup(excitation, duration, 0.)
    native_intervention_distance = float(np.linalg.norm(native_intervention, ord="nuc") / (2 * probability))
    primary_success = float(transition.sum())
    count_marginal = transition.sum(axis=(0, 1, 3))
    lost_tail = float(transition[:, :, 2, :].sum())
    changed_region = sum(float(transition[:, q1, :, q2].sum()) for q1 in range(4) for q2 in range(4) if q1 != q2)
    positive_gap = min(g["excitation_energy"] for g in groups if g["excitation_energy"] > 1e-9)
    charged_gap = min(g["excitation_energy"] for g in groups if g["global_charge"] == 1)
    charged_group_ids = [i for i, g in enumerate(groups) if g["global_charge"] == 1]
    charged_transition = np.zeros((4, 4))
    for row, coefficient in zip(terminals, coefficients):
        charged_transition[row["first_regional_charge"], row["second_regional_charge"]] += float(np.trace(coefficient[charged_group_ids], axis1=1, axis2=2).real.sum())
    charged_probability = float(charged_transition.sum())
    checks = {
        "all12_frozen_clock_states_consumed": len(labels) == 12 and input_error < 1e-12,
        "actual_C_pointer_maps_and_full_export_Kraus_families": chain.map_error < 1e-11 and chain.failure_effect_error < 1e-40,
        "all_first_second_apparatus_stops_retained": abs(primary_success + failure_total - 1) < 2e-12,
        "complete576_ordered_regional_and200448_energy_branches": len(terminals) == 576 and len(terminals) * len(groups) == 200448,
        "terminal_H_Q_probabilities_complete": energy_error < 2e-12,
        "positive_excitation_uses_actual_global_reference": probability > 0 and positive_gap > 0 and abs(charged_gap - calibration["calibration"]["source_anchor_excitation"]) < 1e-12,
        "number_block_compression_within_declared_roundoff": chain.number_block_roundoff < 1e-13,
        "global_charge_commutes_H_regional_charge_does_not": global_commutator < 1e-12 and regional_commutator > .1,
        "measurement_order_changes_actual_poststate": order_distance > .01,
        "positive_excitation_continued_with_complete_counts_and_stops": abs(sum(sum(x["stops"].values()) + sum(x["conclusive_regional_charge_probabilities"]) for x in witness) - probability) < 1e-12,
        "native_intervention_changes_excitation_poststate": native_intervention_distance > 1e-4,
        "noattempt_tail_and_reset_failure_nonzero": count_marginal[0] > .1 and lost_tail > .01 and chain.reset_fail > 0,
    }
    if not all(checks.values()):
        raise ValueError(checks)
    np.savez_compressed(HERE / "particle_states.npz", **chain.arrays)
    paths = [Path(__file__), CLOCK / "CLOCK_CONTRACT.json", CLOCK / "common_clock_states.npz",
             OP / "OPERATIONAL_QUANTUM_REVIEW.json", OP / "operational_states.npz",
             NEXT / "ROW1_NATIVE/routes/quantum_spatial/QUANTUM_SPATIAL.json",
             NEXT / "ROW1_NATIVE/routes/quantum_energy/JOINT_CHARGE_ENERGY.json",
             NEXT / "C11/CHARGED_ENERGY_UNITS.json", NEXT / "C10/spectral_readout.py"]
    pins = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result = {
        "schema": "oph.same_source.particle_chain_result.v1", "checks": checks,
        "source_sha256": pins, "clock_inputs": labels, "first_detector_records": first_rows,
        "continued_clock_and_second_detector_records": second_rows, "terminal_records": terminals,
        "complete_energy_groups": groups, "export_failure_kraus_families": chain.failures,
        "fine_native_readout_records": chain.operational["native_pointer_fine_records"]["C"],
        "reset_success_exact_per_detector": str(chain.reset_exact),
        "reset_timeout_exact_per_detector": str(1 - chain.reset_exact),
        "reset_contract": chain.operational["finite_native_reset_contract"],
        "unconditional": {"primary_conclusive_probability": primary_success, "all_apparatus_stop_probability": failure_total,
                          "regional_count_changed_probability": changed_region,
                          "regional_count_unchanged_probability": primary_success - changed_region,
                          "second_clock_count_probabilities_after_two_detector_successes": count_marginal.tolist(),
                          "detector1_success_weighted_source_energy_change": work_first,
                          "detector2_success_weighted_source_energy_change": work_second,
                          "continued_native_clock_source_energy_change": work_clock},
        "positive_charge_transport": {"global_charge": 1, "unconditional_probability": charged_probability,
                                      "regional_transition_matrix": charged_transition.tolist(),
                                      "regional_count_changed_probability": charged_probability - float(np.trace(charged_transition)),
                                      "regional_count_unchanged_probability": float(np.trace(charged_transition)),
                                      "interpretation": "Ordered regional counts in the positive global-charge sector; unchanged counts include compensating motion, with no tagged-particle or stable-localization inference"},
        "selected_positive_excitation": {"terminal_index": ti, "energy_group_index": gi,
                                         "probability": probability, "energy_group": groups[gi],
                                         "state_archive_key": "selected_excitation_state",
                                         "recorded_regional_charge_before_energy": q_initial,
                                         "selection_rule": "largest unconditional positive global-Q=+1 terminal energy probability; every other primary terminal state retained",
                                         "continued_count_and_regional_records": witness,
                                         "conditional_trace_distance_after_continuation": conditional_state_change,
                                         "native_rate_on_off_conditional_trace_distance": native_intervention_distance},
        "reference": {"energy": model.ground, "allsector_global_spectral_gap_native_units": positive_gap,
                      "charge_one_gap_above_global_reference_native_units": charged_gap,
                      "calibration": calibration["calibration"]},
        "mutations": {"discard_second_count_tail_loses_probability": lost_tail,
                      "success_only_renormalization_hides_stop_probability": failure_total,
                      "regional_energy_order_trace_distance": order_distance,
                      "global_Q_H_commutator_max": global_commutator,
                      "regional_Q_H_commutator_operator_norm": regional_commutator},
        "errors": {"preparation": input_error, "export_family_effect_reconstruction": chain.failure_effect_error,
                   "detector_complete_probability": complete_error, "terminal_energy_probability": energy_error,
                   "number_block_discarded_floating_roundoff_max_entry": chain.number_block_roundoff},
        "state_archive": "particle_states.npz", "state_archive_sha256": sha(HERE / "particle_states.npz"),
        "official_rows_changed": []}
    (HERE / "RESULT.json").write_text(json.dumps(result, indent=2) + "\n")
    contract = {
        "schema": "oph.source_state_clock_particle_contract.v1",
        "particle_instrument_id": "C10_ordered_regional_energy_transport_v1",
        "source_id": chain.clock_contract["source_id"], "source_dimension": 384,
        "source_H_array_sha256_le_float64": chain.clock_contract["source_H_array_sha256_le_float64"],
        "preparation_id": "C10_global_reference_complete_fair_reservoir_v1",
        "preparation_array_sha256_le_complex128": chain.clock_contract["preparation_array_sha256_le_complex128"],
        "reference_id": "actual_C10_allnumber_global_interacting_ground",
        "clock_id": chain.clock_contract["clock_id"], "clock_contract": str((CLOCK / "CLOCK_CONTRACT.json").relative_to(ROOT)),
        "operational_readout_receipt": str((OP / "OPERATIONAL_QUANTUM_REVIEW.json").relative_to(ROOT)),
        "source_sha256": pins,
        "ordered_program": ["complete four-branch reservoir preparation", "common-clock T=.37, lambda=.8, all N=0/1/>=2",
                            "actual operational regional-Q detector C, region 0/1/3", "same common-clock T=.37, lambda=.8, all N=0/1/>=2",
                            "same operational regional-Q detector", "complete grouped commuting global-Q and actual-H energy PVM"],
        "count_classes": ["N=0", "N=1", "N>=2 complete infinite-history complement"],
        "outcome_partition": {"first_detector_stops": "12 inputs x reset/oracle/first-pointer/second-pointer flags",
                              "second_detector_stops": "144 inputs x same four flags",
                              "conclusive_primary": "576 ordered regional/count records x348 grouped global-Q/energy records; zero bins retained",
                              "fine_flags": "RESULT export_failure_kraus_families and fine_native_readout_records; original reset policy retains every timeout history"},
        "state_reconstruction": {
            "clock_inputs": "clock common_clock_states.npz archive keys in RESULT.clock_inputs",
            "first_detector_conclusive": "reset * apply_local_channel(clock_input,C_Qq_local_kraus,coherent_code_indices)",
            "continued": "continued_number_blocks[row]; partition by number_block_offsets and Embed each block at spectral_indices_Nq",
            "second_detector_conclusive": "reset * apply_local_channel(continued,C_Qq_local_kraus,coherent_code_indices)",
            "terminal_energy": "Embed(spectral_indices, spectral_basis @ terminal_energy_coefficients[terminal,group,:d,:d] @ spectral_basis*)",
            "reset_timeout": "(1-reset) input C10 density tensor exact flagged auxiliary timeout history law; source paused",
            "export_failure": "reset * (K_local tensor I48) Encode(input) (K_local* tensor I48), separate flagged space for each archived rectangular abort K; no C10 energy assigned"},
        "measurement_order": {"global_Q_H_commute": True, "regional_Q_H_commute": False,
                              "meaning": "regional results are ordered historical records; terminal energy projection can disturb the preceding regional charge",
                              "energy_apparatus": "Declared grouped spectral Luders instrument, no native C/E spectral compiler claimed"},
        "resources_per_detector": chain.operational["clock_and_resources"],
        "coherent_transfer": chain.operational["apparatus_ancestry"]["coherent_source_transfer"],
        "calibration": calibration["calibration"],
        "claim_boundary": {"finite_measured_regional_transport": True, "positive_global_reference_excitation": True,
                           "actual_selected_excitation_poststate_continued": True,
                           "stable_particle_from_probability_conservation": False, "asymptotic_lifetime_bound": False,
                           "physical_particle_species_identified": False, "predicted_mass": False,
                           "physical_clock_calibrated": False, "continuum_limit_transferred": False},
        "result": "RESULT.json", "result_sha256": sha(HERE / "RESULT.json"),
        "state_archive": result["state_archive"], "state_archive_sha256": result["state_archive_sha256"],
        "official_rows_changed": []}
    (HERE / "PARTICLE_CONTRACT.json").write_text(json.dumps(contract, indent=2) + "\n")
    return result


if __name__ == "__main__":
    result = run()
    print(json.dumps({"checks": result["checks"], "unconditional": result["unconditional"],
                      "reference_gaps": {k: v for k, v in result["reference"].items() if "gap" in k},
                      "mutations": result["mutations"], "errors": result["errors"]}, indent=2))
