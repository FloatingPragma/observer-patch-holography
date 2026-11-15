"""Validate complete outcome partitions and archived selected poststates."""
from pathlib import Path
import copy
import json
import numpy as np

from particle_chain import HERE, CLOCK, ROOT, NEXT, sha, validate_pins


def validate_partition(result, contract, clock_contract):
    for name in ("source_id", "source_dimension", "clock_id", "source_H_array_sha256_le_float64",
                 "preparation_array_sha256_le_complex128"):
        if contract[name] != clock_contract[name]:
            raise ValueError("source/preparation/clock join changed: " + name)
    expected = [key for row in clock_contract["preparation_labels"]
                for key in row["complete_final_count_state_archive_keys"]]
    if [row["clock_archive_key"] for row in result["clock_inputs"]] != expected:
        raise ValueError("missing complete clock input")
    if len(result["first_detector_records"]) != 12 or len(result["continued_clock_and_second_detector_records"]) != 144:
        raise ValueError("missing detector branch")
    terms = result["terminal_records"]
    coverage = {(r["clock_input_index"], r["first_regional_charge"], r["second_count_class"], r["second_regional_charge"]) for r in terms}
    if coverage != {(i, a, n, b) for i in range(12) for a in range(4) for n in range(3) for b in range(4)} or len(terms) != 576:
        raise ValueError("incomplete or relabeled terminal law")
    failure_names = {"reset_timeout", "oracle", "pointer_first", "pointer_second"}
    all_detectors = result["first_detector_records"] + result["continued_clock_and_second_detector_records"]
    if any(set(r["stops"]) != failure_names for r in all_detectors):
        raise ValueError("missing apparatus failure")
    total = sum(r["probability"] for r in terms) + sum(sum(r["stops"].values()) for r in all_detectors)
    if abs(total - 1) > 2e-12:
        raise ValueError("unconditional probability partition changed")
    if contract["measurement_order"]["regional_Q_H_commute"]:
        raise ValueError("regional charge cannot be promoted to commuting energy joint")
    if any(contract["claim_boundary"][key] for key in (
            "stable_particle_from_probability_conservation", "asymptotic_lifetime_bound", "predicted_mass")):
        raise ValueError("finite measurement evidence cannot award physical stability or mass")
    return total


def run():
    result = json.loads((HERE / "RESULT.json").read_text())
    contract = json.loads((HERE / "PARTICLE_CONTRACT.json").read_text())
    clock = json.loads((CLOCK / "CLOCK_CONTRACT.json").read_text())
    validate_pins(result)
    validate_pins(clock)
    if sha(HERE / contract["result"]) != contract["result_sha256"] or sha(HERE / contract["state_archive"]) != contract["state_archive_sha256"]:
        raise ValueError("particle output custody changed")
    total = validate_partition(result, contract, clock)
    archive = np.load(HERE / result["state_archive"], allow_pickle=False)
    coefficient = archive["terminal_energy_coefficients"]
    if coefficient.shape != (576, 348, 2, 2):
        raise ValueError("energy archive shape changed")
    probabilities = np.trace(coefficient, axis1=2, axis2=3).real
    error = float(np.max(abs(probabilities.sum(axis=1) - [r["probability"] for r in result["terminal_records"]])))
    if error > 2e-12 or probabilities.min() < -2e-12:
        raise ValueError("energy archive lost probability")
    transitions = archive["regional_transition_probabilities"]
    for row in result["terminal_records"]:
        position = (row["clock_input_index"], row["first_regional_charge"],
                    row["second_count_class"], row["second_regional_charge"])
        if abs(transitions[position] - row["probability"]) > 1e-13:
            raise ValueError("regional/count record relabeled")
    selected = result["selected_positive_excitation"]
    gi, ti = selected["energy_group_index"], selected["terminal_index"]
    group = result["complete_energy_groups"][gi]
    ids = archive[group["position_indices_archive_key"]]
    basis = archive[group["basis_archive_key"]]
    d = group["degeneracy"]
    reconstructed = np.zeros((384, 384), complex)
    reconstructed[np.ix_(ids, ids)] = basis @ coefficient[ti, gi, :d, :d] @ basis.conj().T
    state = archive[selected["state_archive_key"]]
    state_error = float(np.max(abs(reconstructed - state)))
    if state_error > 1e-12 or abs(np.trace(state).real - selected["probability"]) > 1e-12:
        raise ValueError("selected state replaced by a label")
    if group["global_charge"] != 1 or group["excitation_energy"] <= 0:
        raise ValueError("selected branch is not positive charged excitation")
    source = np.load(NEXT / "C10/joint_source_states.npz", allow_pickle=False)
    eigen_error = float(np.linalg.norm(source["hamiltonian"] @ state - group["energy"] * state))
    charge_error = float(np.linalg.norm(((-source["boundary_B"] - 1)[:, None]) * state))
    if max(eigen_error, charge_error) > 2e-12:
        raise ValueError("selected state has wrong actual H or charge support")
    # Hash-consistent in-memory outcome changes must fail structural checks.
    mutations = {}
    for name, edit in {
        "drop_tail": lambda r, c: r["terminal_records"].pop(),
        "drop_apparatus_failure": lambda r, c: r["first_detector_records"][0]["stops"].pop("oracle"),
        "success_renormalize": lambda r, c: [x.update(probability=x["probability"] / r["unconditional"]["primary_conclusive_probability"]) for x in r["terminal_records"]],
        "wrong_clock": lambda r, c: c.update(clock_id="passive192CAR"),
        "false_simultaneous_regional_energy": lambda r, c: c["measurement_order"].update(regional_Q_H_commute=True),
        "probability_conservation_as_stability": lambda r, c: c["claim_boundary"].update(stable_particle_from_probability_conservation=True),
    }.items():
        r, c = copy.deepcopy(result), copy.deepcopy(contract)
        edit(r, c)
        try:
            validate_partition(r, c, clock)
        except ValueError:
            mutations[name] = True
        else:
            raise ValueError("mutation accepted: " + name)
    return {"schema": "oph.particle_chain.producer_verification.v1", "complete_probability": total,
            "energy_probability_max_error": error, "selected_poststate_reconstruction_error": state_error,
            "selected_actual_H_eigenstate_error": eigen_error, "selected_actual_charge_support_error": charge_error,
            "mutations_rejected": mutations,
            "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in (Path(__file__), HERE / "PARTICLE_CONTRACT.json", HERE / "RESULT.json")}}


if __name__ == "__main__":
    result = run()
    (HERE / "VERIFY.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
